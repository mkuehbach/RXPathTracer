/*
	Copyright Markus Kühbach, Paul Wilhelm Hoffrogge, 2016

	RXPathTracer is an MPI-parallel program for simulating the growth stage during
	primary recrystallization. Its foundations are detailed in 
	
	  ''Efficient Recrystallization Microstructure Modeling by Utilizing Parallel Computation''
	  PhD thesis RWTH Aachen University
	  Successfully defended and accepted for publication August, 14th, 2017
	  Finally open source published in January, 2018
	  
	The source code was developed by Markus Kühbach (as part of his PhD) and Paul Hoffrogge 
	(as part of a bachelor thesis project) under supervision of Luis A. Barrales-Mora and 
	Günter Gottstein at the Institute of Physical Metallurgy and Metal Physics with RWTH Aachen University.  
	
	The authors gratefully acknowledge the financial support from the Deutsche Forschungsgemeinschaft
	(DFG) within the Reinhart Koselleck-Project (GO 335/44-1) and computing time grants kindly provided
	by RWTH Aachen University and the FZ Jülich within the scope of the JARAHPC project RWTH0157 and JARA0076.

	This file is part of RXPathTracer.

	RXPathTracer is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	TopologyTracer is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with SCORE.  If not, see <http://www.gnu.org/licenses/>.
*/


/* 
 * File:   main.cpp
 * Author: Markus
 *
 * Created on 27. Dezember 2013, 13:34
 */

#include "RXPathTracer_Kernel.h"

using namespace std;

#define INPUTFILE	1
#define SIMID		2
//#define MODELTYPE	3


int main(int argc, char** argv) {

//program startup: parameter <input.uds> <simid> <modeltype>
	double starttime = MPI_Wtime();
	if ( argc < 3 ) { cout << "ERROR::Too few input parameter." << endl; return 0; }
	//inputfile
	long simid = atol(argv[SIMID]);
	long modelid = 1; //##MK::atol(argv[MODELTYPE]);

	if ( simid < 1 ) { cout << "ERROR::Invalid simulation ID to identify files." << endl; return 0; }
	if ( modelid != SPHEREMESH_MODEL && modelid != VORONOI_MODEL ) { cout << "ERROR::1 (sphere) or 2 (voronoi) only model." << endl; return 0; }
	
	polyxxP poly = new polyxx;
	poly->simulationid = simid;
	poly->modelType = modelid;
	poly->realStartTimeOfComputing = starttime;
//init parallel environment MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &poly->nRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &poly->myRank);
	//report starting time
	poly->myComputeTimeLog.log("Starting Time [s]", MPI_Wtime());
	if (poly->myRank == MASTER){
	cout << "RXPathTracer_Kernel MPI startup successful." << endl;
	}
	poly->init_adjustPRNGs();
	poly->init_MPIDatatypes();

//Master reading input file, workers prepare to await parameter
	if ( poly->init_readparameter( argv[INPUTFILE] ) != true ) { MPI_Finalize(); delete poly; return 0; }

	//here all master/worker know all parameter, settle out which nuclei are processed by me, the poly->myRank, sytnh defms, apply nucleation models
	poly->init_workpartitioning_myNuclei();
	//poly->init_poissonpointprocess();
	poly->init_large_defms(); //##MK::modification to include grain boundary nucleation in large 3D cuboid aggregate from which a sample is taken
	poly->init_myBoxes();
	//poly->init_myNuclei();
	poly->init_myNuclei_BasedOnLargeDefMS();

	poly->init_disori_lowertriangle();


	//find out minmax values for individual timestepping
	poly->init_checkTimeSteppingOptimPotential();

	//individual timestepping
	if ( poly->modelType == SPHEREMESH_MODEL ) {
		poly->init_MyTimeSteps();
	}

	poly->init_myGrainSizeDistrBinning();

    //poly->init_disori_fast();
//--SIMULATION----------------------------------------------
	if(poly->modelType == SPHEREMESH_MODEL ) {

		poly->sim_ICO_HEDRONMESHING();
		poly->myComputeTimeLog.log("Sphere-Meshing fin[s]", MPI_Wtime());

cout <<"RANK: " << poly->myRank << " Starting sim_ICO_NBORPATHCONSTRUCTION." << endl;
		poly->sim_ICO_NBORPATHCONSTRUCTION();
		poly->myComputeTimeLog.log("Nbor-Path Construction fin[s]", MPI_Wtime());

cout <<"RANK: " << poly->myRank << " Starting sim_ICO_NUCPATHCONSTRUCTION." << endl;
		poly->sim_ICO_NUCPATHCONSTRUCTION();
		poly->myComputeTimeLog.log("Nuc-Path Construction fin[s]", MPI_Wtime());


cout <<"RANK: " << poly->myRank << " Starting sim_ICO_PATHGROWTH." << endl;
		poly->sim_ICO_PATHGROWTH();


cout << "RANK: "<< poly->myRank << " FINISH COMPUTING :-)" << endl;

		poly->myComputeTimeLog.log("Overall Computing fin[s]", MPI_Wtime());
		//wait until postprocessing into kinetics, texture, grain size distribution
		MPI_Barrier ( MPI_COMM_WORLD ); 

		MPI_Reduce(&(poly->mysumWallCollisions), &(poly->allsumWallCollisions), 1 , MPI_LONG, MPI_SUM, MASTER, MPI_COMM_WORLD);

		MPI_Reduce(&(poly->mysumIntegrationSteps), &(poly->allsumIntegrationSteps), 1, MPI_LONG, MPI_SUM, MASTER, MPI_COMM_WORLD);

		MPI_Allreduce ( &(poly->myMaxTime), &(poly->AllMaxTime), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

cout << "RANK: " << poly->myRank << " myMaxTime= " << poly->myMaxTime << " AllMaxTime = " << poly->AllMaxTime<<endl;

		//for each of the global timesteps find a corresponding local timestep
		poly->myComputeTimeLog.log("After 1st Barrier[s]", MPI_Wtime());
		poly->pp_mapMyTimeSteps();
		poly->myComputeTimeLog.log("PP:Timestep-Mapping fin[s]", MPI_Wtime());
		poly->pp_interpKineticsAndTextureData();
		poly->myComputeTimeLog.log("PP:Interpolation fin[s]", MPI_Wtime());
	} 
	else { //== VORONOI_MODEL
		cout << "VORONOI_MODEL currently not validated!" << endl; 
		/*
		poly->sim_VORO_PATHCONSTRUCTION();
		poly->sim_VORO_PATHGROWTH();
		poly->sim_VORO_PRECONDITIONING();
		poly->sim_VORO_CALCVOLUMES();
		*/
	}

	//Prepare GrainVolDist and VolList for writing
	poly->pp_GrainVolDistributionData();
	poly->myComputeTimeLog.log("PP:GrainVolDist fin[s]", MPI_Wtime());
	poly->pp_VolList();
	poly->myComputeTimeLog.log("PP:VolList fin[s]", MPI_Wtime());
//--COMMUNICATION AFTER SIMULATION----------------------------
	//Simulation ended and Data configured for Communication
	if(poly->modelType == SPHEREMESH_MODEL ) {
		//kinetics
                MPI_Barrier( MPI_COMM_WORLD );
                
                //PH::oripool necessarily the same for each polyxx on each process
                long nelements = ( 1 + poly->oripool.size() ) * poly->nRediscrEnsembleRealTimeIntervals; 
		MPI_Reduce( poly->myTimeStepData, poly->ensembleTimeStepData, nelements, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
		//binned grain size histogram
		MPI_Reduce( poly->myGrainVolHistoCounts, poly->ensembleGrainVolHistoCounts, poly->grvolhist.N_Bins, MPI_LONG, MPI_SUM, MASTER, MPI_COMM_WORLD);
		
                //collect in parallel the final volume from all grains
		MPI_Gatherv( &(poly->myNucRXFinalVolume[0]), poly->myIDs.size(), MPI_DOUBLE, &(poly->AllVolumes[0]), poly->WhoHasHowManyNuclei, poly->WhoHasHowManyNucleiCumulated, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
		
	}
	poly->myComputeTimeLog.log("Parallel End[s]", MPI_Wtime());

	if( poly->myRank == MASTER ) {
		poly->out_VolList(); //final volume list
		poly->out_GrainVolDistribution(); //final grain size distribution

		if ( poly->modelType == SPHEREMESH_MODEL ) {
			poly->out_Kinetics();
			poly->out_Texture(); //was switched off
			if(poly->allsumWallCollisions > 0){
				cout <<"WARNING: "<< poly->allsumWallCollisions << " Wall-Collisions occured! "<< endl;
			}
		}

		//DONE:
		cout << "MASTER AND ALL OTHERS ABOUT TO FINISH THE COMPUTATION." << endl;
	}
	poly->myComputeTimeLog.log("End[s]", MPI_Wtime());
	poly->report_RealTimeComputingData();
//-------------------------------------------------
	MPI_Finalize();

	delete poly;
	return 0;
}


/*
	//output logmatrix_rxgBoundingBox - THE BOUNDING BOXES
	stringstream log_rxgbbox_fname;
	ofstream log_rxgbbox_file;

	log_rxgbbox_fname << "START1DCA." << simulationid << ".rxgbbox.csv";
    cout << "File " << log_rxgbbox_fname.str().c_str() << " is opened now" << endl;
	log_rxgbbox_file.open ( log_rxgbbox_fname.str().c_str() );

	//header
	log_rxgbbox_file << "Step;tjmak;Xjmak;LEFT;RIGHT;REAR;FRONT;BOTTOM;TOP;CA=";
	for (long nuc = 0; nuc < nnucpool; nuc++) { log_rxgbbox_file << ";;;;;"; }
	log_rxgbbox_file << endl;

	//data block
	for (long lp = 0; lp < log_rxgbbox_cnt; lp++) {
		long lstep = logpoint_rxgbbox_stepwhenlogged[lp];

		log_rxgbbox_file << lstep << ";" << tjmak[lstep] << ";" << Xjmak[lstep] << ";";
		
		for ( long entry = 0; entry < (SIXDIRECTIONS*nnucpool); entry++ ) {
			log_rxgbbox_file << this->logmatrix_rxgBoundingBox[lp][entry] << ";";
		}
		log_rxgbbox_file << endl;
	}

	log_rxgbbox_file.flush();
	log_rxgbbox_file.close();
cout << "Bounding boxes extremal dimensions written into file." << endl;
	

	cout << "I am = " << poly->myRank << " aware of the fact that there are others... namely " << poly->nRanks << endl;

	delete poly;
	MPI_Finalize();
	return 0;

	char a = 0; //0x00;
	char b = 1; //0x01;
	long zahlen[2] = {11, 22};
	cout << "Zahlen a=" << zahlen[a] << endl;
	cout << "Zahlen b=" << zahlen[b] << endl;

	vector<long> a;
	a.push_back(1); a.push_back(2); a.push_back(3);
	vector<long> b;
	b.push_back(-1); b.push_back( -2 ); b.push_back( -3 );
	return 0;


	
	//I/O with all information in the MASTER
	
	/*
	//create 2D matrix hosting the MPI_SUMs from all nodes in the master
	double** rxgvol_mpisum = NULL;
	double** dgvol_mpisum = NULL;

	if (myRank == MASTER) {
		rxgvol_mpisum = new double*[log_texture_cnt];
		dgvol_mpisum = new double*[log_texture_cnt];
		
		for (long p = 0; p < log_texture_cnt; p++) { 
			rxgvol_mpisum[p] = NULL;
			rxgvol_mpisum[p] = new double[nstandardlagen + PLACE_FOR_RANDOM];
		
			dgvol_mpisum[p] = NULL;
			dgvol_mpisum[p] = new double[nstandardlagen + PLACE_FOR_RANDOM];
		
			for (long ic = 0; ic <= nstandardlagen; ic++) {
				rxgvol_mpisum[p][ic] = 0.0;
				dgvol_mpisum[p][ic] = 0.0;
			} //for all ic
		} //for all containers
	}
	*/
