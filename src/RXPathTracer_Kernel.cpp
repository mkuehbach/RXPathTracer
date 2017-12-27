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


#include <vector>

#include "RXPathTracer_Kernel.h"

inline bool SortByDistanceFromNucSite(const contactinfo &ci1, const contactinfo &ci2) { 
	return ci1.dist < ci2.dist; 
}

inline bool SortByContactTime(const contactinfo &ci1, const contactinfo &ci2) {
	return ci1.when < ci2.when;
}

inline bool SortDblAscending(const double &ai1 , const double &ai2) {
	return ai1 < ai2;
}


polyxx::polyxx()
{
	modelType = SPHEREMESH_MODEL; //spheremesh model
	meshRecursion = 2;

	bcast_ideal = NULL;
	bcast_dg = NULL;
	bcast_rxg = NULL;
	bcast_ori = NULL;
	bcast_Tt = NULL;
	bcast_zener = NULL;
	bcast_nucprintingschedule = NULL;
	AllComputingTimes = NULL;
	AllCollisions = NULL;
	AllIntegrationSteps = NULL;


	meanVol = 0.;
	meanCoordNum = 0.;
	
	standardlagen = NULL;
	nstandardlagen = 0;

	defgpool = NULL;
	ndefgpool = 0;
	rxgpool = NULL;
	nrxgpool = 0;
	LowTrigMobWeightHash = NULL;
	myCA_MobWeightHash = NULL;
	oripool_firstdisjunct_rxg = 0;

	//simulation meta data
	simulationid = 0;
	nCAEnsemble = 0;
	//nNeighborsPerCA = 0; //MK::deprectated
        ppdensity = 0.0;
        ppscaling = 0.0;
	//Vref = 0.0;
	minDist = 0.0;
	rhoMax = RHOMAX_DEFAULT;

	Tiso = 25 + TOFFSET;
	ZenerForceConstant = 0.0;

	defgmedian_rd = 1.0;
	defgmedian_nd = 1.0;
	defgmedian_td = 1.0;

	//##PH::initialize with some default values aluminum at 0K
	
	physConstants.init(27e3,0.,1.,2.86e-4,0.,0.,0.);
	
	LAGBm0 = 1000.0 * 1E12;
	LAGBHact = 1.60 * 1.602 * 1e-19;
	HAGBm0 = 3.0 * 1E12;
	HAGBHact = 1.20 * 1.602 * 1e-19;
	GSm0 = 3.0 * 1E12;
	GSHact = 1.20 * 1.602 * 1e-19;

        RHHAGBm0 = 300.0 * 1E12;
        RHHAGBHact = 2.45 * ECHARGE;
        RHLAGBHAGBCut = 0.0;
        RHLAGBHAGBTrans = 5.0;
        RHLAGBHAGBExponent = 9.0;

	//initialize with some value
	mLAGB = LAGBm0 * exp (-1 * LAGBHact / (kboltzman * (25.0 + TOFFSET) ));
	mHAGB = HAGBm0 * exp (-1 * HAGBHact / (kboltzman * (25.0 + TOFFSET) ));
	mGS = GSm0 * exp (-1 * GSHact / (kboltzman * (25.0 + TOFFSET) ));
	mRHHAGB = RHHAGBm0 * exp(-1 * RHHAGBHact / (kboltzman * (25.0 + TOFFSET) ));
        
	optTimeStepping.myPmax = 0.0;
	optTimeStepping.myPmin = 1.0;
	optTimeStepping.myRhomax = 0.0;

	grvolhist.VolDistBinSize = 1.0;
	grvolhist.VolDistMaxValue = 27.0;
	grvolhist.N_Bins = 27;
        
        largedefmsinfo.nx = 0;
        largedefmsinfo.ny = 0;
        largedefmsinfo.nz = 0;
        largedefmsinfo.nxy = 0;
        largedefmsinfo.nxyz = 0;
        largedefmsinfo.dimx = 1000.0; //in micron!
        largedefmsinfo.dimy = 10.0;
        largedefmsinfo.dimz = 100.0;
        //##MK::set as default here
        largedefmsinfo.gbnucdens2num = 2.0E-4;
        largedefmsinfo.gbnucdrho2dens = 3.5436e-15;
        largedefmsinfo.gbnucmaxscatter = 8.0 / 180.0 * _PI_; //##MK::exactly here 

		largedefmsnuc_treebins = NULL;

        kuhlmann_alpha = 0.0;
	kuhlmann_beta = 0.0;

	//caching variables
	zenerfac = 1.5 * 0.324; //Murr Meyers pure Alu GB energy and Zener-Smith 3/2

	temperature_isothermalonly = false;
	kuhlmann_consider = false;
	zener_consider = false;
        
	mobilitymodel_option = SEBALD_GOTTSTEIN;

	myRank = MASTER;
	nRanks = 1;

	WhoHasHowManyNuclei = NULL;
	WhoHasHowManyNucleiCumulated = NULL;

	myMemGuard = 0.0;
	mysumWallCollisions = 0;
	mysumIntegrationSteps = 0;
	allsumWallCollisions = 0;
	allsumIntegrationSteps = 0;

	mySuccess = false;
	rndlocal.init( -46356 );
	rxpathtracer_input = NULL;
}


polyxx::~polyxx()
{
	delete [] defgpool;
	delete [] rxgpool;
	delete [] LowTrigMobWeightHash;
	delete [] myCA_MobWeightHash;

	delete [] largedefmsnuc_treebins;

	for ( long lf = 0; lf < this->largedefmsnuc_tree.size(); lf++ ) {
		delete [] largedefmsnuc_tree[lf];
	}

	for ( long defms = 0; defms < myDefMS.size(); defms++) {
		delete [] myDefMS[defms];
	}


	for ( long mynuc = 0; mynuc < myIDs.size(); mynuc++) { //from the nucleus to the nbor
		for (long n_f = 0; n_f < numFaces; n_f++) {
			//delete [] this->myIcoPaths[mynuc][n_f].grainregions_outbound;
			//delete [] this->myIcoPaths[mynuc][n_f].regionends_outbound;
			//inbound not filled
			delete [] this->myIcoPaths[mynuc][n_f].grainregions;
		}
	}

	for ( long mynuc = 0; mynuc < myIDs.size(); mynuc++) { //from the nbor to the nucleus
		for (long nbor = 0; nbor < this->myNeighbors[mynuc].size(); nbor++) {

			//delete [] this->myNBorPaths[mynuc][nbor].grainregions_outbound;
			//delete [] this->myNBorPaths[mynuc][nbor].regionends_outbound;
			//delete [] this->myNBorPaths[mynuc][nbor].grainregions_inbound;
			//delete [] this->myNBorPaths[mynuc][nbor].regionends_inbound;
			delete [] this->myNBorPaths[mynuc][nbor].grainregions;
		}
	}

	for ( long i = 0; i < this->myTimeOripoolTextureMatrix.size(); i++ ) {
		free(this->myTimeOripoolTextureMatrix[i]);
	}

	/*for ( long i = 0; i < this->myLastTimeStepOripoolTexture.size(); i++ ) {
		delete [] this->myLastTimeStepOripoolTexture[i];
	}*/

	free(myTimeStepData);
	free(ensembleTimeStepData);
	free(myGrainVolHistoCounts);
	free(ensembleGrainVolHistoCounts);

	free(WhoHasHowManyNuclei);
	free(WhoHasHowManyNucleiCumulated);
}


void polyxx::specificallyDisorientednewOriFromReference( double *bunge, double sigma_rayl_dis2bunge, double *newOri ) //double sigma_rayl_dis2bunge_max, 
{
	//MK::the conventional newOrientationFromReference applies a random unit quaternion on the SO3 on the quaternion represenation
	//of orientation bunge, clearly and in accordance with MacKenzie one expects now that the distribution of disorientation of newOri to bunge
	//generated with such algorithm to follow similar a MacKenzie distribution, however for modeling nucleation spectra this is not desired
	//because the subgrains, at least in those alloys were they are of relevance for nucleation usually have a distinct orientation spread about
	//the average orientation of the grain

	//so it is in addition necessary to filter all possible generated newOris according to their expectation, otherwise
	//a nucleation model always generates practically highly disoriented nuclei

	//##MK::here we utilize an acceptance-rejection sampling
    
//principal idea here is to make an acceptance rejection sampling with particular distributed disorientation angles...
	double matrix[3];	
        matrix[0] = bunge[0];
	matrix[1] = bunge[1];
	matrix[2] = bunge[2];
        
        
        double qmatrix[4];
        euler2quaternion( matrix, qmatrix);
       
	
        double theta, prob, coinFlip;
        double q[4] = {0.0, 0.0, 0.0, 0.0};
        double qrndmisori[4];
        double rotated[4]; 
        double candidate[3];
	//acceptance rejection assuming the subgrains to follow a Rayleigh disorientation distribution 
       
        double th = MAXIMUM_MISORI_FCC;

 
 //cout << "...\t\t prior to while outer" << endl; 
        bool foundvalid = false;
 
	while ( foundvalid == false ) {
                //deviation is internally limited to MINIMUM_ANGULAR_SPREAD if called
                //////////////////////////////////randomMisorientationShoemake( MAXIMUM_MISORI_FCC, qrndmisori ) - ( double theta, double* qr  );
                //get a random quaternion on the SO3, Shoemake's algorithm...
                //find some
            
//cout << "...\t\t\t\t.. while inner" << endl;
                double q[4] = {0.0, 0.0, 0.0, 0.0};
                double qcrit = cos( 0.5 * MAXIMUM_MISORI_FCC);
                if ( th < MINIMUM_ANGULAR_SPREAD ) { //limit to assure the while loop to finish
                        qcrit = cos( 0.5 * MINIMUM_ANGULAR_SPREAD );
                }
                while( q[0] < qcrit ) {
                        double X0 = rndlargedefms.leEcuyer();
                        double X1 = rndlargedefms.leEcuyer();
                        double X2 = rndlargedefms.leEcuyer();

                        double r1 = sqrt(1-X0);
                        double r2 = sqrt(X0);
                        double theta1 = 2 * _PI_ * X1;
                        double theta2 = 2 * _PI_ * X2;

                        q[0] = r1 * sin(theta1); //w
                        q[1] = r1 * cos(theta1); //v
                        q[2] = r2 * sin(theta2); //x
                        q[3] = r2 * cos(theta2); //z

                        double qnorm = sqrt( SQR(q[0]) + SQR(q[1]) + SQR(q[2]) + SQR(q[3]) );

                        //normalize
                        q[0] = q[0] / qnorm;    
                        q[1] = q[1] / qnorm;    
                        q[2] = q[2] / qnorm;    
                        q[3] = q[3] / qnorm;
                        //definately this algorithm samples random on the SO(3) the resulting quaternion however is not necessarily a disorientation!
                }
//cout << "...\t\t\t\t.. while -->" << endl;

                //random quaternion known now
                qrndmisori[0] = q[0];
                qrndmisori[1] = q[1];
                qrndmisori[2] = q[2];
                qrndmisori[3] = q[3];             
      
                //apply this rotation to the matrix orientation
                multiplyQuaternions( qrndmisori , qmatrix, rotated );
         
                double tmpEuler[3];
                quaternion2Euler( rotated, tmpEuler );
                
                candidate[0] = tmpEuler[0];
                candidate[1] = tmpEuler[1];
                candidate[2] = tmpEuler[2];
                                
		theta = this->misorientationCubic( matrix[0], matrix[1], matrix[2], candidate[0], candidate[1], candidate[2] );
                
		prob = 1.0 - exp ( - 0.5 * SQR(theta/sigma_rayl_dis2bunge) );

		//how likely is theta to occur in Rayleigh distribution?, evaluate pdf Rayleigh with sigma_rayl_dis2bunge
		//prob = theta / SQR(sigma_rayl_dis2bunge) * exp ( -0.5 * SQR((theta/sigma_rayl_dis2bunge)) );

		//to avoid the implementation of calculating Lambert's W function that appears when trying to scale the f(sigma) to its maximum value
		//prob = prob / sigma_rayl_dis2bunge_max;

		//acceptance rejectance scheme based on disorientation angle that should be Rayleigh-distributed between nuclei and matrix grain
		coinFlip = rndlargedefms.leEcuyer();
//cout << "theta/prob/coinFlip = " << theta << ";" << prob << ";" << coinFlip << endl;

		if ( prob <= coinFlip )
			foundvalid = true;
	}
 //cout << "...\t\t done a while outer = " << rotated[0] << "--" << rotated[1] << endl;
        
        //get back to Euler angles
        quaternion2Euler( rotated, candidate );
        
	//report accepted orientation
	newOri[0] = candidate[0];
	newOri[1] = candidate[1];
	newOri[2] = candidate[2];
}


long polyxx::get_closest_standardlage( double * quat )
{
	double closest_disori = MAX_FCC;
	long closest_id = RANDOM;
	double q10 = quat[0];
	double q11 = quat[1];
	double q12 = quat[2];
	double q13 = quat[3];

	for (long cand = 0; cand < nstandardlagen; cand++) {
		double disori = misorientationCubicQxQ ( q10, q11, q12, q13,  standardlagen[cand].q0, standardlagen[cand].q1, standardlagen[cand].q2, standardlagen[cand].q3 );

		//scatter within the range and closer as to all other candidates?
		if ( disori <= standardlagen[cand].scatter && disori < closest_disori ) {
			closest_disori = disori;
			closest_id = cand;
		}
	} //for all possible standardlagen candidates

	return closest_id; //cannot be dereferenced!
}


long polyxx::check_disjunctness( double * bunge ) //MK::this function is the only one that is allowed to manage and to change the oripool!
{
	long closestid = -1;
	double disori = 2 * _PI_;
	double closestdisori = DISORI;
	double qbunge[4], qcand[4];

	euler2quaternion( bunge, qbunge );

	//check disorientation to all other components already present in the oripool, 
	//MK::this categorization into a finite set of disjoint orientations in quaternion space is different than what Hirsch Fortunier opted for, but is still random in quaternion space
	//MK::from an implementation point of view it limits the amount of different orientations and at the same time takes care that only relevant differences are handled which can also be verified by experiments!

	for (long cand = 0; cand < oripool.size(); cand++ ) {
		qcand[0] = oripool[cand].q0;
		qcand[1] = oripool[cand].q1;
		qcand[2] = oripool[cand].q2;
		qcand[3] = oripool[cand].q3;
		disori = misorientationCubicQxQ( qbunge[0], qbunge[1], qbunge[2], qbunge[3], qcand[0], qcand[1], qcand[2], qcand[3]);

		if (disori <= closestdisori) {
			closestdisori = disori;
			closestid = cand;
		}
//cout << "TEST\t\t" << cand << "\t\t" << disori << endl;
	}

	if (closestid == -1) { //write entry if closestdisori is larger than DISORI threshold
		struct ori anori;
		anori.bunge1 = bunge[0];
		anori.bunge2 = bunge[1];
		anori.bunge3 = bunge[2];

		anori.q0 = qbunge[0];
		anori.q1 = qbunge[1];
		anori.q2 = qbunge[2];
		anori.q3 = qbunge[3];

		anori.closestideal = get_closest_standardlage( qbunge ); //RANDOM or one of our components

		oripool.push_back( anori );

#ifdef DETAILED_PROMPTS
		cout << "ADDED to pool\t\t" << (oripool.size() - 1) << "\t\t" << anori.bunge1 << "\t\t" << anori.bunge2 << "\t\t" << anori.bunge3 << "\t\t" << anori.q0 << "\t\t" << anori.q1 << "\t\t" << anori.q2 << "\t\t" << anori.q3 << endl;
		cout << "Node;" << this->myRank << ", I have categorized close to " << anori.closestideal << endl;
#endif

		//in fact the input orientation was far enough disjoint from all previously made known orientations
		return (oripool.size() - 1);
	}
	//it was close to an already existent orientation
//cout << "GET\t\t" << closestid << "\t\t" << closestdisori << endl;
	
	return closestid;
}


void polyxx::init_adjustPRNGs ( void )
{
	//MK::ok, as tested for the parkMiller and LEcuyer generators, localized seeds for different results at the local level
	long newseed = ((long) (-1 * myRank) - 1);
	rndlocal.init ( newseed );
        
       
        //all processes generate the same large microstructure from which they sample the nuclei
        rndlargedefms.init( DEFAULT_SEED );

#ifdef DETAILED_PROMPTS
		cout << myRank << " new randomNumber seed is:" << newseed << endl;
#endif

}


void polyxx::init_MPIDatatypes( void )
{
	//MPI_IO_IdealOri_Type
	int ecntsIdeal[2] = {8, 1};
	MPI_Aint displIdeal[2] = { offsetof( MPI_IO_IdealOri, bunge1 ), offsetof( MPI_IO_IdealOri, id ) };
	MPI_Datatype oldIdeal[2] = { MPI_DOUBLE, MPI_LONG };
	MPI_Type_create_struct( 2, ecntsIdeal, displIdeal, oldIdeal, &MPI_IO_IdealOri_Type );

	//MPI_IO_Defgrain_Type
	int ecntsDG[2] = {1, 3};
	MPI_Aint displDG[2] = { offsetof( MPI_IO_Defgrain, ori ), offsetof( MPI_IO_Defgrain, rho0 ) };
	MPI_Datatype oldDG[2] = { MPI_LONG, MPI_DOUBLE };
	MPI_Type_create_struct( 2, ecntsDG, displDG, oldDG, &MPI_IO_Defgrain_Type );

	//MPI_IO_Rxgrain_Type
	int ecntsRXG[2] = {1, 1};
	MPI_Aint displRXG[2] = { offsetof( MPI_IO_Rxgrain, ori ), offsetof( MPI_IO_Rxgrain, tincub ) };
	MPI_Datatype oldRXG[2] = { MPI_LONG, MPI_DOUBLE };
	MPI_Type_create_struct( 2, ecntsRXG, displRXG, oldRXG, &MPI_IO_Rxgrain_Type );

	//MPI_IO_Ori_Type
	int ecntsOri[2] = {7, 1};
	MPI_Aint displOri[2] = { offsetof( MPI_IO_Ori, bunge1 ), offsetof( MPI_IO_Ori, closestideal ) };
	MPI_Datatype oldOri[2] = { MPI_LONG, MPI_DOUBLE };
	MPI_Type_create_struct( 2, ecntsOri, displOri, oldOri, &MPI_IO_Ori_Type );

	//MPI_IO_Temperature_Type
	int ecntsTt[1] = {2};
	MPI_Aint displTt[1] = { offsetof( MPI_IO_Temperature, t ) };
	MPI_Datatype oldTt[1] = { MPI_DOUBLE };
	MPI_Type_create_struct( 1, ecntsTt, displTt, oldTt, &MPI_IO_Temperature_Type );

	//MPI_IO_ZenerDrag_Type
	int ecntsZd[1] = {2};
	MPI_Aint displZd[1] = { offsetof( MPI_IO_ZenerDrag, t ) };
	MPI_Datatype oldZd[1] = { MPI_DOUBLE };
	MPI_Type_create_struct( 1, ecntsZd, displZd, oldZd, &MPI_IO_ZenerDrag_Type );

	int ecntsPSched[1] = {2};
	MPI_Aint displPSched[1] = { offsetof( MPI_IO_PrintSchedule, id) };
	MPI_Datatype oldPSched[1] = { MPI_LONG };
	MPI_Type_create_struct( 1, ecntsPSched, displPSched, oldPSched, &MPI_IO_PrintSchedule_Type);

	//make these usertypes known to me and the MPI world
	MPI_Type_commit(&MPI_IO_IdealOri_Type);
	MPI_Type_commit(&MPI_IO_Defgrain_Type);
	MPI_Type_commit(&MPI_IO_Rxgrain_Type);
	MPI_Type_commit(&MPI_IO_Ori_Type);
	MPI_Type_commit(&MPI_IO_Temperature_Type);
	MPI_Type_commit(&MPI_IO_ZenerDrag_Type);
	MPI_Type_commit(&MPI_IO_PrintSchedule_Type);
}


bool polyxx::init_readparameter( const char * inputfname ) 
{
	bool success = true;
	if ( myRank == MASTER ) {
		std::cout << "Master attempting to read input file" << endl;
		cout<<"inputfname="<<inputfname<<endl;
		if( !open(READ, inputfname, &(rxpathtracer_input)) ) { success = false; std::cout << "ERROR::Reading input file." << endl; } //nothing to close
	}
	MPI_Bcast( &success, 1, MPI_CHAR, MASTER, MPI_COMM_WORLD );

if ( success != true ) { cout << "ERROR::Somewhere during opening the input file." << endl; return false; }

	//MK::no return okay, master reads, other prepare data of fixed sized array in theknown of the input file format 
	success = true;
	if ( myRank == MASTER ) {
		if ( success == true )	success = read_parameter();
		if ( success == true )	success = read_processing();
		if ( success == true )	success = read_zener();
		if ( success == true )	success = read_defgpool();
		if ( success == true )	success = read_rxgpool();
		if ( success == true )  success = read_nucprintingschedule();
		fclose( rxpathtracer_input );
	}
	MPI_Bcast( &success, 1, MPI_CHAR, MASTER, MPI_COMM_WORLD );
if ( success != true ) { cout << "ERROR::Somewhere during reading the input files." << endl; return false; }

	long* pieces = new long[9];
	if ( myRank == MASTER ) {
		pieces[0] = this->nstandardlagen;
		pieces[1] = this->ndefgpool;
		pieces[2] = this->nrxgpool;
		pieces[3] = this->noripool;
		pieces[4] = this->time.size();
		pieces[5] = this->zenertime.size();
		pieces[6] = this->oripool_firstdisjunct_rxg;
		pieces[7] = this->oripool.size();
		pieces[8] = numberOfNucsToPrint;
	}
	MPI_Bcast( pieces, 9, MPI_LONG, MASTER, MPI_COMM_WORLD );
	
	//allocate memory
	bcast_ideal = new MPI_IO_IdealOri[ pieces[0] ];
	bcast_dg = new MPI_IO_Defgrain[ pieces[1] ];
	bcast_rxg = new MPI_IO_Rxgrain[ pieces[2] ];
	bcast_ori = new MPI_IO_Ori[ pieces[3] ];
	bcast_Tt = new MPI_IO_Temperature[ pieces[4] ];
	bcast_zener = new MPI_IO_ZenerDrag[ pieces[5] ];
	//6 and 7 initialized by workers themselves
	bcast_nucprintingschedule = new MPI_IO_PrintSchedule[ pieces[8] ];
	if ( MasterParameterBcast() != true )  { cout << myRank << "-t hrank throws during Bcast parameter." << endl; return false; }
	if ( AllPrepareForInputBcast() != true ) { cout << myRank << "-th rank throws during Bcast prepare." << endl; return false; }
	
	MPI_Bcast( bcast_ideal, pieces[0], MPI_IO_IdealOri_Type, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( bcast_dg, pieces[1], MPI_IO_Defgrain_Type, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( bcast_rxg, pieces[2], MPI_IO_Rxgrain_Type, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( bcast_ori, pieces[3], MPI_IO_Ori_Type, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( bcast_Tt, pieces[4], MPI_IO_Temperature_Type, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( bcast_zener, pieces[5], MPI_IO_ZenerDrag_Type, MASTER, MPI_COMM_WORLD );

	if(pieces[8] > 0){
	MPI_Bcast( bcast_nucprintingschedule, pieces[8], MPI_IO_PrintSchedule_Type , MASTER, MPI_COMM_WORLD );
	}
	if ( WorkerInterpretInputBcast( pieces ) != true ) { cout << myRank << "-th rank throws during postprocessing Bcasts." << endl; return false; }
	
	delete [] bcast_ideal;
	delete [] bcast_dg;
	delete [] bcast_rxg;
	delete [] bcast_ori;
	delete [] bcast_Tt;
	delete [] bcast_zener;
	delete [] bcast_nucprintingschedule;
	bcast_ideal = NULL;
	bcast_dg = NULL;
	bcast_rxg = NULL;
	bcast_ori = NULL;
	bcast_Tt = NULL;
	bcast_zener = NULL;
	delete [] pieces;

	//now that all data have been read let the master distribute them
	/*if ( myRank == MASTER ) {
		setup_DataArraysfor_MPI_Send();
		setup_ParameterArraysfor_MPI_Send();
	} 
	else {
		setup_ParameterArraysfor_MPI_Recv();
	}


//--MPI COMMUNICATION------------------------------------
	//Broadcast input-parameters to other Nodes
	MPI_Bcast(poly->DoubleParameter,NDOUBLEPARAMETER,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
	MPI_Bcast(poly->LongParameter,NLONGPARAMETER,MPI_LONG,MASTER,MPI_COMM_WORLD);
	//Save parameter from array to certain place and alloc for dynamic-data
	poly->init_Parameters_from_MPI_data();
	//Now broadcast all dynamic-sized Data e.g. Temperatures, oripool...
	if( poly -> myRank == MASTER){
		MPI_Bcast(&(poly->MPI_DoubleDatav[0]),poly->nDoubleDataList,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
		MPI_Bcast(&(poly->MPI_LongDatav[0]),poly->nLongDataList,MPI_LONG,MASTER,MPI_COMM_WORLD);
		poly->MPI_DoubleDatav.clear();
		poly->MPI_LongDatav.clear();
	} else {
		MPI_Bcast(poly->MPI_DoubleP,poly->nDoubleDataList,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
		MPI_Bcast(poly->MPI_LongP,poly->nLongDataList,MPI_LONG,MASTER,MPI_COMM_WORLD);
		poly->init_DynamicArrays_from_MPI_data();
	}

cout<<"RANK: "<<poly->myRank<<" RECEIVED MODEL DATA "<<endl;
cout<<"RANK: "<<poly->myRank<<" nTemp = "<<poly->nTemp<<" nZener = "<<poly->nZener<<" ndefgpool = "<<poly->ndefgpool<<" nrxgpool = "<<poly->nrxgpool<<" noripool="<<poly->oripool.size()<<endl;
*/
	
	init_someData();
	
	if (myRank == MASTER) { cout << "All parameter were sucessfully read and broadcast." << endl; }

	return true;
}


bool polyxx::MasterParameterBcast( void ) 
{
	//bcast them all ##MK can be optimized to reduce latencies ... units are correctly interpreted here
	MPI_Bcast( &meshRecursion, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &minDist, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &nCAEnsemble, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );
	//MPI_Bcast( &nNeighborsPerCA, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &grainPrintingEnabled, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );
	//MPI_Bcast( &Vref, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &ppdensity, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
        MPI_Bcast( &ppscaling, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &boxx, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &boxy, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &boxz, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &defgmedian_rd, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &defgmedian_nd, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &defgmedian_td, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );

	long nCoefs;
	if (myRank == MASTER) nCoefs = physConstants.get_NumberOfCoefficients();
	MPI_Bcast( &nCoefs, 1, MPI_LONG, MASTER, MPI_COMM_WORLD);

	double *phys_coefficients = new double [nCoefs];
	if (myRank == MASTER) physConstants.get_Coefficients(phys_coefficients);
	MPI_Bcast( phys_coefficients, nCoefs, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	if (myRank != MASTER) { physConstants.init(phys_coefficients);	}
	delete [] phys_coefficients;

	MPI_Bcast( &LAGBm0, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &LAGBHact, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &HAGBm0, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &HAGBHact, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &GSm0, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &GSHact, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );

	//rollett humphreys model
        MPI_Bcast( &RHHAGBm0, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
        MPI_Bcast( &RHHAGBHact, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
        MPI_Bcast( &RHLAGBHAGBCut, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
        MPI_Bcast( &RHLAGBHAGBTrans,  1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
        MPI_Bcast( &RHLAGBHAGBExponent, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );

	MPI_Bcast( &rvoption, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &kuhlmann_alpha, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
	MPI_Bcast( &kuhlmann_beta, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );

	MPI_Bcast( &zeneroption, 1, MPI_LONG, MASTER, MPI_COMM_WORLD ); 
	MPI_Bcast( &zenerfac, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD ); 
	
        MPI_Bcast( &mobilitymodel_option, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );
        
        MPI_Bcast( &nRediscrEnsembleRealTimeIntervals, 1, MPI_LONG, MASTER, MPI_COMM_WORLD ); 

	MPI_Bcast( &(grvolhist.VolDistBinSize), 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD ); 
	MPI_Bcast( &(grvolhist.VolDistMaxValue), 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD ); 

	return true;
}


bool polyxx::AllPrepareForInputBcast( void ) 
{
	// pack information
	if ( myRank != MASTER ) {
		//worker have nothing to do they are waiting for master, which
	}
	else {
		//packs ideals...
		for ( long i = 0; i < nstandardlagen; i++ ) {
			bcast_ideal[i].bunge1 = standardlagen[i].bunge1;
			bcast_ideal[i].bunge2 = standardlagen[i].bunge2;
			bcast_ideal[i].bunge3 = standardlagen[i].bunge3;
			bcast_ideal[i].q0 = standardlagen[i].q0;
			bcast_ideal[i].q1 = standardlagen[i].q1;
			bcast_ideal[i].q2 = standardlagen[i].q2;
			bcast_ideal[i].q3 = standardlagen[i].q3;
			bcast_ideal[i].scatter = standardlagen[i].scatter;
			bcast_ideal[i].id = standardlagen[i].id;
		}

		//pack defg...
		for ( long dg = 0; dg < ndefgpool; dg++) {
			bcast_dg[dg].ori = defgpool[dg].ori;
			bcast_dg[dg].rho0 = defgpool[dg].rho0;
			bcast_dg[dg].rho = defgpool[dg].rho;
			bcast_dg[dg].kuhlmann_tcrit = defgpool[dg].kuhlmann_tcrit;
		}

		//pack rxg
		for ( long rxg= 0; rxg < nrxgpool; rxg++ ) {
			bcast_rxg[rxg].ori = rxgpool[rxg].ori;
			bcast_rxg[rxg].tincub = rxgpool[rxg].tincub;
		}

		//pack ori
		for ( long orr = 0; orr < noripool; orr++ ) {
			bcast_ori[orr].bunge1 = oripool[orr].bunge1;
			bcast_ori[orr].bunge2 = oripool[orr].bunge2;
			bcast_ori[orr].bunge3 = oripool[orr].bunge3;
			bcast_ori[orr].q0 = oripool[orr].q0;
			bcast_ori[orr].q1 = oripool[orr].q1;
			bcast_ori[orr].q2 = oripool[orr].q2;
			bcast_ori[orr].q3 = oripool[orr].q3;
			bcast_ori[orr].closestideal = oripool[orr].closestideal;
		}

		//pack Tt
		for ( long tt = 0; tt < time.size(); tt++ ) {
			bcast_Tt[tt].t = time[tt];
			bcast_Tt[tt].T = temperature[tt];
		}

		//pack Zener
		for ( long zz = 0; zz < zenertime.size(); zz++ ) {
			bcast_zener[zz].t = zenertime[zz];
			bcast_zener[zz].pz = zenerforce[zz];
		}
		//pack NucPrinting
		for (long p = 0; p < numberOfNucsToPrint; p++ ) {
			bcast_nucprintingschedule[p].id = NucIDsforPrinting[p];
			bcast_nucprintingschedule[p].period = NucPrintingPeriod[p];
		}
	}
	
	return true;
}


bool polyxx::WorkerInterpretInputBcast( long * elements )
{
	if ( myRank != MASTER ) { //interpret the bcasts
		//ideals
		standardlagen = NULL;
		nstandardlagen = elements[0];
		standardlagen = new ideal[nstandardlagen];
		if ( standardlagen == NULL ) { cout << "ERROR::During memory allocation in workerinterpret bcast standardlagen" << endl; return false; }
		this->myMemGuard += nstandardlagen * sizeof(ideal);
		worker_interp_ideal();
	

		//defgpool
		ndefgpool = elements[1];
		defgpool = new defg[ndefgpool];
		if ( defgpool == NULL ) { cout << "ERROR::During memory allocation in workerinterpret bcast defgp" << endl; return false; }
		this->myMemGuard += ndefgpool * sizeof(defg);
		worker_interp_defgpool();
		oripool_firstdisjunct_rxg = elements[6];
		

		//rxgpool
		nrxgpool = elements[2];
		rxgpool = NULL;
		rxgpool = new rxg[nrxgpool];
		if ( rxgpool == NULL ) { cout << "ERROR::During memory allocation in workerinterpret bcast rxgp" << endl; return false; }
		this->myMemGuard += nrxgpool * sizeof(rxg);
		worker_interp_rxgpool();
		noripool = elements[7];

		
		//pack ori
		worker_interp_oripool();
		if ( noripool != elements[3] ) { cout << "ERROR::Datapool inconsistency" << endl; return false; }

		//processing
		nTemp = elements[4];
		this->worker_interp_processing();

		//Zener
		nZener = elements[5];
		worker_interp_zener();
		//Printing
		numberOfNucsToPrint = elements[8];
		worker_interp_printing();
	}
	//else {}	//nothing to do for the master

	MPI_Barrier ( MPI_COMM_WORLD );

	return true;
}
	


void polyxx::init_someData( void )
{
	kuhlmann_consider = false;
	if (rvoption == 1) kuhlmann_consider = true;

	zener_consider = false;
	if (zeneroption == 1) zener_consider = true;

	temperature_isothermalonly = false;
	if ( quickisothermal_option == 1 ) temperature_isothermalonly = true;

	if ( temperature_isothermalonly == true ) {
		if (nTemp < 1) { cout << "ERROR::Isothermal processing desired but no data provided in Processing uds file." << endl; return; }
		Tiso = temperature[0];
	}
}


bool polyxx::read_parameter ( void ) 
{
	dataBlockP parameter = readDataBlock("ParameterList", rxpathtracer_input);
	
	meshRecursion = geTInt("MeshRecursion", parameter);
	minDist = geTReal("MinimumDistanceTraveled", parameter); //MK::now is in MICROMETER, was * 1e-6;
	if ( minDist <= MINDIST_DEFAULT ) { cout << "ERROR::minDist too small!" << endl; return false; }

	nCAEnsemble = geTInt("NucleiTested", parameter);
	if ( nCAEnsemble <= 0 ) { cout << "ERROR::no nuclei devised!" << endl; return false; }

	//nNeighborsPerCA = geTInt("NumberOfNeighborsPerCA", parameter);
	//if ( nNeighborsPerCA <= MINNEIGHBORS_DEFAULT ) { cout << "ERROR::too less neighbors devised!" << endl; return false; }
	grainPrintingEnabled = geTInt("GrainPrintingEnabled", parameter);
	//Vref = geTReal("VReference", parameter) * CUBE(MICROMETER_PER_METER);
	//if ( Vref <= 0.0 ) { cout << "ERROR::Reference volume negative or zero!" << endl; return false; }

	ppdensity = geTReal("NucleationDensity",parameter) * CUBE(1e-6); //in 1/m^3 transform to 1/micron^3
        ppscaling = geTReal("PoissonCubeScaling", parameter); //in micron
        
        //old approach::##MK::double boxSize = 1E6* pow((nNeighborsPerCA + 1) / ppdensity, ONETHIRD);//because density was in 1/m^3 boxsize in microns
        //cout << "Old box size = " << (1E6* pow((125 + 1) / ppdensity, ONETHIRD)) << endl;
        //new approach
	double boxSize = geTReal("BoxSize", parameter); //in MICRON
	boxx = boxSize;
	boxy = boxSize;
	boxz = boxSize;

	//boxx = geTReal("LocalBoxDimXRD", parameter); //now is in MICROMETER was in * 1e-6;
	//boxy = geTReal("LocalBoxDimYND", parameter);
	//boxz = geTReal("LocalBoxDimZTD", parameter);
	if (boxx <= 0.0 || boxx >= MAXBOXSIZE_DEFAULT || boxy <= 0.0 || boxy >= MAXBOXSIZE_DEFAULT || boxz <= 0.0 || boxz >= MAXBOXSIZE_DEFAULT ) { cout << "ERROR::Invalid box dimensions!" << endl; return false; }
cout << "Cubic Box-Size = " << boxSize << endl;

	//##MK::20140211PARA all CA initially have the same size
	defgmedian_rd = geTReal("MedianGrainSizeinRD", parameter); //in MICRON; // RD - x
	defgmedian_nd = geTReal("MedianGrainSizeinND", parameter); // ND - y
	defgmedian_td = geTReal("MedianGrainSizeinTD", parameter); // TD - z
	if ( defgmedian_rd < MINIMUM_DEFGRAIN_MEDIAN || defgmedian_nd < MINIMUM_DEFGRAIN_MEDIAN || defgmedian_td < MINIMUM_DEFGRAIN_MEDIAN ) { cout << "ERROR::Invalid deformed grain size!" << endl; return false; }

	double G0, dGdT, Tmelt;
	double b0, dbdT_C, dbdT_a, dbdT_b;
	G0 = geTReal("ShearModulusAt0K", parameter) * 1E-12; //is in [GPa] = [1E9 Pa = 1E9 N/m^2] transform to Pa/micron^2)
	dGdT = geTReal("ShearModulusTempCoeff", parameter) * 1E-12; //is in [GPa/K] = [1E3 kg/(K s^2 um)]
	Tmelt = geTReal("MeltingTemperature", parameter) + TOFFSET; //given in Degree Celsius
	
	b0 = geTReal("BurgersVectorAt0DegCelsius", parameter) * 1E6; //given in m so multiply with 1E6 to get in micron
	dbdT_C = geTReal("LatticeExpansionC", parameter) * 1e-6 ;
	dbdT_a = geTReal("LatticeExpansiona", parameter) ; // in [1/K]
	dbdT_b = geTReal("LatticeExpansionb", parameter) ; // in [1/K^2] 

	physConstants.init(G0,dGdT,Tmelt,b0,dbdT_C,dbdT_a,dbdT_b);

	LAGBm0 = geTReal("LAGBm0", parameter) * 1E18;//given in m^4/Js = m^4/(Nms) = m^2s/kg in [um^2 s/kg], given in [m^4/Js] so multiply with 1E12
	LAGBHact = geTReal("LAGBHact", parameter) * ECHARGE;
	HAGBm0 = geTReal("HAGBm0", parameter) * 1E18;//m0 in [um^2 s/kg], given in [m^4/Js] so multiply with 1E12
	HAGBHact = geTReal("HAGBHact", parameter) * ECHARGE;
	GSm0 = geTReal("GSm0", parameter) * 1E18;//m0 in [um^2 s/kg], given in [m^4/Js] so multiply with 1E12
	GSHact = geTReal("GSHact", parameter) * ECHARGE;
	
        RHHAGBm0 = geTReal("RHModelHAGBm0", parameter) * 1E18;
        RHHAGBHact = geTReal("RHModelHAGBHact", parameter) * ECHARGE;
        RHLAGBHAGBCut = geTReal("RHModelLAGBHAGBCut", parameter);
        RHLAGBHAGBTrans = geTReal("RHModelLAGBHAGBTrans", parameter);
        RHLAGBHAGBExponent = geTReal("RHModelLAGBHAGBExponent", parameter);
        
        rvoption = geTInt("RVKuhlmannConsider", parameter);
	kuhlmann_alpha = geTReal("RVKuhlmannAlpha", parameter) * SQR(1E-6); //is in 1/m^2 transform into 1/µm^2
	kuhlmann_beta = geTReal("RVKuhlmannBeta", parameter);

	zeneroption = geTInt("ZenerConsider", parameter);
	zenerfac = geTReal("ZenerFactor", parameter);//in [J/m^2]=[kg/s^2] no size dependence
	nRediscrEnsembleRealTimeIntervals = geTInt("NumberOfTimeIntervals",parameter);
	if ( nRediscrEnsembleRealTimeIntervals < 2 || nRediscrEnsembleRealTimeIntervals < MINIMUM_REDISCR_TIMEINTERVALS ) { cout << "Invalid time rediscretization!" << endl; return false; }
	grvolhist.VolDistBinSize = geTReal("VolumeDistributionBinSize",parameter);
	grvolhist.VolDistMaxValue = geTReal("VolumeDistributionMaxValue",parameter);
	
	binarymatrix_outputoption = geTInt("BinaryMatrixInformation", parameter);
	if (binarymatrix_outputoption != NUC_VOLUME_LONG && binarymatrix_outputoption != NUC_VOLUME_DOUBLE && binarymatrix_outputoption != NUC_SURFACE_LONG) { cout << "Invalid arguments in output option, double is used instead." << endl; binarymatrix_outputoption = NUC_VOLUME_DOUBLE; }

	quickisothermal_option = geTInt("OnlyQuickIsothermal", parameter );

        
        mobilitymodel_option = SEBALD_GOTTSTEIN;
        if ( geTInt("MobilityModel", parameter) == 2 ) mobilitymodel_option = ROLLETT_HUMPHREYS;
	//##MK::solute drag still missing


	//import ideal texture components
	dataBlockP idealcomp = readDataBlock( "IdealComponents", rxpathtracer_input);
	dataLineP iline = idealcomp->first;
	nstandardlagen = idealcomp->lineCount;
	standardlagen = new ideal[nstandardlagen];

	this->myMemGuard += nstandardlagen * sizeof(ideal);

	for ( long ii = 0; ii < nstandardlagen; ii++) {
		double eul[3], qbunge[4], scatt;

		eul[0] = getReal( iline, 2 ) / 180.0 * _PI_;
		eul[1] = getReal( iline, 3 ) / 180.0 * _PI_;
		eul[2] = getReal( iline, 4 ) / 180.0 * _PI_;

		scatt = getReal( iline, 5 ) / 180.0 * _PI_;

		euler2quaternion ( eul, qbunge );

		standardlagen[ii].bunge1 = eul[0];
		standardlagen[ii].bunge2 = eul[1];
		standardlagen[ii].bunge3 = eul[2];
		standardlagen[ii].q0 = qbunge[0];
		standardlagen[ii].q1 = qbunge[1];
		standardlagen[ii].q2 = qbunge[2];
		standardlagen[ii].q3 = qbunge[3];
		standardlagen[ii].scatter = scatt;

#ifdef DETAILED_PROMPTS
cout << "Idealcomponent imported::bunge123;q0123;scatt;" << standardlagen[ii].bunge1 << ";" << standardlagen[ii].bunge2 << ";" << standardlagen[ii].bunge3 << ";" << standardlagen[ii].q0 << ";" << standardlagen[ii].q1 << ";" << standardlagen[ii].q2 << ";"<< standardlagen[ii].q3 << ";" << standardlagen[ii].scatter << endl;
#endif
		iline = iline->next;
	}

	return true;
}


bool polyxx::read_processing ( void )
{
	//WHAT IS READ HERE:
	//time/temperature vectors with SAMESIZE, contain size information (long)
	dataBlockP timetemp = readDataBlock( "ProcessingSchedule", rxpathtracer_input );
	long nts = timetemp->lineCount;
	dataLineP line = timetemp->first;

	//MK::pairs an option but scanning for current data is only in time so utilize cacheline wider
	time.reserve( nts );
	temperature.reserve ( nts );

	double ti, Ti;

	for (long ts = 0; ts < nts; ts++) {
		ti = getReal(line, 1);
		Ti = getReal(line, 2) + TOFFSET; //celsius degree into Kelvin

		//###MK::could add consistency checks

		time.push_back( ti );
		temperature.push_back( Ti );

#ifdef DETAILED_PROMPTS
cout << "Node;" << this->myRank << ";ti;Ti;" << ti << ";" << Ti << "in the vector instead time[ts];temperature[ts];" << "\ts=" << ts << "\t" << time[ts] << ";" << temperature[ts] << endl;
#endif
		line = line->next;
	}
	if ( time.size() != temperature.size() ) { cout << "ERROR::information inconstency during reading input tprocessing." << endl; return false; }

	nTemp = time.size();

	//##only isothermal processing desired?
#ifdef DETAILED_PROMPTS
	cout << "Processing profile successfully imported." << endl;
#endif

	return true;
}


bool polyxx::read_zener ( void )
{
	//WHAT IS READ HERE?
	//zenertime, zenerforce vectors (DOUBLE), SAMESIZE, size-information in vector
	//ADDITIONAL:
	//NOTHING
	dataBlockP timefr = readDataBlock( "ZenerDrag", rxpathtracer_input );
	long ntsz = timefr->lineCount;
	dataLineP line = timefr->first;

	zenertime.reserve ( ntsz );
	zenerforce.reserve ( ntsz );

	double ti, pzi;

	for (long ts = 0; ts < ntsz; ts++) { //##MK::transform in micron
		ti = getReal(line, 1);
		pzi = getReal(line, 2);
		pzi *= zenerfac;

		//###MK::could add consistency checks
		zenertime.push_back( ti );
		zenerforce.push_back( pzi );

#ifdef DETAILED_PROMPTS
cout << "Node;" << this->myRank << ";ti;pzi;" << ti << ";" << pzi << "in the vector instead zenertime[ts];zenerforce[ts];" << "\ts=" << ts << "\t" << zenertime[ts] << ";" << zenerforce[ts] << endl;
#endif
		line = line->next;
	}
	if ( zenertime.size() != zenerforce.size() ) { cout << "ERROR::information inconstency during reading input Zener forces." << endl; return false; }

	nZener=zenertime.size();

#ifdef DETAILED_PROMPTS
	cout << "Zener drag evolution successfully imported." << endl;
#endif

	return true;
}


bool polyxx::read_defgpool (void )
{
	//WHAT IS READ HERE?
	//ndefgpool (long/int)-> size of defgpool -> storage table
	//defgpool (defg-type: 3*double, long, 3*double) of size ndefgpool
	//oripool_firstdisjunct_rxg (long) -> storage table
	dataBlockP defgblock = readDataBlock("DeformedGrainsPool", rxpathtracer_input);
	ndefgpool = defgblock->lineCount;
	dataLineP line = defgblock->first;

	//MK::this pool is static and shared by all nodes
	oripool.reserve(ndefgpool);

	//denotes deformed grains that differ in their properties such as deformation degree, dislocation density and orientation
	defgpool = NULL;
	defgpool = new defg[ndefgpool];
	if ( defgpool == NULL ) { cout << "ERROR::During memory allocation in init_defgpool" << endl; return false; }
	this->myMemGuard += ndefgpool * sizeof(defg);

	long id;
	double bunge[3];
	double q[4] = {1.0, 0.0, 0.0, 0.0};

	for (long dg = 0; dg < ndefgpool; dg++) {
		bunge[0] = getReal(line, 1) / 180.0 * _PI_;
		bunge[1] = getReal(line, 2) / 180.0 * _PI_;
		bunge[2] = getReal(line, 3) / 180.0 * _PI_;

		id = check_disjunctness( bunge );

		defgpool[dg].ori = id;
		defgpool[dg].rho0 = getReal(line, 4) * 1E-12; //given in 1/m^2 = 1E-12/um^2 so multiply with 1E-12
		defgpool[dg].rho = defgpool[dg].rho0; //MK::initially rho = rho0

		if ( defgpool[dg].rho > rhoMax) {
			rhoMax = defgpool[dg].rho;
		}

		defgpool[dg].kuhlmann_tcrit = INFINITE;
		if (this->kuhlmann_consider == true) {
			defgpool[dg].kuhlmann_tcrit = ( exp( -1 * ( KUHLMANNFINALRHO - defgpool[dg].rho0 ) / kuhlmann_alpha) - 1 ) / kuhlmann_beta;
		}

		line = line->next;
	} // for all grains in the pool

	//now all disjunct deformed grain orientations are known, mark where the rx grains dereferencing starts!
	oripool_firstdisjunct_rxg = oripool.size();

	return true;
}


bool polyxx::read_rxgpool ( void )
{
	//WHAT IS READ HERE?
	//nrxgpool (long): size of rxgpool
	//rxgpool: Object of rxg-type of nrxgpool
	//rxg-type: long, double
	//+++++++++++//
	//ADDITIONAL: NOTHING
	dataBlockP rxgblock = readDataBlock("RecrystallizingGrainsPool", rxpathtracer_input);
	nrxgpool = rxgblock->lineCount;
	dataLineP line = rxgblock->first;

	//MK::there are already many orientations in the oripool most likely some are not significantly disjunct from the already existing ones so no reserve again!
	rxgpool = NULL;
	rxgpool = new rxg[nrxgpool];
	if ( rxgpool == NULL ) { cout << "ERROR::During memory allocation in init_defgpool" << endl; return false; }
	this->myMemGuard += nrxgpool * sizeof(rxg);

	long id;
	double bunge[3];
	double q[4] = {1.0, 0.0, 0.0, 0.0};

	for (long rg = 0; rg < nrxgpool; rg++) {
		bunge[0] = getReal(line, 1) / 180.0 * _PI_;
		bunge[1] = getReal(line, 2) / 180.0 * _PI_;
		bunge[2] = getReal(line, 3) / 180.0 * _PI_;

		id = check_disjunctness( bunge );

		rxgpool[rg].ori = id;
		rxgpool[rg].tincub = getReal(line, 4);

#ifdef DETAILED_PROMPTS
		cout << "id;e1;e2;e3;q0;q1;q2;q3;ori;tincub;" << rg << ";" << oripool[id].bunge1 << ";" << oripool[id].bunge2 << ";" << oripool[id].bunge3 << ";" << oripool[id].q0 << ";" << oripool[id].q1 << ";" << oripool[id].q2 << ";" << oripool[id].q3 << ";" << rxgpool[rg].ori << ";" << rxgpool[rg].tincub << endl;
#endif

		line = line->next;
	} // for all grains in the pool

	noripool = oripool.size();

	cout << nrxgpool << " recrystallizing grain orientations have been successfully imported." << endl;

	return true;
}


bool polyxx::read_nucprintingschedule ( void ) {
	if(grainPrintingEnabled == 1){
	dataBlockP nucprintblock = readDataBlock("GrainPrinting", rxpathtracer_input);
	numberOfNucsToPrint = nucprintblock->lineCount;
	dataLineP line = nucprintblock->first;
	NucIDsforPrinting.resize(numberOfNucsToPrint);
	NucPrintingPeriod.resize(numberOfNucsToPrint);
	
	for (long pr = 0; pr < numberOfNucsToPrint; pr ++){
		NucIDsforPrinting[pr] = getInt(line, 1);
		NucPrintingPeriod[pr] = getInt(line, 2);
		line = line->next;
	}
	cout << numberOfNucsToPrint << " Nuclei will be printed time resolved." << endl;
	} else {
	numberOfNucsToPrint = 0;
	}
	return true;
}


void polyxx::worker_interp_ideal( void ) 
{
	for ( long ii = 0; ii < nstandardlagen; ii++) {
		standardlagen[ii].bunge1 = bcast_ideal[ii].bunge1;
		standardlagen[ii].bunge2 = bcast_ideal[ii].bunge2;
		standardlagen[ii].bunge3 = bcast_ideal[ii].bunge3;
		standardlagen[ii].q0 = bcast_ideal[ii].q0;
		standardlagen[ii].q1 = bcast_ideal[ii].q1;
		standardlagen[ii].q2 = bcast_ideal[ii].q2;
		standardlagen[ii].q3 = bcast_ideal[ii].q3;
		standardlagen[ii].scatter = bcast_ideal[ii].scatter;
		standardlagen[ii].id = bcast_ideal[ii].id;
	}
}


void polyxx::worker_interp_defgpool( void )
{
	for (long dg = 0; dg < ndefgpool; dg++) {
		defgpool[dg].ori = bcast_dg[dg].ori;
		defgpool[dg].rho0 = bcast_dg[dg].rho0;
		defgpool[dg].rho = bcast_dg[dg].rho;

		if ( defgpool[dg].rho > rhoMax) {
			rhoMax = defgpool[dg].rho;
		}

		defgpool[dg].kuhlmann_tcrit = INFINITE;
		if ( this->kuhlmann_consider == true ) { 
			defgpool[dg].kuhlmann_tcrit = bcast_dg[dg].kuhlmann_tcrit;
		}
	}
}


void polyxx::worker_interp_rxgpool ( void )
{
	for (long rg = 0; rg < nrxgpool; rg++) {
		rxgpool[rg].ori = bcast_rxg[rg].ori;
		rxgpool[rg].tincub = bcast_rxg[rg].tincub;
	}
}


void polyxx::worker_interp_oripool ( void )
{
	oripool.reserve( noripool );

	for ( long orr = 0; orr < noripool; orr++ ) {
		struct ori anori;

		anori.bunge1 = bcast_ori[orr].bunge1;
		anori.bunge2 = bcast_ori[orr].bunge2;
		anori.bunge3 = bcast_ori[orr].bunge3;
		anori.q0 = bcast_ori[orr].q0;
		anori.q1 = bcast_ori[orr].q1;
		anori.q2 = bcast_ori[orr].q2;
		anori.q3 = bcast_ori[orr].q3;
		anori.closestideal = bcast_ori[orr].closestideal;

		oripool.push_back ( anori );
	}
}


void polyxx::worker_interp_processing ( void )
{
	time.reserve( nTemp );
	temperature.reserve( nTemp );
	
	for (long ti = 0; ti < nTemp; ti++) {
		this->time.push_back( bcast_Tt[ti].t );
		this->temperature.push_back( bcast_Tt[ti].T );
	}

	if ( time.size() != temperature.size() ) { cout << "ERROR::information inconstency during reading input tprocessing." << endl; return; }
}


void polyxx::worker_interp_zener( void )
{
	zenertime.reserve( nZener );
	zenerforce.reserve( nZener );
	
	for (long zi = 0; zi < nZener; zi++) {
		this->zenertime.push_back( bcast_zener[zi].t );
		this->zenerforce.push_back( bcast_zener[zi].pz );
	}

	if ( zenertime.size() != zenerforce.size() ) { cout << "ERROR::information inconstency during reading input ZenerDrag." << endl; return; }
}


void polyxx::worker_interp_printing( void ){
	
	if(numberOfNucsToPrint > 0){
		NucIDsforPrinting.resize( numberOfNucsToPrint );
		NucPrintingPeriod.resize( numberOfNucsToPrint );
	
		for (long p = 0; p < numberOfNucsToPrint; p++) {
			NucIDsforPrinting[p] = bcast_nucprintingschedule[p].id;
			NucPrintingPeriod[p] = bcast_nucprintingschedule[p].period;
		}
	}
}

/*

void polyxx::setup_DataArraysfor_MPI_Send( void ){
	init_rwMPIData( NO );
	nDoubleDataList = MPI_DoubleDatav.size();
	nLongDataList = MPI_LongDatav.size();
	
	cout << "MASTER: Double Data Array for transfer has size = " << nDoubleDataList << endl;
	cout << "MASTER: Long   Data Array for transfer has size = " << nLongDataList << endl;
}
void polyxx::setup_ParameterArraysfor_MPI_Send( void){
	rw_MPIDoubleParameter(NO);
	rw_MPILongParameter(NO);
}
void polyxx::setup_ParameterArraysfor_MPI_Recv( void )
{
	DoubleParameter = (double*) malloc( NDOUBLEPARAMETER * sizeof(double) );
	LongParameter = (long*) malloc( NLONGPARAMETER * sizeof(long) );
}
void polyxx::init_Parameters_from_MPI_data( void ){
	if(myRank != MASTER){
		rw_MPIDoubleParameter(YES);
		rw_MPILongParameter(YES);
		//all parameter saved at right place
		//Now Allocate for Data Block-Data
		MPI_DoubleP = (double*) malloc(nDoubleDataList * sizeof(double));
		MPI_LongP = (long*) malloc(nLongDataList * sizeof(long));
		//
		standardlagen = new ideal[nstandardlagen];
		//
		time.resize(nTemp);
		temperature.resize(nTemp);
		//
		zenertime.resize(nZener);
		zenerforce.resize(nZener);
		//
		defgpool = new defg[ndefgpool];
		rxgpool = new rxg[nrxgpool];
		oripool.resize(noripool);
	}
	//
	free(DoubleParameter);
	free(LongParameter);

};
void polyxx::init_DynamicArrays_from_MPI_data( void ){
	double *dp = MPI_DoubleP;
	long *lp = MPI_LongP;
	//Pointer values will change during execution of init_rw...,
	//so deallocation won't be possible anymore 
	init_rwMPIData(YES);
	free(dp);
	free(lp);
};
void polyxx::rw_MPIDoubleParameter(bool init){
	if(init == NO) DoubleParameter = (double*) malloc(NDOUBLEPARAMETER*sizeof(double));
	for (long i = 0; i < NDOUBLEPARAMETER; i++){
		rw_DoubleStorage(i,DoubleParameter[i],init);	
	}
}
void polyxx::rw_MPILongParameter(bool init){
	if(init == NO) LongParameter = (long*) malloc(NLONGPARAMETER*sizeof(long));
	for (long i = 0; i < NLONGPARAMETER; i++){
		rw_LongStorage(i,LongParameter[i],init);	
	}
}
void polyxx::rw_MPIDoubleData( double &value, bool init) 
{
	//##MK::is it ever executed with yes ?, here we have way too much branching?
	if ( init == YES ) {
		//cout<<"*MPI_DoubleP="<<*MPI_DoubleP<<endl;
		value = *MPI_DoubleP;
		MPI_DoubleP++;
	} 
	else {
		MPI_DoubleDatav.push_back(value);
	}
}
void polyxx::rw_MPILongData(long &value, bool init){
	if(init==YES){
		value = *MPI_LongP;
		MPI_LongP ++;
	} else {
		MPI_LongDatav.push_back(value);
	}
}
bool polyxx::rw_LongStorage ( long p_ID, long &value, bool init ){
	long *p;
	switch(p_ID){
	case 0 : p = &modelType; break;
	case 1 : p = &meshRecursion; break;
	case 2 : p = &nCAEnsemble; break;
	case 3 : p = &nNeighborsPerCA; break;
	case 4 : p = &rvoption; break;
	case 5 : p = &zeneroption; break;
	case 6 : p = &binarymatrix_outputoption; break;
	case 7 : p = &quickisothermal_option; break;
	case 8 : p = &oripool_firstdisjunct_rxg; break;
	case 9 : p = &nstandardlagen; break;
	case 10 : p = &nTemp; break;
	case 11 : p = &nZener; break;
	case 12 : p = &ndefgpool; break;
	case 13 : p = &nrxgpool; break;
	case 14 : p = &noripool; break;
	case 15 : p = &nDoubleDataList; break;
	case 16 : p = &nLongDataList; break;
	case 17 : p = &nRediscrEnsembleRealTimeIntervals; break;
	default : return NO;
	}
	if( init == YES){
	*p = value;
	} else {
	value = *p;
	}
	return YES;
}
bool polyxx::rw_DoubleStorage ( long p_ID, double &value, bool init ){
	double *p;
	switch(p_ID){
	case 0 : p = &minDist; break;
	case 1 : p = &Vref; break;
	case 2 : p = &boxx; break;
	case 3 : p = &boxy; break;
	case 4 : p = &boxz; break;
	//
	case 5 : p = &defgmedian_rd; break;
	case 6 : p = &defgmedian_nd; break;
	case 7 : p = &defgmedian_td; break;
	//
	case 8 : p = &defgscatter_rd; break;
	case 9 : p = &defgscatter_nd; break;
	case 10 : p = &defgscatter_td; break;
	//
	case 11 : p = &G; break;
	case 12 : p = &b; break;
	case 13 : p = &LAGBm0; break;
	case 14 : p = &LAGBHact; break;
	case 15 : p = &HAGBm0; break;
	case 16 : p = &HAGBHact; break;
	case 17 : p = &GSm0; break;
	case 18 : p = &GSHact; break;
	//
	case 19 : p = &kuhlmann_alpha; break;
	case 20 : p = &kuhlmann_beta; break;
	//
	case 21 : p = &zenerfac; break;
	case 22 : p = &rhoMax; break;
	//
	case 23 : p = &VolDistBinSize; break;
	case 24 : p = &VolDistMaxValue; break;
	default: return NO;
	}
	if( init == YES){
	*p = value;
	} else {
	value = *p;
	}
	return YES;
} 
void polyxx::init_rwMPIData( bool init )
{
	//MK::must not be called by myRank > MASTER!
	init_rwStandardlagenMPIData(init);
	init_rwTemperatureMPIData(init);
	init_rwZenerMPIData(init);
	init_rwDefgPoolMPIData(init);
	init_rwRxgPoolMPIData(init);
	init_rwOripoolMPIData(init);
}
void polyxx::init_rwStandardlagenMPIData(bool init)
{
	for ( long i = 0; i < nstandardlagen; i++) {
		//
		rw_MPIDoubleData(standardlagen[i].bunge1, init);
		rw_MPIDoubleData(standardlagen[i].bunge2, init);
		rw_MPIDoubleData(standardlagen[i].bunge3, init);
		//
		rw_MPIDoubleData(standardlagen[i].q0, init);
		rw_MPIDoubleData(standardlagen[i].q1, init);
		rw_MPIDoubleData(standardlagen[i].q2, init);
		rw_MPIDoubleData(standardlagen[i].q3, init);
		//
		rw_MPIDoubleData(standardlagen[i].scatter, init);
		//##MK::what about the id?
	}
}
void polyxx::init_rwTemperatureMPIData(bool init)
{
	for (long i = 0; i < nTemp; i++){
		rw_MPIDoubleData(time[i], init);
		rw_MPIDoubleData(temperature[i], init);
	} 
}
void polyxx::init_rwZenerMPIData(bool init)
{
	for (long i = 0; i < nZener; i++){
		rw_MPIDoubleData(zenertime[i], init);
		rw_MPIDoubleData(zenerforce[i], init);
	} 
}
void polyxx::init_rwDefgPoolMPIData(bool init)
{
	for (long dg = 0; dg < ndefgpool; dg++) {
		rw_MPILongData(defgpool[dg].ori, init);
		//
		rw_MPIDoubleData(defgpool[dg].rho, init);
		rw_MPIDoubleData(defgpool[dg].rho0, init);
		rw_MPIDoubleData(defgpool[dg].kuhlmann_tcrit, init);
	} 
}
void polyxx::init_rwRxgPoolMPIData(bool init)
{
	for (long rg = 0; rg < nrxgpool; rg++) {
		rw_MPILongData(rxgpool[rg].ori, init);
		//
		rw_MPIDoubleData(rxgpool[rg].tincub, init);
	}
}
void polyxx::init_rwOripoolMPIData(bool init){
	for (long i = 0; i < noripool; i++) {
		rw_MPIDoubleData(oripool[i].bunge1, init);
		rw_MPIDoubleData(oripool[i].bunge2, init);
		rw_MPIDoubleData(oripool[i].bunge3, init);
		//
		rw_MPIDoubleData(oripool[i].q0, init);
		rw_MPIDoubleData(oripool[i].q1, init);
		rw_MPIDoubleData(oripool[i].q2, init);
		rw_MPIDoubleData(oripool[i].q3, init);
		//
		rw_MPILongData(oripool[i].closestideal, init);
	}
}

*/


void polyxx::init_myvectors( long wpieces )
{
	//memory management
	myIDs.reserve( wpieces );
	myNuclei.reserve( wpieces );
	myBoxes.reserve ( wpieces );
	myDefMS.reserve ( wpieces );
	myNeighbors.reserve( wpieces );
	myNBorPaths.reserve ( wpieces );
	myNucRXFinalVolume.reserve ( wpieces );
	myLastTimeStepOripoolTexture.reserve ( wpieces );
	myCeaseTimeSteps.reserve( wpieces );
	myCoordinationNumbers.reserve ( wpieces );

	WhoHasWhichNucleus.reserve( nCAEnsemble );

	WhoHasHowManyNuclei = (int*) calloc( nRanks, sizeof(int) );
	WhoHasHowManyNucleiCumulated = (int*) calloc( nRanks, sizeof(int) );
}


void polyxx::init_workpartitioning_myNuclei( void )
{
	//##MK::currently all automata are of equal size so simple spread work on the nodes as equally as possible
	//MK::as a modulo operation is utilized ALL NODES come to the same partitioning scheme, therefore we do not need global synchronization up to this point
	long luckyNode;

	//smartly prepare for many things to come
	if ( myRank == MASTER ) {
		if ( (nCAEnsemble % nRanks) != 0 ) { cout << "WARNING::Choose number of nuclei as a multiple number of nRanks to prevent excessive microstopping before barriers." << endl; }
	}


	//how much work per node, rounding off via long
	long workPiece = (nCAEnsemble / nRanks);

	init_myvectors( workPiece );

	//MK::now distribute the work classical MPI modulo way
	for (long nuc = 0; nuc < nCAEnsemble; nuc++) {

		luckyNode = nuc % this->nRanks; //number % x with x = 0 is undefined, but nRanks >= 1
		WhoHasHowManyNuclei[luckyNode] = WhoHasHowManyNuclei[luckyNode] + 1;

		if ( myRank == luckyNode ) {
			myIDs.push_back ( nuc );
			//myNuclei, myBoxes, myNeighbors, myNBorPaths handled locally later
			long* defms = NULL;
			myDefMS.push_back ( defms );

#ifdef DETAILED_PROMPTS
cout << "Node;" << myRank << " I have been assigned to work on CA " << nuc << endl;
#endif
		}

		//bookkeeping
		WhoHasWhichNucleus.push_back( luckyNode ); //all know this info as modulo is univocal

	} //partition all nuclei
	int pos = 0;
	for (int Rank = 0; Rank < nRanks; Rank++){
		WhoHasHowManyNucleiCumulated[Rank] = pos;
		pos += WhoHasHowManyNuclei[Rank];
	}

	cout << "RANK: " << myRank <<" has " << WhoHasHowManyNuclei[myRank]<<" Nuclei. workPiece = " << workPiece <<" myPos = "<< WhoHasHowManyNucleiCumulated[myRank] << endl;
#ifdef DETAILED_PROMPTS
	cout << "Node;" << myRank << " has been assigned " << myIDs.size() << " work pieces." << endl;
#endif
	if (myRank == MASTER) { cout << "All nuclei desired were partitioned." << endl; }
}


void polyxx::init_poissonpointprocess( void )
{
    //on a large window, ppdensity is in nuclei per m^3
    double PoissonCubeEdgeLength = ppscaling; //micron^3
    double PoissonCubeRVE = CUBE( PoissonCubeEdgeLength ); //micron^3
    long NPointsSampled = PoissonCubeRVE * ppdensity;
    
    //generate that process
    for ( long p = 0; p < NPointsSampled; p++ ) {
        struct point ap;
        
        ap.x = rndlocal.leEcuyer() * PoissonCubeEdgeLength;
        ap.y = rndlocal.leEcuyer() * PoissonCubeEdgeLength;
        ap.z = rndlocal.leEcuyer() * PoissonCubeEdgeLength;
        
        pointprocess.push_back( ap );        
    }
    
    cout << "Point process sampled with edge/vol/npoints/pp.size() = " << PoissonCubeEdgeLength << ";" << PoissonCubeRVE << ";" << NPointsSampled << ";" << pointprocess.size() << endl;
}


double polyxx::DisoriTwoGrains ( long o1, long o2 )
{
	double q01 = oripool[o1].q0;
	double q11 = oripool[o1].q1;
	double q21 = oripool[o1].q2;
	double q31 = oripool[o1].q3;

	double q02 = oripool[o2].q0;
	double q12 = oripool[o2].q1;
	double q22 = oripool[o2].q2;
	double q32 = oripool[o2].q3;

	return misorientationCubicQxQ( q01, q11, q21, q31,  q02, q12, q22, q32 );
}


long polyxx::gbfacenucleationmodel( bool onlyat_hagb, long dg1, long dg2, double scale_dens2lamb, double scale_nuc2numb, double bndarea_micron )
{
	//implements the physics of the nucleation at grain boundary faces
	long gbnuc = NO_NUCLEI;

	//##MK::DEBUG can later be deleted
	if (dg1 < 0 || dg1 >= this->largedefms.size() ) { cout << "ACCESSVIOLATION AGAINST DEFGPOOL!" << endl; return 0; }
	if (dg2 < 0 || dg2 >= this->largedefms.size() ) { cout << "ACCESSVIOLATION AGAINST DEFGPOOL!" << endl; return 0; }

	//rhoDifference?
	double rho1 = largedefms[dg1].rho0 * 1.0e12;
	long ori1 = largedefms[dg1].ori;
	double rho2 = largedefms[dg2].rho0 * 1.0e12;
	long ori2 = largedefms[dg2].ori;

	double disori = DisoriTwoGrains( ori1, ori2 );

	//no nucleation at low-angle grain boundaries desired? then no nuclei!
	if ( onlyat_hagb == true && disori < LAGB_TO_HAGB_TRANSITION )
		return gbnuc;

	//number density is scaled according to Poisson-distribution the probability of seeding more nuclei the higher the stronger and the large drho and the large the boundary area
	//larger and higher the dislocation density difference
	double drho =  fabs(rho2 - rho1);

	//too low a difference? no nuclei here!
	if ( drho <= MINIMUM_DRHO ) 
		return gbnuc;

	//at the moment it is not important which grain the nucleus seeds as we are only interested in the position in space
	/*
	if ( rho2 > rho1 ) drho = rho2 - rho1;
	if ( rho1 < rho2 ) drho = rho1 - rho2;
	*/

	double lamb = drho * scale_dens2lamb; //MK:: was until 20150720 ( drho / drhomax) * SCALING_LAMBDA;

//cout << "gbfacenuc;dg1;dg2;disori;rho1;rho2;scale;lamb =" << dg1 << ";" << dg2 << ";" << disori << ";" << rho1 << ";" << rho2 << ";" << scale_dens2lamb << "\t\t" << lamb << endl;

	double pcum[POISSON_CUMSUM_TABLE];
	double pkcumsum = 0.0;
	for ( long k = 0; k <= POISSON_CUMSUM_CUTOFF; k++ ) {
		pcum[k] = this->poissondistribution( lamb, k );

//cout << "\t\t---" << k << "\t\t" << pcum[k] << "\t\t" << pow(lamb, k) << "\t\t" << fac(k) << "\t\t" << pkcumsum << endl;

		pkcumsum = pkcumsum + pcum[k];
	}

	if ( pkcumsum < DOUBLE_ACCURACY ) { cout << "ERROR::pkcumsum is unexpectedly close to zero " << pkcumsum << endl; return NO_NUCLEI; }

	//form cumulated distribution
	double pk = 0.0;
	for ( long k = 0; k <= POISSON_CUMSUM_CUTOFF; k++) {
		pk = pk + pcum[k];
		pcum[k] = pk / pkcumsum;

//cout << "\t\t" << k << "\t\t" << pcum[k] << endl;

	}

	//pick random number
	double luckynumber = rndlargedefms.leEcuyer();

	long j = 0; //which number to take?
	while ( (luckynumber > pcum[j]) && (j < POISSON_CUMSUM_CUTOFF) ) {
		j++;
	}

	double nnuc = j;
        
//cout << "nnuc;scale_nuc2numb;bndarea = " << nnuc << ";" << scale_nuc2numb << ";" << bndarea_micron << endl;
	nnuc = nnuc * scale_nuc2numb * bndarea_micron; //number does not translate in so many nuclei

        
        
	//cast down
	gbnuc = (long) nnuc; //truncates trailing decimals

	//therefore, are nuclei or lost or not?
	nnuc = nnuc - (double) gbnuc;

	//add another one
	if ( nnuc >= rndlargedefms.leEcuyer() ) 
		gbnuc++;

//cout << luckynumber << "\t\t" << j << "\t\t" << pcum[j] << "\t\t" << bndarea_micron << " gbnuc to place = " << gbnuc << " disori = " << setprecision(6) << disori << endl;

	return gbnuc;
}


void polyxx::init_gbface_nucleation_yz( bool onlyhagb ) 
{
        long nnx = largedefmsinfo.nx;
        long nny = largedefmsinfo.ny;
        long nnz = largedefmsinfo.nz;
        long nnxy = largedefmsinfo.nxy;
        double ddx = largedefmsinfo.dimx;
        double ddy = largedefmsinfo.dimy;
        double ddz = largedefmsinfo.dimz;
    
    //YZ planes along X
	double faceareayz_micron = ddy * ddz;
//##cout << "faceareayz_micron=" << ((dgrsizey/dgrsizex) * dgrsizex2micron) << "\t\t" << ((dgrsizez / dgrsizex) * dgrsizex2micron) << endl;
	double xpos, ypos, zpos;
	long N;
	long dgr1, dgr2, dgr1id, dgr2id, dgr_minrho_oid;
	double ix, iy, iz;
	double scale_rhodens2lamb = largedefmsinfo.gbnucdrho2dens; //dens2num
	double scale_nuc2dens = largedefmsinfo.gbnucdens2num;
        double seedmin_dgav0 = largedefmsinfo.gbnucmaxscatter;
     
	for ( long z = 0; z < nnz; z++ ) {
		zpos = z;
		zpos = zpos * ddz;

		for ( long y = 0; y < nny; y++) {
			ypos = y;
			ypos = ypos * ddy;

			for ( long x = 0; x <= nnx; x++) { //MK::ndgrx inclusive!, position of the boundary, assuming inside unit cube, left and right part are becoming cut..., xpos left edge, ypos bottom edge, zpos front edge
				xpos = x;
				xpos = xpos * ddx;

				dgr1 = x - 1;
				dgr2 = x;

				//assume deformed grain aggregate periodic
				if ( dgr1 < 0 )			dgr1 += nnx;
				if ( dgr1 >= nnx )	dgr1 -= nnx;
				if ( dgr2 < 0 )			dgr2 += nnx;
				if ( dgr2 >= nnx )	dgr2 -= nnx;

//cout << z << ";" << y << ";" << x << "\t\t" << dgr1 << "\t\t" << dgr2 << "--" << ndgrx << ";" << ndgry << ";" << ndgrz << endl;

                                
                                dgr1id = dgr1+(y*nnx)+(z*nnxy);
                                dgr2id = dgr2+(y*nnx)+(z*nnxy);
                                //minimum rho
                                dgr_minrho_oid = dgr1id;
                                if ( largedefms[dgr2id].rho0 > largedefms[dgr1id].rho0 )
                                    dgr_minrho_oid = dgr2id;
                                
                                
				N = this->gbfacenucleationmodel( onlyhagb, dgr1id, dgr2id , scale_rhodens2lamb, scale_nuc2dens, faceareayz_micron );

                                double seedminOri[3];
                                seedminOri[0] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge1;
                                seedminOri[1] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge2;
                                seedminOri[2] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge3;
                                
				//place the nuclei critical size apart
				for ( unsigned long n = 0; n < N; n++ ) {
					ix = xpos;
					iy = ypos + rndlargedefms.leEcuyer() * ddy;
					iz = zpos + rndlargedefms.leEcuyer() * ddz;

                                        struct rxgo anuc;
                                        anuc.cx = ix;
                                        anuc.cy = iy;
                                        anuc.cz = iz;
                                        
                                        double obunge[3];
                                        //##MK::specific orientation from reference
//cout << "...placing nuclei n/N/ix/iy/iz = " <<  n << ";" << N << ";" << ix << ";" << iy << ";" << iz << " SO3 rnd now ... " << seedminOri[0] << "--" << seedminOri[1] << "--" << seedminOri[2] << endl;
                                        specificallyDisorientednewOriFromReference( seedminOri, seedmin_dgav0, obunge );
//cout << "...passed SO3 rnd obunge123 = " << obunge[0] << ";" << obunge[1] << ";" << obunge[2] << endl;
                                    
                                        anuc.bunge1 = obunge[0];
                                        anuc.bunge2 = obunge[1];
                                        anuc.bunge3 = obunge[2];
                                        
                                        anuc.oid = UNKNOWN; //will be assigned later in the code
                                        anuc.rxgpoolid = UNKNOWN;
                                        
                                        largedefmsnuc.push_back( anuc );
				} //for all gbnuclei

			} //the next boundaries
		}
	}
cout << this->myRank << "-th populated the YZ boundaries, largedefmsnuc.size() = " << largedefmsnuc.size() << endl;
}


void polyxx::init_gbface_nucleation_xz( bool onlyhagb ) 
{
        long nnx = largedefmsinfo.nx;
        long nny = largedefmsinfo.ny;
        long nnz = largedefmsinfo.nz;
        long nnxy = largedefmsinfo.nxy;
        double ddx = largedefmsinfo.dimx;
        double ddy = largedefmsinfo.dimy;
        double ddz = largedefmsinfo.dimz;
        
	//XZ planes along Y
	double faceareaxz_micron = ddx * ddz;
//##cout << "faceareaxz_micron=" << ((dgrsizex/dgrsizex) * dgrsizex2micron) << "\t\t" << ((dgrsizez / dgrsizex) * dgrsizex2micron) << endl;
	double xpos, ypos, zpos;
	long N;
	long dgr1, dgr2, dgr1id, dgr2id, dgr_minrho_oid;
	double ix, iy, iz;
	double scale_rhodens2lamb = largedefmsinfo.gbnucdrho2dens; //dens2num
	double scale_nuc2dens = largedefmsinfo.gbnucdens2num;
        double seedmin_dgav0 = largedefmsinfo.gbnucmaxscatter;
        
	for ( long z = 0; z < nnz; z++ ) {
		zpos = z;
		zpos = zpos * ddz;

		for ( long x = 0; x < nnx; x++) {
			xpos = x;
			xpos = xpos * ddx;

			for ( long y = 0; y <= nny; y++) {
				ypos = y;
				ypos = ypos * ddy;

				dgr1 = y - 1;
				dgr2 = y;

				if ( dgr1 < 0 )			dgr1 += nny;
				if ( dgr1 >= nny )	dgr1 -= nny;
				if ( dgr2 < 0 )			dgr2 += nny;
				if ( dgr2 >= nny )	dgr2 -= nny;

//cout << z << ";" << y << ";" << x << "\t\t" << dgr1 << "\t\t" << dgr2 << "--" << ndgrx << ";" << ndgry << ";" << ndgrz << endl;

                                dgr1id = x+(dgr1*nnx)+(z*nnxy);
                                dgr2id = x+(dgr2*nnx)+(z*nnxy);
                                //minimum rho
                                dgr_minrho_oid = dgr1id;
                                if ( largedefms[dgr2id].rho0 > largedefms[dgr1id].rho0 )
                                    dgr_minrho_oid = dgr2id;
                                
                                
				N = this->gbfacenucleationmodel( onlyhagb, dgr1id, dgr2id , scale_rhodens2lamb, scale_nuc2dens, faceareaxz_micron );

                                double seedminOri[3];
                                seedminOri[0] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge1;
                                seedminOri[1] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge2;
                                seedminOri[2] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge3;

				//place the nuclei critical size apart
				for ( unsigned long n = 0; n < N; n++ ) {
					ix = xpos + rndlargedefms.leEcuyer() * ddx;
					iy = ypos;
					iz = zpos + rndlargedefms.leEcuyer() * ddz;

					struct rxgo anuc;
                                        anuc.cx = ix;
                                        anuc.cy = iy;
                                        anuc.cz = iz;
                                        
                                        double obunge[3];
                                        //##MK::specific orientation from reference
                                        specificallyDisorientednewOriFromReference( seedminOri, seedmin_dgav0, obunge );
                                       
                                        anuc.bunge1 = obunge[0];
                                        anuc.bunge2 = obunge[1];
                                        anuc.bunge3 = obunge[2];
                                        
                                        anuc.oid = UNKNOWN; //will be assigned later in the code
                                        anuc.rxgpoolid = UNKNOWN;
                                        
                                        largedefmsnuc.push_back( anuc );
				} //for all gbnuclei

			} //the next boundaries
		}
	}
cout << this->myRank << "-th populated the XZ boundaries, largedefmsnuc.size() = " << largedefmsnuc.size() << endl;
}


void polyxx::init_gbface_nucleation_xy( bool onlyhagb ) 
{
        long nnx = largedefmsinfo.nx;
        long nny = largedefmsinfo.ny;
        long nnz = largedefmsinfo.nz;
        long nnxy = largedefmsinfo.nxy;
        double ddx = largedefmsinfo.dimx;
        double ddy = largedefmsinfo.dimy;
        double ddz = largedefmsinfo.dimz;	
        //XY planes along Z
	double faceareaxy_micron = ddx * ddy;
//##cout << "faceareaxy_micron=" << ((dgrsizex/dgrsizex) * dgrsizex2micron) << "\t\t" << ((dgrsizey / dgrsizex) * dgrsizex2micron) << endl;
	double xpos, ypos, zpos;
	long N;
	long dgr1, dgr2, dgr1id, dgr2id, dgr_minrho_oid;
	double ix, iy, iz;
       double scale_rhodens2lamb = largedefmsinfo.gbnucdrho2dens; //dens2num
	double scale_nuc2dens = largedefmsinfo.gbnucdens2num;
         double seedmin_dgav0 = largedefmsinfo.gbnucmaxscatter;
         
	for ( long y = 0; y < nny; y++ ) {
		ypos = y;
		ypos = ypos * ddy;

		for ( long x = 0; x < nnx; x++) {
			xpos = x;
			xpos = xpos * ddx;

			for ( long z = 0; z <= nnz; z++) {
				zpos = z;
				zpos = zpos * ddz;

				dgr1 = z - 1;
				dgr2 = z;

				if ( dgr1 < 0 )			dgr1 += nnz;
				if ( dgr1 >= nnz )	dgr1 -= nnz;
				if ( dgr2 < 0 )			dgr2 += nnz;
				if ( dgr2 >= nnz )	dgr2 -= nnz;

//cout << z << ";" << y << ";" << x << "\t\t" << dgr1 << "\t\t" << dgr2 << "--" << ndgrx << ";" << ndgry << ";" << ndgrz << endl;

                                
                                dgr1id = x+(y*nnx)+(dgr1*nnxy);
                                dgr2id = x+(y*nnx)+(dgr2*nnxy);
                                //minimum rho
                                dgr_minrho_oid = dgr1id;
                                if ( largedefms[dgr2id].rho0 > largedefms[dgr1id].rho0 )
                                    dgr_minrho_oid = dgr2id;
                                
				N = this->gbfacenucleationmodel( onlyhagb, dgr1id, dgr2id , scale_rhodens2lamb, scale_nuc2dens, faceareaxy_micron );

                                double seedminOri[3];
                                seedminOri[0] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge1;
                                seedminOri[1] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge2;
                                seedminOri[2] = oripool[defgpool[largedefms[dgr_minrho_oid].defgpoolid].ori].bunge3;

				//place the nuclei critical size apart
				for ( unsigned n = 0; n < N; n++ ) {
					ix = xpos + rndlargedefms.leEcuyer() * ddx; //sx, sy, and sz are already truncated for boundary contact
					iy = ypos + rndlargedefms.leEcuyer() * ddy;
					iz = zpos;

					struct rxgo anuc;
                                        anuc.cx = ix;
                                        anuc.cy = iy;
                                        anuc.cz = iz;
                                        
                                        double obunge[3];
                                        //##MK::specific orientation from reference
                                        specificallyDisorientednewOriFromReference( seedminOri, seedmin_dgav0, obunge );
                                        
                                        anuc.bunge1 = obunge[0];
                                        anuc.bunge2 = obunge[1];
                                        anuc.bunge3 = obunge[2];
                                        
                                        anuc.oid = UNKNOWN; //will be assigned later in the code
                                        anuc.rxgpoolid = UNKNOWN;
                                        
                                        largedefmsnuc.push_back( anuc );
				} //for all gbnuclei

			} //the next boundaries
		}
	}

cout << this->myRank << "-th populated the XY boundaries, largedefmsnuc.size() = " << largedefmsnuc.size() << endl;
}


void polyxx::init_large_defms_nucleation( void )
{
        init_gbface_nucleation_yz( true ); //true - nucleation only at HAGB (dg>15°) otherwise on each boundary
        init_gbface_nucleation_xz( true );
	init_gbface_nucleation_xy( true );
        
        //register all nuclei that where added and categorize in sufficiently disjoint orientations
        //this limits the total amount of disjoint orientations that are characterized in the simulation
        
        //##MK::register orientation in oripool, ALL PROCESSES HAVE THE SAME DEFORMATION STRUCTURE AND HENCE THE SAME ORIPOOL IMAGE!
        double obunge[3];
        long oi;
        
        for (long nuc = 0; nuc < largedefmsnuc.size(); nuc++ ) {
            obunge[0] = largedefmsnuc[nuc].bunge1;
            obunge[1] = largedefmsnuc[nuc].bunge2;
            obunge[2] = largedefmsnuc[nuc].bunge3;
            
            oi = check_disjunctness ( obunge );
            largedefmsnuc[nuc].oid = oi;
        }
        //registrate all nuclei in the myrxgpool
        long oldSize = nrxgpool;
        long newSize = oldSize + largedefmsnuc.size();

		for (long nn = 0; nn < largedefmsnuc.size(); nn++ ) {
			largedefmsnuc[nn].rxgpoolid = oldSize + nn;
		}

        rxgP newrxgpool = NULL;
        newrxgpool = new rxg[newSize];
        if ( newrxgpool == NULL ) { cout << "Unsuccessful assignment in init_large_defms_nucleation!" << endl; }
        //copy old entries
        for ( long rr = 0; rr < oldSize; rr++ ) {
            newrxgpool[rr].ori = rxgpool[rr].ori;
            newrxgpool[rr].tincub = rxgpool[rr].tincub;
        }
        //add new entries
        for ( long nw = oldSize; nw < newSize; nw++ ) {
            newrxgpool[nw].ori = largedefmsnuc[nw-oldSize].oid;
            newrxgpool[nw].tincub = SITE_SATURATED;
        }
        
        //delete old pool and link in new pool
        delete [] rxgpool;
        rxgpool = newrxgpool;
        
        //update size of the pool
        nrxgpool = newSize;
        
cout << this->myRank << " init_large_defms_nucleation, my nrxgpool = " << nrxgpool << endl;
}


void polyxx::init_large_defms( void )
{
    //##MK::initializes a large cuboid 3D aggregate of grains picked randomly from the deformed structure
    
    //##DEBUG
    //set properties of this deformed structure, see also polyxx constructor
    largedefmsinfo.nx = 1; //5 grains a 1000micron in RD, large system (5000micron)^3
    largedefmsinfo.ny = 100; //500 grains a 10micron in ND
    largedefmsinfo.nz = 10; //50 grains a 100micron in TD
    largedefmsinfo.nxy = largedefmsinfo.nx * largedefmsinfo.ny;
    largedefmsinfo.nxyz = largedefmsinfo.nx * largedefmsinfo.ny * largedefmsinfo.nz;
    largedefmsinfo.dimx = defgmedian_rd;
    largedefmsinfo.dimy = defgmedian_nd;
    largedefmsinfo.dimz = defgmedian_td;
    
    //##MK::set as default here which was utilized in the PhDMK study
	largedefmsinfo.gbnucdens2num = 3.9E-4; //3.9E-4; //2.0E-4;
    largedefmsinfo.gbnucdrho2dens = 3.5436e-15;
    largedefmsinfo.gbnucmaxscatter = 8.0 / 180.0 * _PI_;
    
    //sample this structure from the available list of deformed grains, orientations of which where already introduced
    long ndg = largedefmsinfo.nxyz;
    long cand = 0;
    for ( long dg = 0; dg < ndg; dg++ ) {
        struct defglean adg;
        
        cand = rndlargedefms.leEcuyer() * ndefgpool; //[0,1) assures g to be smaller than ndefgpool
        
        adg.defgpoolid = cand;
        adg.ori = defgpool[cand].ori;
        adg.rho0 = defgpool[cand].rho0;     
        
        largedefms.push_back ( adg );
    }
    //order 3D implicit first positive x stack in y stack in xy slices positively in z

    //identify all nuclei in this structure, ##MK::NECESSARY BECAUSE IF THE NUCLEATION SITE IS NOT KNOWN WE CANNOT CENTER THE BOX...
    init_large_defms_nucleation();
    
	cout << "Init_large_defms after nucleation = " << ( MPI_Wtime() - this->realStartTimeOfComputing ) << " seconds." << endl;
  
     //##MK::oripool is fixed and known from now on!
     noripool = oripool.size();
     
cout << this->myRank << " success to generate largedefms with = " << largedefms.size() << " defgs, " << largedefmsnuc.size() << " nuclei, and an oripool of " << oripool.size() << " elements." << endl;
}


void polyxx::init_myBoxes( void )
{
	//init myNuclei handled later in the nucleation model init_myNuclei_nucleation
	if ( (boxx * boxy * boxz) <= DOUBLE_ACCURACY ) { cout << "WARNING::Potential dbl precision inaccuracy!" << endl; }

	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
		struct box bo;

		bo.xrd = boxx; //##MK::micron
		bo.ynd = boxy;
		bo.ztd = boxz;
		bo.xyz = boxx * boxy * boxz;

		//initialize defgr translation vector with NO_RELATIVE_TRANSLATION
		bo.tx = 0.0;
		bo.ty = 0.0;
		bo.tz = 0.0;

		//initially assume monocrystal filling the box
		bo.ngrx = 1;
		bo.ngry = 1;
		bo.ngrz = 1;
		bo.ngrxy = 1;
		bo.ngrxyz = 1;

		myBoxes.push_back( bo ); //copy constructor

		//enable printing for particular nuclei
		nuc_printing.push_back(0);
		for (long i = 0; i < numberOfNucsToPrint; i++){
			if( WhoHasHowManyNucleiCumulated[myRank] + mynuc != NucIDsforPrinting[i] ) continue;
			nuc_printing[mynuc] = NucPrintingPeriod[i];
		}
	}

	//paths are constructed later...

	if (myRank == MASTER) { std::cout << "All internal memory ready." << endl; }
}


void polyxx::init_myNuclei_BasedOnLargeDefMS( void )
{
	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
                //now deformed grain structure and nuclei are fixed
                //so identify in the largedefms one nucleus that is far enough away from the boundary, here no periodic boundary conditions
                //myBoxes[*].xrd is in micron
            
            bool found = false;
            //coordinates in largedefms run from [0,maxsize]^3 in micron, the nuclei are inside this all inside this region
            double box_x = myBoxes[mynuc].xrd * 0.5;
            double box_y = myBoxes[mynuc].ynd * 0.5;
            double box_z = myBoxes[mynuc].ztd * 0.5;
            double owin_extx = (double) largedefmsinfo.nx * largedefmsinfo.dimx;
            double owin_exty = (double) largedefmsinfo.ny * largedefmsinfo.dimy;
            double owin_extz = (double) largedefmsinfo.nz * largedefmsinfo.dimz;
            
            long refnucleus = UNKNOWN_CANDIDATE;
            
            while ( found == false ) {
                refnucleus = rndlocal.leEcuyer() * largedefmsnuc.size();
                
                if ( (largedefmsnuc[refnucleus].cx - box_x) < 0.0 )             continue;
                if ( (largedefmsnuc[refnucleus].cx + box_x) > owin_extx )      continue;
                if ( (largedefmsnuc[refnucleus].cy - box_y) < 0.0 )             continue;
                if ( (largedefmsnuc[refnucleus].cy + box_y) > owin_exty )       continue;
                if ( (largedefmsnuc[refnucleus].cz - box_z) < 0.0 )             continue;
                if ( (largedefmsnuc[refnucleus].cz + box_z) > owin_extz )      continue;
                
                //not continued? inside!
                found = true;
                //refnucleus can be utilized
            }
            
            if ( refnucleus < 0) cout << "ERROR::Unable to identify a reference nucleus in myRank/myNuc = " << this->myRank << "--" << mynuc << endl;
            
            //box center coordinates
            double bcx = largedefmsnuc[refnucleus].cx;
            double bcy = largedefmsnuc[refnucleus].cy;
            double bcz = largedefmsnuc[refnucleus].cz;
            //box is inside the largedefms but where is the next grain boundary?
            
            //how many (in units of entire grains) is the box located in the largedefms?
            long igx = (bcx - box_x) / defgmedian_rd;
            long igy = (bcy - box_y) / defgmedian_nd;
            long igz = (bcz - box_z) / defgmedian_td;
			//for instance left box boundary inside igx-th deformed grain in x direction...
            
               //determine first the relative displacement of the defgr aggregate to get the necessary number of defgrains which fully enclose the simulation volume

		myBoxes[mynuc].tx = (igx * defgmedian_rd) - (bcx - box_x);
                myBoxes[mynuc].ty = (igy * defgmedian_nd) - (bcy - box_y);
		myBoxes[mynuc].tz = (igz * defgmedian_td) - (bcz - box_z);
                
                //##DEBUG, have to be negative or zero in which case the box is at the boundary!
                if ( myBoxes[mynuc].tx > 0.0) cout << "myBoxes[" <<  mynuc << "] invalid tx = " << myBoxes[mynuc].tx << endl;
                if ( myBoxes[mynuc].ty > 0.0) cout << "myBoxes[" <<  mynuc << "] invalid ty = " << myBoxes[mynuc].ty << endl;
                if ( myBoxes[mynuc].tz > 0.0) cout << "myBoxes[" <<  mynuc << "] invalid tz = " << myBoxes[mynuc].tz << endl;
                
		//then how many grains in either direction are necessary?, tx, ty, tz are negative!
		myBoxes[mynuc].ngrx = (( myBoxes[mynuc].xrd - myBoxes[mynuc].tx ) / defgmedian_rd); //int rounds off so...
		myBoxes[mynuc].ngrx++;
		if ( myBoxes[mynuc].ngrx == 0 ) myBoxes[mynuc].ngrx++;

		myBoxes[mynuc].ngry = (( myBoxes[mynuc].ynd - myBoxes[mynuc].ty ) / defgmedian_nd);
		myBoxes[mynuc].ngry++;
		if ( myBoxes[mynuc].ngry == 0 ) myBoxes[mynuc].ngry++;

		myBoxes[mynuc].ngrz = (( myBoxes[mynuc].ztd - myBoxes[mynuc].tz ) / defgmedian_td);
		myBoxes[mynuc].ngrz++;
		if ( myBoxes[mynuc].ngrz == 0 ) myBoxes[mynuc].ngrz++;

		myBoxes[mynuc].ngrxy = myBoxes[mynuc].ngrx * myBoxes[mynuc].ngry;
		myBoxes[mynuc].ngrxyz = myBoxes[mynuc].ngrxy * myBoxes[mynuc].ngrz;

		//allocate memory
		long n = myBoxes[mynuc].ngrxyz;
		long* defms = NULL;
		defms = new long[n];
		if (defms == NULL) { cout << myRank << "ERROR::Allocating memory during init_myNuclei!" << endl; }
		this->myMemGuard += n * sizeof(long);

//cout << myRank << ";" << mynuc << ";txyz;" << myBoxes[mynuc].tx << ";" << myBoxes[mynuc].ty << ";" << myBoxes[mynuc].tz << ";ngrxyz;" << myBoxes[mynuc].ngrx << ";" << myBoxes[mynuc].ngry << ";" << myBoxes[mynuc].ngrz << "\tn=" << n << endl;

		//alignment of the deformed grains is as in +x aligned, stacked in +y into slabs stacked in +z, so implicit 3D linear array

		//###MK::copy the deformed grains that are located inside the box from the largedefms
                //the offsets are igx, igy, igz,
                long locdefms = 0;
                long gg = 0;
                long gz_pbc, gy_pbc, gx_pbc; //periodic boundary conditions may apply
                
                for ( long gz = 0; gz < myBoxes[mynuc].ngrz; gz++ ) { //order is the same implicit one
                    gz_pbc = gz + igz;
                    if ( gz_pbc > largedefmsinfo.nz ) gz_pbc = gz_pbc - largedefmsinfo.nz; //% largedefmsinfo.nz;
                    
                    for ( long gy = 0; gy < myBoxes[mynuc].ngry; gy++ ) {
                        gy_pbc = gy + igy;
                        if ( gy_pbc > largedefmsinfo.ny ) gy_pbc = gy_pbc - largedefmsinfo.ny; //% largedefmsinfo.ny;
                        
                        for (long gx = 0; gx < myBoxes[mynuc].ngrx; gx++) {
                            gx_pbc = gx + igx;
                            if ( gx_pbc > largedefmsinfo.nx ) gx_pbc = gx_pbc - largedefmsinfo.nz; //% largedefmsinfo.nx;
                            
                            locdefms = gx_pbc + (gy_pbc * largedefmsinfo.nx) + (gz_pbc * largedefmsinfo.nxy);
                            
                            //##DEBUG
                            if ( locdefms >= largedefms.size() ) cout << "Locdefms invalid!" << endl;
                            if ( largedefms[locdefms].defgpoolid >= this->ndefgpool ) cout << "Largedefms id invalid!" << endl;
                            //##DEBUG
                            
                            defms[gg] = largedefms[locdefms].defgpoolid; //references which of the template grains is located here!
                            gg++;
			  
//cout << "mynuc|bcxyz|igxyz|txyz|nxyz--|gxyz_pbc|locdefms||ngrxyz//"<<mynuc<<";"<<bcx<<";"<<bcy<<";"<<bcz<<"|"<<igx<<";"<<igy<<";"<<igz<<"|"<<myBoxes[mynuc].tx<<";"<<myBoxes[mynuc].ty<<";"<<myBoxes[mynuc].tz<<"|"<<largedefmsinfo.nx<<"-"<<largedefmsinfo.ny<<"-"<<largedefmsinfo.nz<<"|"<<gx_pbc<<";"<<gy_pbc<<";"<<gz_pbc<<"|"<< locdefms <<"|"<<myBoxes[mynuc].ngrx<<";"<<myBoxes[mynuc].ngry<<";"<<myBoxes[mynuc].ngrz << endl;
                        }
                    }
                } //structure constructed

		myDefMS[mynuc] = defms; //null was already pushed_back before

		//reference nucleus and its neighbors from largedefms
		double tentry = MPI_Wtime();

        init_myNuclei_nucleation_BasedOnLargeDefMS( mynuc, refnucleus );

		double texit = MPI_Wtime();
cout << mynuc << " time spent to identify neighbors = " << (texit - tentry) << endl;

        }
}


void polyxx::init_myNuclei ( void )
{
	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
		//determine first the relative displacement of the defgr aggregate to get the necessary number of defgrains which fully enclose the simulation volume

		myBoxes[mynuc].tx = -1.0 * rndlocal.leEcuyer() * defgmedian_rd;
		myBoxes[mynuc].ty = -1.0 * rndlocal.leEcuyer() * defgmedian_nd;
		myBoxes[mynuc].tz = -1.0 * rndlocal.leEcuyer() * defgmedian_td;

		//then how many grains in either direction are necessary?, tx, ty, tz are negative!
		myBoxes[mynuc].ngrx = (( boxx - myBoxes[mynuc].tx ) / defgmedian_rd); //int rounds off so...
		myBoxes[mynuc].ngrx++;
		if ( myBoxes[mynuc].ngrx == 0 ) myBoxes[mynuc].ngrx++;

		myBoxes[mynuc].ngry = (( boxy - myBoxes[mynuc].ty ) / defgmedian_nd);
		myBoxes[mynuc].ngry++;
		if ( myBoxes[mynuc].ngry == 0 ) myBoxes[mynuc].ngry++;

		myBoxes[mynuc].ngrz = (( boxz - myBoxes[mynuc].tz ) / defgmedian_td);
		myBoxes[mynuc].ngrz++;
		if ( myBoxes[mynuc].ngrz == 0 ) myBoxes[mynuc].ngrz++;

		myBoxes[mynuc].ngrxy = myBoxes[mynuc].ngrx * myBoxes[mynuc].ngry;
		myBoxes[mynuc].ngrxyz = myBoxes[mynuc].ngrxy * myBoxes[mynuc].ngrz;

		//allocate memory
		long n = myBoxes[mynuc].ngrxyz;
		long* defms = NULL;
		defms = new long[n];
		if (defms == NULL) { cout << myRank << "ERROR::Allocating memory during init_myNuclei!" << endl; }
		this->myMemGuard += n * sizeof(long);

//cout << myRank << ";" << mynuc << ";txyz;" << myBoxes[mynuc].tx << ";" << myBoxes[mynuc].ty << ";" << myBoxes[mynuc].tz << ";ngrxyz;" << myBoxes[mynuc].ngrx << ";" << myBoxes[mynuc].ngry << ";" << myBoxes[mynuc].ngrz << "\tn=" << n << endl;

		//alignment of the deformed grains is as in +x aligned, stacked in +y into slabs stacked in +z, so implicit 3D linear array

		//###MK::assign random uncorrelated placement of GIA grains
		for (long g = 0; g < n; g++) {
			defms[g] = rndlocal.leEcuyer() * ndefgpool; //[0,1) assures g to be smaller than ndefgpool
		}

		myDefMS[mynuc] = defms; //null was already pushed_back before
	}

        //place the nuclei and neighbors
        init_myNuclei_nucleation();
        
	/*
	//##VERBOSE but sloppy and inefficient...
	MPI_Barrier(MPI_COMM_WORLD);
	for (int r = MASTER; r < this->nRanks; r++) {

		if ( r == this->myRank ) {
			//##MK::GIA-suffices that describes aggregate of cuboid grains which can have different but locally homogeneous properties = ori, rho
			stringstream log_defms_fname;
			ofstream log_defms_file;
			log_defms_fname << "RXPATHTRACER.DefMS." << this->simulationid << ".DefMS." << myRank << ".csv";
			log_defms_file.open ( log_defms_fname.str().c_str() );
			//log_defms_file << "Rank(Seed);ngrx;ngry;ngrz;ngrxy;ngrxyz;tx;ty;tz;xrd;ynd;ztd;AllGrainID" << endl;
			log_defms_file << "Rank;Seed;VolA;VolB;VolC;TotalVol" << endl;

			double totalvol = 0.0;
			double volabc[3] = {0.0, 0.0, 0.0};
			log_defms_file << this->myRank << ";" << (((long) -1 * this->myRank)-1); // << ";" << myBoxes[mn].ngrx << ";" << myBoxes[mn].ngry << ";" << myBoxes[mn].ngrz << ";" << myBoxes[mn].ngrxy << ";" << myBoxes[mn].ngrxyz << ";" << myBoxes[mn].xrd << ";" << myBoxes[mn].ynd << ";" << myBoxes[mn].ztd << ";" << myBoxes[mn].tx << ";" << myBoxes[mn].ty << ";" << myBoxes[mn].tz;

			for ( long mn = 0; mn < myIDs.size(); mn++) {

				for (long z = 0; z < myBoxes[mn].ngrz; z++) {
					for (long y = 0; y < myBoxes[mn].ngry; y++) {
						for (long x = 0; x < myBoxes[mn].ngrx; x++) {

							long id = x + y * (myBoxes[mn].ngrx) + (z * myBoxes[mn].ngrx * myBoxes[mn].ngry);

							if ( id >= myBoxes[mn].ngrxyz )	cout << "ERROR::Access violation in mn " << mn << endl; 

							double xmi = x;		xmi *= defgmedian_rd;	xmi += myBoxes[mn].tx;
							double xmx = x+1;	xmx *= defgmedian_rd;	xmx += myBoxes[mn].tx;
							double ymi = y;		ymi *= defgmedian_nd;	ymi += myBoxes[mn].ty;
							double ymx = y+1;	ymx *= defgmedian_nd;	ymx += myBoxes[mn].ty;
							double zmi = z;		zmi *= defgmedian_td;	zmi += myBoxes[mn].tz;
							double zmx = z+1;	zmx *= defgmedian_td;	zmx += myBoxes[mn].tz;
							if ( xmi < 0.0 )				xmi = 0.0;
							if ( xmx > myBoxes[mn].xrd )	xmx = myBoxes[mn].xrd;
							if ( ymi < 0.0 )				ymi = 0.0;
							if ( ymx > myBoxes[mn].ynd )	ymx = myBoxes[mn].ynd;
							if ( zmi < 0.0 )				zmi = 0.0;
							if ( zmx > myBoxes[mn].ztd )	zmx = myBoxes[mn].ztd;

							double volgr = (xmx-xmi)*(ymx-ymi)*(zmx-zmi);
							if ( myDefMS[mn][id] == 0 ) volabc[0] = volabc[0] + volgr;
							if ( myDefMS[mn][id] == 1 ) volabc[1] = volabc[1] + volgr;
							if ( myDefMS[mn][id] == 2 ) volabc[2] = volabc[2] + volgr;
							totalvol = totalvol + volgr;
						}
					}
				}
				//for (unsigned long g = 0; g < myBoxes[mn].ngrxyz; g++) 	log_defms_file << ";" << myDefMS[mn][g];

			} //all boxes
			log_defms_file << ";" << (volabc[0]) << ";" << (volabc[1]) << ";" << (volabc[2]) << ";" << (totalvol) << endl;
			cout << myRank << "\t\t" << totalvol << endl;

			log_defms_file.flush();
			log_defms_file.close();
		} //worker ioed

		MPI_Barrier(MPI_COMM_WORLD);
	}
	//##VERBOSE
	*/
}


void polyxx::init_myNeighbors ( long towhichnuc ) 
{
	vector<nucsite> localenvironment;

	double tmpx = this->myNuclei[towhichnuc].x;
	double tmpy = this->myNuclei[towhichnuc].y;
	double tmpz = this->myNuclei[towhichnuc].z;
	double bx = this->myBoxes[towhichnuc].xrd;
	double by = this->myBoxes[towhichnuc].ynd;
	double bz = this->myBoxes[towhichnuc].ztd;

#ifdef DETAILED_PROMPTS
	cout << myRank << ";myRank;tmpx;tmpy;tmpz;" << tmpx << ";" << tmpy << ";" << tmpz << endl;
#endif
        //#############requires changes anyway!

	for (long nbor = 0; nbor < 125; nbor++) {
		//##MK::first show, ASSUMING CSR OF THE NEIGHBORS NUCLEATION SITES, equals mynuc realizations of CSR point processes with nRanks seeds on the parkMiller generator

		//##MK::change here to cut an arbitrary section from a R3 process that is predefined
		double ddd = 0.0;
		double xx, yy, zz;
		
		do { //search for a disjoint site EPSILONBALL farther apart from tmpx,y,z the testnucleus in simulation towhichnuc
			xx = rndlocal.leEcuyer() * bx;
			yy = rndlocal.leEcuyer() * by;
			zz = rndlocal.leEcuyer() * bz;

			ddd = SQR((tmpx - xx)) + SQR((tmpy - yy)) + SQR((tmpz - zz));

			ddd = pow (ddd, 0.5);
			
		} while ( ddd <= EPSILONBALL );

		//neighbor found
		//##MK::more sophisticated neighbor site model necessary!
		struct nucsite neighborsite;

		//##MK:::neighbor position
		neighborsite.x = xx;
		neighborsite.y = yy;
		neighborsite.z = zz;

		//##MK::neighbor grain/character orientation
		//long rgid = nbor % nrxgpool;

		long rgid = rndlocal.leEcuyer() * nrxgpool;

		neighborsite.rxgid = rgid;

		//###neighbor incubation time - debug - site-saturated
		neighborsite.tincub = SITE_SATURATED;

//#ifdef DETAILED_PROMPTS
	//cout << myRank << ";myRank;nbor;xx;yy;zz;ddd=" << nbor << ";" << xx << ";" << yy << ";" << zz << ";" << ddd << endl;
//#endif

		localenvironment.push_back( neighborsite ); //copy constructor
	}

	//if ( myNeighbors.size() != towhichnuc ) { cout << myRank << " ERROR:Inconsistent construction of neighborhood in init_myNeighbor" << endl; return; }

	//##MK::vector of vector not most efficient vector of pointer to struct, however paths are always local and as such the cache miss overhead amortized at least to some extend...
	this->myNeighbors.push_back( localenvironment ); //MK::thus the function MUST NOT HAVE AN ARGUMENT

//cout << myRank << " myNeighbors[mynuc].size = " << myNeighbors[towhichnuc].size() << endl;
}


void polyxx::init_myNeighbors_BasedOnPointProcess ( long towhichnuc ) 
{
	vector<nucsite> localenvironment;

	double tmpx = this->myNuclei[towhichnuc].x;
	double tmpy = this->myNuclei[towhichnuc].y;
	double tmpz = this->myNuclei[towhichnuc].z;
	double bx = this->myBoxes[towhichnuc].xrd;
	double by = this->myBoxes[towhichnuc].ynd;
	double bz = this->myBoxes[towhichnuc].ztd;

#ifdef DETAILED_PROMPTS
	cout << myRank << ";myRank;tmpx;tmpy;tmpz;" << tmpx << ";" << tmpy << ";" << tmpz << endl;
#endif
	//pick one location from the pointprocess about which a box of dimensions [-bi;+bi] does not protrude beyond the window
	//the pointprocess lives on [0, ppscaling]^3 micron

	long npointstotal = pointprocess.size();
	long domaincenter;
	double px, py, pz;
	double dlim = ppscaling;
	bool inside = false;

	while ( inside == false ) {
		//pick a location at random
		domaincenter = rndlocal.leEcuyer() * npointstotal;

		px = pointprocess[domaincenter].x;
		py = pointprocess[domaincenter].y;
		pz = pointprocess[domaincenter].z;

		//dismiss if protruding beyond a boundary of the pointprocess universe [0, ppscaling]^3
		if ( px < (bx/2.0) )			continue;
		if ( px > (dlim - (bx/2.0)) )	continue;
		if ( py < (by/2.0) )			continue;
		if ( py > (dlim - (by/2.0)) )	continue;
		if ( pz < (bz/2.0) )			continue;
		if ( pz > (dlim - (bz/2.0)) )	continue;

		//inside
		inside = true;
	}

	//establish limits in the coordinate system of the point process universe
	//##MK::sequential optimization divide a priori by 2...
	double xmi = px - (bx/2.0);
	double xmx = px + (bx/2.0);
	double ymi = py - (by/2.0);
	double ymx = py + (by/2.0);
	double zmi = pz - (bz/2.0);
	double zmx = pz + (bz/2.0);
	double ddd = 0.0;

	//now draw all points which are in this window
	long nbor = 0;
	for ( long lower = 0; lower < domaincenter; lower++ ) {
		if ( pointprocess[lower].x < xmi )	continue;
		if ( pointprocess[lower].x > xmx )	continue;
		if ( pointprocess[lower].y < ymi )	continue;
		if ( pointprocess[lower].y > ymx )	continue;
		if ( pointprocess[lower].z < zmi )	continue;
		if ( pointprocess[lower].z > zmx )	continue;

		//not continued, so inside!
		//##MK::more sophisticated neighbor site model necessary!
		struct nucsite neighborsite;

		//##MK:::neighbor position, affine transformation from pointprocess coorsys in box coorsys
		neighborsite.x = pointprocess[lower].x - xmi;
		neighborsite.y = pointprocess[lower].y - ymi;
		neighborsite.z = pointprocess[lower].z - zmi;

		//warning if distance to center grain is less than EPSILONBALL
		ddd = SQR((tmpx - neighborsite.x)) + SQR((tmpy - neighborsite.y)) + SQR((tmpz - neighborsite.z));
		if ( ddd <= SQR(EPSILONBALL) ) cout << "WARNING::Nucleus " << lower << " which becomes " << nbor << " too close to center!" << endl;

		//##MK::neighbor grain/character orientation
		//long rgid = nbor % nrxgpool;

		long rgid = rndlocal.leEcuyer() * nrxgpool;

		neighborsite.rxgid = rgid;

		//###neighbor incubation time - debug - site-saturated
		neighborsite.tincub = SITE_SATURATED;

//#ifdef DETAILED_PROMPTS
//cout << myRank << "lower;nbor;x;y;z;nbx;nby;nbz;ddd = " << lower << ";" << nbor << "--" << pointprocess[lower].x << ";" << pointprocess[lower].y << ";" << pointprocess[lower].z << "--" << neighborsite.x << ";" << neighborsite.y << ";" << neighborsite.z << "--" << ddd << endl;
//#endif

		localenvironment.push_back( neighborsite ); //copy constructor

		nbor++;
	}

	//jump over domain center to avoid the infrequent event that the iterator is the domaincenter
	for ( long upper = domaincenter + 1; upper < npointstotal; upper++ ) {
		if ( pointprocess[upper].x < xmi )	continue;
		if ( pointprocess[upper].x > xmx )	continue;
		if ( pointprocess[upper].y < ymi )	continue;
		if ( pointprocess[upper].y > ymx )	continue;
		if ( pointprocess[upper].z < zmi )	continue;
		if ( pointprocess[upper].z > zmx )	continue;

		//not continued, so inside!
		//##MK::more sophisticated neighbor site model necessary!
		struct nucsite neighborsite;

		//##MK:::neighbor position, affine transformation from pointprocess coorsys in box coorsys
		neighborsite.x = pointprocess[upper].x - xmi;
		neighborsite.y = pointprocess[upper].y - ymi;
		neighborsite.z = pointprocess[upper].z - zmi;

		//warning if distance to center grain is less than EPSILONBALL
		ddd = SQR((tmpx - neighborsite.x)) + SQR((tmpy - neighborsite.y)) + SQR((tmpz - neighborsite.z));
		if ( ddd <= SQR(EPSILONBALL) ) cout << "WARNING::Nucleus " << upper << " which becomes " << nbor << " too close to center!" << endl;
		//##MK::neighbor grain/character orientation
		//long rgid = nbor % nrxgpool;

		long rgid = rndlocal.leEcuyer() * nrxgpool;

		neighborsite.rxgid = rgid;

		//###neighbor incubation time - debug - site-saturated
		neighborsite.tincub = SITE_SATURATED;

//#ifdef DETAILED_PROMPTS
//cout << myRank << "upper;nbor;x;y;z;nbx;nby;nbz;ddd = " << upper << ";" << nbor << "--" << pointprocess[upper].x << ";" << pointprocess[upper].y << ";" << pointprocess[upper].z << "--" << neighborsite.x << ";" << neighborsite.y << ";" << neighborsite.z << "--" << ddd << endl;
//#endif

		localenvironment.push_back( neighborsite ); //copy constructor

		nbor++;
	}

	//if ( myNeighbors.size() != towhichnuc ) { cout << myRank << " ERROR:Inconsistent construction of neighborhood in init_myNeighbor" << endl; return; }

	//##MK::vector of vector not most efficient vector of pointer to struct, however paths are always local and as such the cache miss overhead amortized at least to some extend...
	this->myNeighbors.push_back( localenvironment ); //MK::thus the function MUST NOT HAVE AN ARGUMENT

//cout << myRank << " myNeighbors[mynuc].size = " << myNeighbors[towhichnuc].size() << endl;
}


void polyxx::init_myNeighbors_BasedOnLargeDefMS( long towhichnuc, long rid )
{
	vector<nucsite> localenvironment;

	double bcx = largedefmsnuc[rid].cx; //this->myNuclei[towhichnuc].x;
	double bcy = largedefmsnuc[rid].cy; //this->myNuclei[towhichnuc].y;
	double bcz = largedefmsnuc[rid].cz; //this->myNuclei[towhichnuc].z;
        //box with centers bcx, bcy, bcz, and dimensions bx, by, bz is not protruding outside largedefms structure hence find all neighbors excluding rid inside the box
	double bx = this->myBoxes[towhichnuc].xrd;
	double by = this->myBoxes[towhichnuc].ynd;
	double bz = this->myBoxes[towhichnuc].ztd;

#ifdef DETAILED_PROMPTS
	##cout << myRank << ";myRank;tmpx;tmpy;tmpz;" << tmpx << ";" << tmpy << ";" << tmpz << endl;
#endif
	//##MK::find all nuclei that are not rid and inside the box, naive approach here...
        
	double xmi = bcx - (bx/2.0);
	double xmx = bcx + (bx/2.0);
	double ymi = bcy - (by/2.0);
	double ymx = bcy + (by/2.0);
	double zmi = bcz - (bz/2.0);
	double zmx = bcz + (bz/2.0);
	double ddd = 0.0;

	//now draw all points which are in this window
        long npointstotal = largedefmsnuc.size();
        
	long nbor = 0;
	for ( long lower = 0; lower < rid; lower++ ) {
		if ( largedefmsnuc[lower].cx < xmi )	continue;
		if ( largedefmsnuc[lower].cx > xmx )	continue;
		if ( largedefmsnuc[lower].cy < ymi )	continue;
		if ( largedefmsnuc[lower].cy > ymx )	continue;
		if ( largedefmsnuc[lower].cz < zmi )	continue;
		if ( largedefmsnuc[lower].cz > zmx )	continue;

		//not continued, so inside!
		//##MK::more sophisticated neighbor site model necessary!
		struct nucsite neighborsite;

		//##MK:::neighbor position, affine transformation from pointprocess coorsys in box coorsys
		neighborsite.x = largedefmsnuc[lower].cx - xmi;
		neighborsite.y = largedefmsnuc[lower].cy - ymi;
		neighborsite.z = largedefmsnuc[lower].cz - zmi;

		//warning if distance to center grain is less than EPSILONBALL
		ddd = SQR((bcx - neighborsite.x)) + SQR((bcy - neighborsite.y)) + SQR((bcz - neighborsite.z));
		if ( ddd <= SQR(EPSILONBALL) ) cout << "WARNING::Nucleus " << lower << " which becomes " << nbor << " too close to center!" << endl;

		//##MK::neighbor grain/character orientation
		//long rgid = nbor % nrxgpool;
		//long rgid = rndlocal.leEcuyer() * nrxgpool;

		neighborsite.rxgid = largedefmsnuc[lower].rxgpoolid;

		//##MK::neighbor incubation time - debug - site-saturated
		neighborsite.tincub = SITE_SATURATED;

//#ifdef DETAILED_PROMPTS
//cout << myRank << "lower;nbor;x;y;z;nbx;nby;nbz;ddd = " << lower << ";" << nbor << "--" << pointprocess[lower].x << ";" << pointprocess[lower].y << ";" << pointprocess[lower].z << "--" << neighborsite.x << ";" << neighborsite.y << ";" << neighborsite.z << "--" << ddd << endl;
//#endif

		localenvironment.push_back( neighborsite ); //copy constructor

		nbor++;
	}

	//jump over domain center to avoid the infrequent event that the iterator is the domaincenter
	for ( long upper = rid + 1; upper < npointstotal; upper++ ) {
		if ( largedefmsnuc[upper].cx < xmi )	continue;
		if ( largedefmsnuc[upper].cx > xmx )	continue;
		if ( largedefmsnuc[upper].cy < ymi )	continue;
		if ( largedefmsnuc[upper].cy > ymx )	continue;
		if ( largedefmsnuc[upper].cz < zmi )	continue;
		if ( largedefmsnuc[upper].cz > zmx )	continue;

		//not continued, so inside!
		//##MK::more sophisticated neighbor site model necessary!
		struct nucsite neighborsite;

		//##MK:::neighbor position, affine transformation from pointprocess coorsys in box coorsys
		neighborsite.x = largedefmsnuc[upper].cx - xmi;
		neighborsite.y = largedefmsnuc[upper].cy - ymi;
		neighborsite.z = largedefmsnuc[upper].cz - zmi;

		//warning if distance to center grain is less than EPSILONBALL
		ddd = SQR((bcx - neighborsite.x)) + SQR((bcy - neighborsite.y)) + SQR((bcz - neighborsite.z));
		if ( ddd <= SQR(EPSILONBALL) ) cout << "WARNING::Nucleus " << upper << " which becomes " << nbor << " too close to center!" << endl;
		//##MK::neighbor grain/character orientation
		//long rgid = nbor % nrxgpool;

		//long rgid = rndlocal.leEcuyer() * nrxgpool;
		//neighborsite.rxgid = rgid;
                
                neighborsite.rxgid = largedefmsnuc[upper].rxgpoolid;

		//###neighbor incubation time - debug - site-saturated
		neighborsite.tincub = SITE_SATURATED;

//#ifdef DETAILED_PROMPTS
//cout << myRank << "upper;nbor;x;y;z;nbx;nby;nbz;ddd = " << upper << ";" << nbor << "--" << pointprocess[upper].x << ";" << pointprocess[upper].y << ";" << pointprocess[upper].z << "--" << neighborsite.x << ";" << neighborsite.y << ";" << neighborsite.z << "--" << ddd << endl;
//#endif

		localenvironment.push_back( neighborsite ); //copy constructor

		nbor++;
	}

	//##MK::vector of vector not most efficient vector of pointer to struct, however paths are always local and as such the cache miss overhead amortized at least to some extend...
	this->myNeighbors.push_back( localenvironment ); //MK::thus the function MUST NOT HAVE AN ARGUMENT

cout << "Neighbors to nucleus " << towhichnuc << " identified  a total of = " << localenvironment.size() << endl;
//cout << myRank << " myNeighbors[mynuc].size = " << myNeighbors[towhichnuc].size() << endl;    
}


void polyxx::init_myNuclei_nucleation ( void )
{
	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
		struct nucsite ncs;

		//##MK::nucleus is by default placed in the center of the simulation domain
		ncs.x = 0.5 * myBoxes[mynuc].xrd;
		ncs.y = 0.5 * myBoxes[mynuc].ynd;
		ncs.z = 0.5 * myBoxes[mynuc].ztd;

		//##MK::incubation time model - debug - site saturated
		ncs.tincub = SITE_SATURATED;

		//##MK::nucleation texture - debug - random, //##MK::possibly MPI is necessary to pass particular nucleus orientations to some nodes
		ncs.rxgid = mynuc % nrxgpool;

		//ncs.rxgid = rndlocal.leEcuyer() * nrxgpool; //[0,1) guarantees smaller than nrxgpool

		myNuclei.push_back( ncs ); //copy-constructor

		//MK::now place appropriate neighbors around the nucleus - old method
		//init_myNeighbors( mynuc );
		//MK::now place appropriate neighbors around the nucleus - new method
		init_myNeighbors_BasedOnPointProcess( mynuc );
	}
	cout << "RANK: " << myRank << " selected Nuclei = " << myNuclei.size() << endl;
}


void polyxx::init_myNuclei_nucleation_BasedOnLargeDefMS( long towhichnuc, long refid )
{
		struct nucsite ncs;
         
                //##MK::nucleus is by default placed in the center of the simulation domain
  		ncs.x = 0.5 * myBoxes[towhichnuc].xrd;
		ncs.y = 0.5 * myBoxes[towhichnuc].ynd;
		ncs.z = 0.5 * myBoxes[towhichnuc].ztd;

		//##MK::incubation time model - debug - site saturated
		ncs.tincub = SITE_SATURATED;

                //all nucleus orientations were already registered
                //identify the nucleus in the box center from the large microstructure largdefmsnuc
                ncs.rxgid = largedefmsnuc[refid].rxgpoolid;
                
                //ncs.rxgid = rndlocal.leEcuyer() * nrxgpool; //[0,1) guarantees smaller than nrxgpool

		myNuclei.push_back( ncs ); //copy-constructor

		//MK::now place appropriate neighbors around the nucleus - old method
		//init_myNeighbors( mynuc );
		//MK::now place appropriate neighbors around the nucleus - new method
		//init_myNeighbors_BasedOnPointProcess( mynuc );
                init_myNeighbors_BasedOnLargeDefMS( towhichnuc, refid );

                cout << "RANK: " << myRank << " selected Nuclei = " << myNuclei.size() << endl;       
}


void polyxx::init_disori_lowertriangle ( void )
{
	//###further optimization, do this in parallel, currently all nodes know all orientations, thus all possible disorientations

    double timer = MPI_Wtime();
#ifdef DETAILED_PROMPTS
	cout << "Node " << this->myRank << " " << oripool.size() << " disjunct orientations are currently in the orientation pool, counting noripool " << noripool << endl;
#endif

	//double q1[4], q2[4];calculate disorientation between all components
	double q1[4], q2[4], qdis[4];

	double weightedMob;
	double maxDev40_111 = 0.174532925199433; //10 / 180 * _PI_;
	double _sqrt3 = 1.0 / sqrt( 3.0 );
	double oneNinth = 1.0 / 9.0; 
	//the quaternion that describes a 40deg<111> misorientation
	double m40_111[4] = { cos( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ) };

	//set up lower-triangle matrix
	long nentries = 0.5 * oripool.size() * (oripool.size() + 1);

cout << "Disori matrix has = " << nentries << endl;

	//set a realistic threshold above which disorientations are not precached but rather failure occurs suggesting to 
	//caching is beneficial if disorientation is often to be dereferenced, but overall size of matrix stays in reasonable memory limit
	if ( (( nentries * sizeof(double) + sizeof(double*)) / CUBE(1024)) > DISORICACHE_MEMLIMIT ) { cout << "ERROR::Disorientation lower-triangle matrix too large." << endl; return; }


	long entry;
	LowTrigMobWeightHash = NULL;
	LowTrigMobWeightHash = new double[nentries];
	if (LowTrigMobWeightHash == NULL) { cout << "ERROR::During allocation of memory for lower triangle matrix in init_disori" << endl; return; }
	this->myMemGuard += nentries * sizeof(double);

	//###MK::debug - parallelize, otherwise a potential bottleneck!
	for (long i = 0; i < nentries; i++) {
		LowTrigMobWeightHash[i] = 0.0;
	}
	//cout<<"Ori-ID | Ori [phi1,PHI,phi2]"<<endl;
	/*
	for (long row = 0; row < oripool.size(); row++){
		cout<<row<<" | ["<< oripool[row].bunge1*RADTODEG <<" , " << oripool[row].bunge2*RADTODEG <<" , "<< oripool[row].bunge3*RADTODEG<<" ]"<<endl;
	}
	*/
	double Pmax = -2.0;
	//cout<<"Mobility Lower Triag Matrix Output:"<<endl;
	//cout<<"___________________________________"<<endl;
	for (long row = 0; row < oripool.size(); row++ ) {
		q1[0] = oripool[row].q0;
		q1[1] = oripool[row].q1;
		q1[2] = oripool[row].q2;
		q1[3] = oripool[row].q3;

		for (long col = 0; col <= row; col++) {
			//row >= col is assured, here assuming switching symmetry!
			q2[0] = oripool[col].q0;
			q2[1] = oripool[col].q1;
			q2[2] = oripool[col].q2;
			q2[3] = oripool[col].q3;

			weightedMob = 0;

			misorientationQuaternionCubic( q1, q2, qdis );

			double theta = qdis[0];
			if( theta > 1.0 ) {
				theta = (double) (int) theta;
			}
			theta = 2*acos(theta);

			entry = (( 0.5 * row * (row + 1) ) + col);

			if ( mobilitymodel_option == SEBALD_GOTTSTEIN ) {
				//assume LAGB
				if ( theta <= 0.261799387799149 ) { //15deg/180*PI
					weightedMob = -1.0;
					LowTrigMobWeightHash[entry] = weightedMob;
					//cout<<weightedMob<<"(LAGB)\t";
					continue;
				}

				//okay obviously HAGB
				weightedMob = 0.0;

				//close to 40deg<111> boundary?

				//##MK::even though it works I do not know why in particular it is possible to apply the Grimmer disorientation approach to two misorientation quaternions to get dev_40_111
				double dev_40_111 = misorientationCubicQxQ( qdis[0], qdis[1], qdis[2], qdis[3], m40_111[0], m40_111[1], m40_111[2], m40_111[3] );

				//###MK::dev_40_111 is consistent with MTex3.3.1 as well as all others quaternions but 59.0, 36.7, 63.4 vs 73.3033, 67.1493, 78.1088 even though 40deg<111> rot is not dev40_111 = 0.0 !

				if( dev_40_111 <  maxDev40_111 ) {
					weightedMob = SQR( cos( 0.5 * _PI_ * dev_40_111 / maxDev40_111 ) );
				}

				LowTrigMobWeightHash[entry] = weightedMob;

				continue;
			}

			//not continued, alternative model ROLLETT_HUMPHREYS
			weightedMob = 1.0 - (RHLAGBHAGBCut * exp( -1.0 * RHLAGBHAGBTrans * pow( (theta/0.261799387799149), RHLAGBHAGBExponent ) ));
			LowTrigMobWeightHash[entry] = weightedMob;

//cout << "weigthedMob/theta/row/col/oprow/opcol------>" << weightedMob << ";" << theta << "--" << row << ";" << col << "--" << oripool[row].bunge1/_PI_*180.0 << ";" << oripool[col].bunge1/_PI_*180.0 << "--" << get_lowtrig_mobilityweight( row, col ) << endl;

			//double coor = get_lowtrig_mobilityweight( row , col );
			//cout<<weightedMob<<"\t";
		} //for all columns
		//cout<<endl;
	} //fill the row
	//cout<<"___________________________________"<<endl;
	if (myRank == MASTER) { cout << "Lower triangle matrix calculated in " << setprecision(6) << ((double) MPI_Wtime() - timer) << " seconds " << endl; }
}


void polyxx::init_checkTimeSteppingOptimPotential ( void ) 
{
	//MK/PH::not in all automata we do have all the same orientations and mobilities, this might allow the dynamical time stepping to be optimized separately in the nodes

	//principal design choice here was like that, we would like to use as many nodes as possible to be fast, 
	//thus in a first pure MPI implementation as few mynuc-nbors containers as possible, these form a local collection of nuclei that in order to simplify the postprocessing
	//of kinetics, texture and grain size run in the same time pattern.

	//##However this is problematic when strong heterogeneities suggest potential for significant local optimized time stepping but then requiring global rediscretization for all nuclei
	//as well as collecting, because each differently tuned single solution adds necessary counters for intermediate quants instead of being able to aggregate data
	//the apparent target conflict of minimizing the total number of interations per nucleus/path here obviously contradicts not a priori decidable with increasiong memory and thus rediscretization effort 

	//find min/max P-Value over all of my nuclei-deformed grains combinations, maxRho-Value similarly
	optTimeStepping.myPmax = 0.0;
	optTimeStepping.myPmin = 1.0;
	optTimeStepping.myRhomax = 0.0;

	double Pcur, rho_cur;
	long rxgoid;

	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
		rxgoid = rxgpool[myNuclei[mynuc].rxgid].ori;

		for (long defg = 0; defg < myBoxes[mynuc].ngrxyz; defg++){
			Pcur = get_lowtrig_mobilityweight( rxgoid, defgpool[myDefMS[mynuc][defg]].ori );

			rho_cur = get_rho( myDefMS[mynuc][defg] ); //##MK::static, without recovery

			if( Pcur > optTimeStepping.myPmax ) 
				optTimeStepping.myPmax = Pcur;
			if( Pcur > -0.5 && Pcur < optTimeStepping.myPmin) 
				optTimeStepping.myPmin = Pcur; //PH::if Pcur=-1 (LAGB) don't consider as Pmin
			if(rho_cur > optTimeStepping.myRhomax) 
				optTimeStepping.myRhomax = rho_cur;
		}
	}
	cout <<"RANK: " << myRank << " myRhomax = " << optTimeStepping.myRhomax << " myPmax = " << optTimeStepping.myPmax << " myPmin = " << optTimeStepping.myPmin << endl;
}


void polyxx::init_MyTimeSteps ( void ) 
{
	//##THIS WILL NOT WORK AS SOON THERE IS A DISTRIBUTION OF NUCLEI!?
	double tincub_min = myNuclei[0].tincub; //init
	double tincub_cur;
	long nbor;
	double T;
	//double r_mean = 2 * pow( ( boxx*boxy*boxz*3 / ((nNeighborsPerCA+1)*4*_PI_)) , ONETHIRD ); //MKPH::average distance

	double r_mean = (1.0 / this->ppdensity) / (4.0 * _PI_) * 3.0; //ppdensity is now in micron^3
        r_mean = 2 * pow( r_mean, ONETHIRD );
cout << "init_MyTimeSteps with r_mean = " << r_mean << " micron precached." << endl;

        myTimeAlloc = ( (long) TSTEPMINALLOC * r_mean / minDist ) + 1; //MK::added +1, ##implicit cast suffices?

cout << "init_MyTimeSteps with r_mean/myTimeAlloc = " << r_mean << " micron precached with " << myTimeAlloc << " preallocated slots." << endl;	
        
	//myTimeSteps.reserve( myTimeAlloc ); //PH::one page 4KB
	//myTemperatures.reserve( myTimeAlloc );
	
	//find min tincub
	for ( long mynuc = 0; mynuc < myIDs.size(); mynuc++){
		if( myNuclei[mynuc].tincub < tincub_min ) 
			tincub_min = myNuclei[mynuc].tincub;
                
		for (nbor = 0; nbor < myNeighbors[mynuc].size(); nbor++ ) {
			tincub_cur = myNeighbors[mynuc][nbor].tincub;
			if( tincub_cur < tincub_min ) 
				tincub_min = tincub_cur;
		}
	}

	//PB::set first processing point to (tincub_min, temperature T(tincub_min)) , this avoids scanning for each path and nucleus in the local ensemble the cache misses on evaluating the processing profile
	T = get_temperature( tincub_min );
	//populate with the first pair of data
	myTimeSteps.push_back( tincub_min );
	myTemperatures.push_back(T);
	
	//calculate timestep-sizes dependent on max velocity 
	init_precomputeMyTimeTemperature(tincub_min, (myTimeAlloc - 1) ); //-1 first already placed
}



void polyxx::init_precomputeMyTimeTemperature( double t_0, long ntsteps )
{
	double t_cur = t_0; //continue calculating from t zero

	double T_cur, neg_kT;
	double v1, v2, vLAGB, vmax;

	T_cur = get_temperature( t_cur ); //always delivers a temperature if time is beyond processing scheme, it delivers the last time or the Tiso if quick isothermal

	for ( long j = 0; j < ntsteps; j++ ) {
		neg_kT = (NEG_KBOLTZMANN / T_cur);
		physConstants.update_Constants(T_cur);
                
		if ( mobilitymodel_option == SEBALD_GOTTSTEIN ) {
			mGS = GSm0 * exp( GSHact * neg_kT );
			mLAGB = LAGBm0 * exp( LAGBHact * neg_kT );
			mHAGB = HAGBm0 * exp( HAGBHact * neg_kT );
			//find out max possible velocity for this timestep, condensation effect is considered
			vLAGB = mLAGB * physConstants.get_halfG_b2() * optTimeStepping.myRhomax;

			vmax = vLAGB;
			v1 = (optTimeStepping.myPmin * mGS + (1 - optTimeStepping.myPmin) * mHAGB) * physConstants.get_halfG_b2() * optTimeStepping.myRhomax;
			if (v1 > vmax) 
				vmax = v1;
			v2 = (optTimeStepping.myPmax * mGS + (1 - optTimeStepping.myPmax) * mHAGB) * physConstants.get_halfG_b2() * optTimeStepping.myRhomax;
			if (v2 > vmax) 
				vmax = v2;
		}
		if ( mobilitymodel_option == ROLLETT_HUMPHREYS ) {
                    mRHHAGB = RHHAGBm0 * exp( RHHAGBHact * neg_kT );
                    vmax = mRHHAGB * physConstants.get_halfG_b2() * optTimeStepping.myRhomax;
		}

		//set timestep-size according to max velocity and minimal distance
		t_cur += minDist / vmax;
		T_cur = get_temperature( t_cur );

//cout << "tcur/Tcur/minDist/mGS/mRHHAGB/vmax/G/b/Gb2/rhoMax = " << t_cur << ";" << T_cur << ";" << minDist << ";" << mGS << ";" << mRHHAGB << ";" << vmax << ";" << physConstants.get_G() << ";" << physConstants.get_b() << ";" << physConstants.get_halfG_b2() << ";" << optTimeStepping.myRhomax << endl;

		myTimeSteps.push_back( t_cur );
		myTemperatures.push_back( T_cur );
	}
}


void polyxx::init_disori_fast ( void )
{
	//###MK::utilize mathMethods to categorize orientations, all nodes know all orientations, thus all possible disorientations
	//MK::but not all MPI nodes need to calculate all values, if the number of disjoint oris in the oripool is small best thing to to n x n via MPI or OpenMP, if not each node does his own subset to save memory and time!

#ifdef DETAILED_PROMPTS
	cout << "Node " << this->myRank << " " << oripool.size() << " disjunct orientations are currently in the orientation pool, counting noripool " << noripool << endl;
#endif

	//double q1[4], q2[4];calculate disorientation between all my nuclei orientations an all deformation components
	double q1[4], q2[4], qdis[4];

	double weightedMob;
	double maxDev40_111 = 0.174532925199433; //10 / 180 * _PI_;
	double _sqrt3 = 1 / sqrt( 3.0 );
	double oneNinth = 1.0 / 9.0; 
	//the quaternion that describes a 40Â°<111> misorientation
	double m40_111[4] = { cos( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ), _sqrt3 * sin( oneNinth * _PI_ ) };

	///allocate an array of size myCAsomany * oripool_firstdisjunct_rxg doubles for mobility weights
	long nentries = myIDs.size() * oripool_firstdisjunct_rxg;

cout << "myRank/nentries for debug::" << myRank << "\t\t" << nentries << endl;

	myCA_MobWeightHash = NULL;
	myCA_MobWeightHash = new double[nentries];
	if (myCA_MobWeightHash == NULL) { cout << myRank << "ERROR::During allocation of memory for myCA_MobWeightHash in init_disori_fast" << endl; return; }
	this->myMemGuard += nentries * sizeof(double);
	//###set entire array to zero not really necessary but good approach to detect potential failed P values then the simulation considers all boundaries to be high angle!
	for (long i = 0; i < nentries; i++) { myCA_MobWeightHash[i] = 0.0; }

	double local_Pmax = -2.0;

	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
		long nuc = myNuclei[mynuc].rxgid; //get id to oripool of the mynuc-th nucleus

		//data in the MobWeightHash is organized implicit with stripes of length oripool_firstdisjunct_rxg...
		long offset = mynuc * oripool_firstdisjunct_rxg;

		q1[0] = oripool[nuc].q0;
		q1[1] = oripool[nuc].q1;
		q1[2] = oripool[nuc].q2;
		q1[3] = oripool[nuc].q3;

		//cache all potential misorientation scenarios that could be encountered
		//MK::201402PARA potential for even higher reduction in memory consumption when hashtab is created on the fly...
		for ( long dgid = 0; dgid < oripool_firstdisjunct_rxg; dgid++ ) {
			q2[0] = oripool[dgid].q0;
			q2[1] = oripool[dgid].q1;
			q2[2] = oripool[dgid].q2;
			q2[3] = oripool[dgid].q3;

			weightedMob = 0.0;

			misorientationQuaternionCubic( q1, q2, qdis );

			double theta = qdis[0];
			if( theta > 1.0 ) { theta = (double) (int) theta; }
			theta = 2*acos(theta);

			//assign categorical intrinsic boundary mobility to disorientation among the two orientations
			if( theta <= 0.261799387799149 ) //15deg
				weightedMob = -1.0;
			else {
				weightedMob = 0.0;

				//how close is the orientation to a 40deg<111> in quaternion space?
				double dev_40_111 = misorientationCubicQxQ( qdis[0], qdis[1], qdis[2], qdis[3], m40_111[0], m40_111[1], m40_111[2], m40_111[3] );

				//oh close enough, apply Sebald-Gottstein-type-model mobility attentuation
				if( dev_40_111 <  maxDev40_111 ) {
					weightedMob = SQR( cos( 0.5 * _PI_ * dev_40_111 / maxDev40_111 ) );
				}
			}

			myCA_MobWeightHash[offset + dgid] = weightedMob;
			//if (weightedMob >= local_Pmax) { local_Pmax = weightedMob; }


cout << "MYNUC=" << mynuc << "\tDGID=" << dgid << "\tP\t" << weightedMob << "\toffset\t" << offset << "\tentry\t" << (offset + dgid) << "\tmyCA_MobWeightHash[entry]\t" << myCA_MobWeightHash[offset + dgid] << endl;
cout << q1[0] << ";" << q1[1] << ";" << q1[2] << ";" << q1[3] << "\t" << q2[0] << ";" << q2[1] << ";" << q2[2] << ";" << q2[3] << endl;
		} //for all disjunct defg orientatons
	} //for all myNuclei

	if (myRank == MASTER) {cout << "A myCA_MobWeightHashed, maximum local P found in my Node/me is = " << local_Pmax << endl;}
}


void polyxx::sim_pathsegmentation_nbors( pathP pl, long mn ) //, long nb, double nxx, double nyy, double nzz )
{
	vector<double> dtilde;
	dtilde.reserve ( myBoxes[mn].ngrxyz ); //dtilde is auto-freed after function exit
	bool newentry = true;

	//find interceptions with x, y, z planes by calculating dtilde of the path parameterization (ox,oy,oz)**T + dtilde * (u,v,w)**T
	double plen = pl->totallength;

	//tx,ty,tz shift of the grain structure is always < 1.0, ox, oy, oz is the origin of the path!
	//XPlanes first
	double uu = pl->direcx;
	if ( fabs(uu) > EPS_INPLANE_ACCURACY ) { //fabs to avoid /= 0.0 and to catch inplane alignment
		long nplxmax = myBoxes[mn].ngrx;

		for ( long nplx = 1; nplx <= nplxmax; nplx++ ) { //only if nplx strictly smaller so no boundary at only one grain in x direction, MK::however will be catched by < plen anyway
			double dtt = (defgmedian_rd * ((double) nplx) ) + myBoxes[mn].tx; //+ because txyz < 0, oxyz > 0
			dtt = dtt - pl->ox;
			dtt /= uu;

			//only when dtt is positive, < pl->totallength, and the value to EPS_SEGMENTS_ELIMDOUBLE accuracy not already there is a valid interception, otherwise wrong direction (dtt < 0), plane cutting at ox ( dtt == 0.0 ), outside the box
			newentry = true;
			for ( unsigned int i = 0; i < dtilde.size(); i++) { if ( (SQR(dtt - dtilde[i])) < SQR(EPS_SEGMENTS_ELIMDOUBLE) ) newentry = false; }

			//negative dtt cannot be accepted as we consider only the path as a path segment starting at nucsite and pointing in the direction of the nbors in the local coordinate system
			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); } //MK::could also be made <= plen
		} //all planes
	} //else, dtilde is still left empty if path lays in the yz plane

	//YPlanes
	double vv = pl->direcy;
	if ( fabs(vv) > EPS_INPLANE_ACCURACY ) {
		long nplymax = myBoxes[mn].ngry;

		for ( long nply = 1; nply <= nplymax; nply++ ) {
			double dtt = (defgmedian_nd * ((double) nply) ) + myBoxes[mn].ty;
			dtt = dtt - pl->oy;
			dtt /= vv;

			newentry = true;
			for ( unsigned i = 0; i < dtilde.size(); i++) { if ( (SQR(dtt - dtilde[i])) < SQR(EPS_SEGMENTS_ELIMDOUBLE) ) newentry = false; }

			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); }
		}
	}

	//ZPlanes
	double ww = pl->direcz;
	if ( fabs(ww) > EPS_INPLANE_ACCURACY ) {
		long nplzmax = myBoxes[mn].ngrz;

		for ( long nplz = 1; nplz <= nplzmax; nplz++ ) {
			double dtt = (defgmedian_td * ((double) nplz)) + myBoxes[mn].tz;
			dtt = dtt - pl->oz;
			dtt /= ww; //now division is never by 0.0

			newentry = true;
			for ( unsigned i = 0; i < dtilde.size(); i++) { if ( SQR((dtt - dtilde[i])) < SQR(EPS_SEGMENTS_ELIMDOUBLE) ) newentry = false; }

			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); }
		}
	}

	//MK::sort this list in ascending order
	std::sort( dtilde.begin(), dtilde.end(), SortDblAscending );

	//do not forget to put the totallength, as for in particular necessary when the raytracing process did not find any interception wiut a defgr boundary
	dtilde.push_back( plen );

	//MK::elimination of reduced dtilde values for interceptions at intersection lines of planes is now obsolete 
	//because there are cant be any doubles as the trace is parameterized and dtilde the only free parameter but the trace
	//is confied along the direction direcx,y,z in the box, in the case that the ray cuts exactly in a corner of a grain
	//it intersects both with a pair xy, yz or xz und thus is exactly one time accounted for!

	//allocate memory for local information
	int nseg = dtilde.size();
	pl->nsegments = nseg;

	//##MK::in the case of the SPHEREMESH_MODEL the outbound path is not necessary...
	regionpairP defms = NULL;
	defms = new regionpair[nseg];
	if (defms == NULL) { cout << "ERROR::During memory allocation for defms path in rank,nuc,nbor " << myRank << ";" << mn << endl; return; }
	this->myMemGuard += nseg * sizeof(regionpair);

	//get the coordinates of the points halfway in each section to assign the ID of the corresponding grains 
	int ix, iy, iz, ixyz;

	//MK::only INBOUND, solver always iterates forward nbors towards nucleus!
	double dtraveled_inb;

	if ( dtilde.size() == 1 ) {
		dtraveled_inb = dtilde[0] * 0.5;
		ix = ( (pl->ox + (dtraveled_inb * uu) - myBoxes[mn].tx) / defgmedian_rd);
		iy = ( (pl->oy + (dtraveled_inb * vv) - myBoxes[mn].ty) / defgmedian_nd);
		iz = ( (pl->oz + (dtraveled_inb * ww) - myBoxes[mn].tz) / defgmedian_td); //############was ix! n00b
		ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);

		//cout << "ixyz;out;" << ixyz << "\t";

		defms[0].dgid = this->myDefMS[mn][ixyz];

	}
	if ( dtilde.size() > 1 ) {
		//add the first
		dtraveled_inb = dtilde[0] * 0.5;
		ix = ((pl->ox + (dtraveled_inb * uu) - myBoxes[mn].tx) / defgmedian_rd); //rounds off
		iy = ((pl->oy + (dtraveled_inb * vv) - myBoxes[mn].ty) / defgmedian_nd);
		iz = ((pl->oz + (dtraveled_inb * ww) - myBoxes[mn].tz) / defgmedian_td);
		ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);

		//cout << "ixyz;out;" << ixyz << "\t";

		defms[0].dgid = this->myDefMS[mn][ixyz];

		//if dtilde.size > 1 then dtilde.size is at least 2, so access to element id 1 is valid and in that case following for loop executed once
		for ( long inbound = 1; inbound < dtilde.size(); inbound++ ) {
			dtraveled_inb = dtilde[inbound-1] + ((dtilde[inbound] - dtilde[inbound-1]) * 0.5); //dtilde strictly positive and in ascending order at least ELIMDOUBLE apart
			ix = ((pl->ox + (dtraveled_inb * uu) - myBoxes[mn].tx) / defgmedian_rd); 
			iy = ((pl->oy + (dtraveled_inb * vv) - myBoxes[mn].ty) / defgmedian_nd); 
			iz = ((pl->oz + (dtraveled_inb * ww) - myBoxes[mn].tz) / defgmedian_td); 
			ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);

			//cout << "ixyz;out;" << ixyz << "\t";

			defms[inbound].dgid = this->myDefMS[mn][ixyz];
		}
	}

	//now, construct the regionends for the forward iterate for INBOUND (simple the dtilde values)
	for ( unsigned long inb = 0; inb < dtilde.size(); inb++ ) {
		defms[inb].regionends = dtilde[inb];
	}


	pl->grainregions = defms;

	dtilde.clear();
}


void polyxx::sim_pathsegmentation_nucleus( pathP pl, long mn )
{
	vector<double> dtilde;
	dtilde.reserve ( myBoxes[mn].ngrxyz ); //dtilde is auto-freed after function exit
	bool newentry = true;

	//find interceptions with x, y, z planes by calculating dtilde of the path parameterization (ox,oy,oz)**T + dtilde * (u,v,w)**T
	double plen = pl->totallength;

	//tx,ty,tz shift of the grain structure is always < 1.0
	//XPlanes first
	double uu = pl->direcx;
	if ( fabs(uu) > EPS_INPLANE_ACCURACY ) { //to prevent /= 0.0 and to capture in plane alignment
		long nplxmax = myBoxes[mn].ngrx;

		for ( long nplx = 1; nplx <= nplxmax; nplx++ ) { //only if nplx strictly smaller so no boundary at only one grain in x direction, MK::however will be catched by < plen anyway
			double dtt = (defgmedian_rd * ((double) nplx) ) + myBoxes[mn].tx; //+ because txyz < 0, oxyz > 0, MK::tx is in absolute values ie Micron not relative!
			dtt = dtt - pl->ox;
			dtt /= uu; 

			//only when dtt is positive, < pl->totallength, and the value to EPS_SEGMENTS_ELIMDOUBLE accuracy not already there is a valid interception, otherwise wrong direction (dtt < 0), plane cutting at ox ( dtt == 0.0 ), outside the box
			newentry = true;
			for ( unsigned int i = 0; i < dtilde.size(); i++) { if ( (SQR(dtt - dtilde[i])) < SQR(EPS_SEGMENTS_ELIMDOUBLE) ) newentry = false; }

			//negative dtt cannot be accepted as we consider only the path as a path segment starting at nucsite and pointing in the direction of the nbors in the local coordinate system
			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); } //MK::could also be made <= plen
		} //all planes
	} //else, dtilde is still left empty if path lays in the yz plane

	//YPlanes
	double vv = pl->direcy;
	if ( fabs(vv) > EPS_INPLANE_ACCURACY ) {
		long nplymax = myBoxes[mn].ngry;

		for ( long nply = 1; nply <= nplymax; nply++ ) {
			double dtt = (defgmedian_nd * ((double) nply) ) + myBoxes[mn].ty;
			dtt = dtt - pl->oy;
			dtt /= vv;

			newentry = true;
			for ( unsigned i = 0; i < dtilde.size(); i++) { if ( (SQR(dtt - dtilde[i])) < SQR(EPS_SEGMENTS_ELIMDOUBLE) ) newentry = false; }

			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); }
		}
	}

	//ZPlanes
	double ww = pl->direcz;
	if ( fabs(ww) > EPS_INPLANE_ACCURACY ) {
		long nplzmax = myBoxes[mn].ngrz;

		for ( long nplz = 1; nplz <= nplzmax; nplz++ ) {
			double dtt = (defgmedian_td * ((double) nplz)) + myBoxes[mn].tz;
			dtt = dtt - pl->oz;
			dtt /= ww;

			newentry = true;
			for ( unsigned i = 0; i < dtilde.size(); i++) { if ( (SQR(dtt - dtilde[i])) < SQR(EPS_SEGMENTS_ELIMDOUBLE) ) newentry = false; }

			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); }
		}
	}

	//sort this list in ascending order
	std::sort( dtilde.begin(), dtilde.end(), SortDblAscending );

	//do not forget to put the totallength, as for in particular necessary when the raytracing process did not find any interception wiut a defgr boundary
	dtilde.push_back( plen );

	//MK::elimination of reduced dtilde values for interceptions at intersection lines of planes is now obsolete 
	//because there are cant be any doubles as the trace is parameterized and dtilde the only free parameter but the trace
	//is confied along the direction direcx,y,z in the box, in the case that the ray cuts exactly in a corner of a grain
	//it intersects both with a pair xy, yz or xz und thus is exactly one time accounted for!

	//get to know how many intercepts there are, make sure to take care of the value dtilde == path.totallength where the nbor is located

	//allocate memory for local information
	int nseg = dtilde.size();
	pl->nsegments = nseg;

	//##MK::in the case of the SPHEREMESH_MODEL the outbound path is not necessary...
	regionpairP defms = NULL;
	defms = new regionpair[nseg];
	if (defms == NULL) { cout << "ERROR::During memory allocation for defms path in rank,nuc,nbor " << myRank << ";" << mn << endl; return; }
	this->myMemGuard += nseg * sizeof(regionpair);

	//get the coordinates of the points halfway in each section to assign the ID of the corresponding grains 
	int ix, iy, iz, ixyz;

	//only OUTBOUND
	double dtraveled_outb; //MK::there is at least one element, namely the monocrystal!

	//if along the path there is no internal structure, i.e. nuc and nbor in the same homogeneous deformed grain!
	if ( dtilde.size() == 1 ) {
		dtraveled_outb = dtilde[0] * 0.5;
		ix = ((pl->ox + (dtraveled_outb * uu) - myBoxes[mn].tx) / defgmedian_rd); //rounds off
		iy = ((pl->oy + (dtraveled_outb * vv) - myBoxes[mn].ty) / defgmedian_nd);
		iz = ((pl->oz + (dtraveled_outb * ww) - myBoxes[mn].tz) / defgmedian_td);
		ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);

		//cout << "ixyz;out;" << ixyz << "\t"; //VERBOSE

		defms[0].dgid = this->myDefMS[mn][ixyz];
	}
	if ( dtilde.size() > 1 ) {
		//add the first
		dtraveled_outb = dtilde[0] * 0.5;
		ix = ((pl->ox + (dtraveled_outb * uu) - myBoxes[mn].tx) / defgmedian_rd); //rounds off
		iy = ((pl->oy + (dtraveled_outb * vv) - myBoxes[mn].ty) / defgmedian_nd);
		iz = ((pl->oz + (dtraveled_outb * ww) - myBoxes[mn].tz) / defgmedian_td);
		ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);

		//cout << "ixyz;out;" << ixyz << "\t";

		defms[0].dgid = this->myDefMS[mn][ixyz];

		//if dtilde.size > 1 then dtilde.size is at least 2, so access to element id 1 is valid and in that case following for loop executed once
		for ( long outbound = 1; outbound < dtilde.size(); outbound++ ) {
			dtraveled_outb = dtilde[outbound-1] + ((dtilde[outbound] - dtilde[outbound-1]) * 0.5); //list is in ascending order ELIMDOUBLE apart
			ix = ((pl->ox + (dtraveled_outb * uu) - myBoxes[mn].tx) / defgmedian_rd); 
			iy = ((pl->oy + (dtraveled_outb * vv) - myBoxes[mn].ty) / defgmedian_nd); 
			iz = ((pl->oz + (dtraveled_outb * ww) - myBoxes[mn].tz) / defgmedian_td); 
			ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);

			//cout << "ixyz;out;" << ixyz << "\t";

			defms[outbound].dgid = this->myDefMS[mn][ixyz];
		}
	}

	//now, construct the regionends for the forward iterate for OUTBOUND (simple the dtilde values)
	for ( unsigned long outb = 0; outb < dtilde.size(); outb++ ) {
		defms[outb].regionends = dtilde[outb];
	}

	pl->grainregions = defms;


	/*VERBOSE ...
	//##debug check content of the grains
	cout << "grainregions_out::" << endl;
	for ( unsigned k = 0; k < nseg; k++) {
		cout << "(" << pl->grainregions[k].dgid << "|" << pl->grainregions[k].regionends << ")";
	}
	cout << endl << "grainregions_in::" << endl;
	*/


	dtilde.clear();
}


/*
//##MK::NOT WORKING CORRECTLY!
void polyxx::sim_pathsegmentation( pathP pl, long mn, long nb )
{
	vector<double> dtilde;
	dtilde.reserve ( myBoxes[mn].ngrxyz ); //dtilde is auto-freed after function exit
	bool newentry = true;

	//find interceptions with x, y, z planes by calculating dtilde of the path parameterization (ox,oy,oz)**T + dtilde * (u,v,w)**T
	double plen = pl->totallength;

	//tx,ty,tz shift of the grain structure is always < 1.0
	//XPlanes first
	double uu = pl->direcx;
	if ( fabs(uu) > EPS_INPLANE_ACCURACY ) {
		long nplxmax = myBoxes[mn].ngrx;

		for ( long nplx = 1; nplx <= nplxmax; nplx++ ) { //only if nplx strictly smaller so no boundary at only one grain in x direction, MK::however will be catched by < plen anyway
			double dtt = (defgmedian_rd * ((double) nplx) ) + myBoxes[mn].tx; //+ because txyz < 0, oxyz > 0
			dtt = dtt - pl->ox;
			dtt /= uu; //now division is never by 0.0

			//only when dtt is positive, < pl->totallength, and the value to EPS_SEGMENTS_ELIMDOUBLE accuracy not already there is a valid interception, otherwise wrong direction (dtt < 0), plane cutting at ox ( dtt == 0.0 ), outside the box
			newentry = true;
			for ( unsigned int i = 0; i < dtilde.size(); i++) { if (( dtt - dtilde[i]) < EPS_SEGMENTS_ELIMDOUBLE ) newentry = false; }

			//negative dtt cannot be accepted as we consider only the path as a path segment starting at nucsite and pointing in the direction of the nbors in the local coordinate system
			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); } //MK::could also be made <= plen
		} //all planes
	} //else, dtilde is still left empty if path lays in the yz plane

	//YPlanes
	double vv = pl->direcy;
	if ( fabs(vv) > EPS_INPLANE_ACCURACY ) {
		long nplymax = myBoxes[mn].ngry;

		for ( long nply = 1; nply <= nplymax; nply++ ) {
			double dtt = (defgmedian_nd * ((double) nply) ) + myBoxes[mn].ty;
			dtt = dtt - pl->oy;
			dtt /= vv; //now division is never by 0.0

			newentry = true;
			for ( unsigned i = 0; i < dtilde.size(); i++) { if (( dtt - dtilde[i]) < EPS_SEGMENTS_ELIMDOUBLE ) newentry = false; }

			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); }
		}
	}

	//ZPlanes
	double ww = pl->direcz;
	if ( fabs(ww) > EPS_INPLANE_ACCURACY ) {
		long nplzmax = myBoxes[mn].ngrz;

		for ( long nplz = 1; nplz <= nplzmax; nplz++ ) {
			double dtt = (defgmedian_td * ((double) nplz)) + myBoxes[mn].tz;
			dtt = dtt - pl->oz;
			dtt /= ww; //now division is never by 0.0

			newentry = true;
			for ( unsigned i = 0; i < dtilde.size(); i++) { if (( dtt - dtilde[i]) < EPS_SEGMENTS_ELIMDOUBLE ) newentry = false; }

			if ( dtt > 0.0 && dtt < plen && newentry == true ) { dtilde.push_back ( dtt ); }
		}
	}

	//sort this list in ascending order
	std::sort( dtilde.begin(), dtilde.end(), SortDblAscending );

	//do not forget to put the totallength, as for in particular necessary when the raytracing process did not find any interception wiut a defgr boundary
	dtilde.push_back( plen );

	//MK::elimination of reduced dtilde values for interceptions at intersection lines of planes is now obsolete 
	//because there are cant be any doubles as the trace is parameterized and dtilde the only free parameter but the trace
	//is confied along the direction direcx,y,z in the box, in the case that the ray cuts exactly in a corner of a grain
	//it intersects both with a pair xy, yz or xz und thus is exactly one time accounted for!

	//get to know how many intercepts there are, make sure to take care of the value dtilde == path.totallength where the nbor is located
	//VERBOSE:
	//cout << myRank << "\t\t" << mn << "\t\t" << nb << "\t\t" << dtilde.size() << "\t\t" << plen << "=??=" << dtilde[dtilde.size() - 1] << "nbor at=" << myNeighbors[mn][nb].x << ";" << myNeighbors[mn][nb].y << ";" <<  myNeighbors[mn][nb].z << endl;
	//
	7VERBOSE
	//for ( unsigned i = 0; i < dtilde.size(); i++) {
	//	cout << dtilde[i] << ";";
	//}
	cout << endl;
	//VERBOSE END
	//allocate memory for local information
	int nseg = dtilde.size();
	pl->nsegments = nseg;

	//##MK::in the case of the SPHEREMESH_MODEL the outbound path is not necessary...
	long* defms = NULL;
	defms = new long[nseg];
	if (defms == NULL) { cout << "ERROR::During memory allocation for defms path in rank,nuc,nbor " << myRank << ";" << mn << ";" << nb << endl; return; }
	this->myMemGuard += nseg * sizeof(long);
	//mirror order
	long* _defms = NULL;
	_defms = new long[nseg];
	if ( _defms == NULL) { cout << "ERROR::During memory allocation for defms path in rank,nuc,nbor " << myRank << ";" << mn << ";" << nb << endl; return; }
	this->myMemGuard += nseg * sizeof(long);
	
	//get the coordinates of the points halfway in each section to assign the ID of the corresponding grains 
	int ix, iy, iz, ixyz;
	//first OUTBOUND
	double dtraveled_outb; //MK::there is at least one element, namely the monocrystal!

	//if along the path there is no internal structure, i.e. nuc and nbor in the same homogeneous deformed grain!
	if ( dtilde.size() == 1 ) {
		dtraveled_outb = dtilde[0] * 0.5;
		ix = ((pl->ox + (dtraveled_outb * uu) - myBoxes[mn].tx) / defgmedian_rd); //rounds off
		iy = ((pl->oy + (dtraveled_outb * vv) - myBoxes[mn].ty) / defgmedian_nd);
		iz = ((pl->oz + (dtraveled_outb * ww) - myBoxes[mn].tz) / defgmedian_td);
		ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);
		//VERBOSE:
		//cout << "ixyz;out;" << ixyz << "\t";
		//
		defms[0] = this->myDefMS[mn][ixyz];
	}
	if ( dtilde.size() > 1 ) {
		//add the first
		dtraveled_outb = dtilde[0] * 0.5;
		ix = ((pl->ox + (dtraveled_outb * uu) - myBoxes[mn].tx) / defgmedian_rd); //rounds off
		iy = ((pl->oy + (dtraveled_outb * vv) - myBoxes[mn].ty) / defgmedian_nd);
		iz = ((pl->oz + (dtraveled_outb * ww) - myBoxes[mn].tz) / defgmedian_td);
		ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);
		//VERBOSE:
		//cout << "ixyz;out;" << ixyz << "\t";
		//VERBOSE
		defms[0] = this->myDefMS[mn][ixyz];

		//if dtilde.size > 1 then dtilde.size is at least 2, so access to element id 1 is valid and in that case following for loop executed once
		for ( long outbound = 1; outbound < dtilde.size(); outbound++ ) {
			dtraveled_outb = dtilde[outbound-1] + ((dtilde[outbound] - dtilde[outbound-1]) * 0.5); //list is in ascending order ELIMDOUBLE apart
			ix = ((pl->ox + (dtraveled_outb * uu) - myBoxes[mn].tx) / defgmedian_rd); 
			iy = ((pl->oy + (dtraveled_outb * vv) - myBoxes[mn].ty) / defgmedian_nd); 
			iz = ((pl->oz + (dtraveled_outb * ww) - myBoxes[mn].tz) / defgmedian_td); 
			ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);
			//VERBOSE:
			//cout << "ixyz;out;" << ixyz << "\t";
			//VERBOSE
			defms[outbound] = this->myDefMS[mn][ixyz];
		}
	}
	pl->grainregions_outbound = defms;

	//second INBOUND
	//MK::solver always iterates forward!
	double dtraveled_inb;

	if ( dtilde.size() == 1 ) {
		dtraveled_inb = dtilde[0] * 0.5;
		ix = ( (pl->ox + (dtraveled_inb * uu) - myBoxes[mn].tx) / defgmedian_rd);
		iy = ( (pl->oy + (dtraveled_inb * vv) - myBoxes[mn].ty) / defgmedian_nd);
		iz = ( (pl->oz + (dtraveled_inb * ww) - myBoxes[mn].tz) / defgmedian_td); //############was ix! n00b
		ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);
		//VERBOSE
		//cout << "ixyz;out;" << ixyz << "\t";
		//
		_defms[0] = this->myDefMS[mn][ixyz];

//###if ( _defms[0] > this->ndefgpool ) {cout << "ERRTilde1::" << dtilde.size() << " dtilde[0] = " << dtilde[0] << " nseg " << nseg << "\t\t" << _defms[0] << "|ixyz|" << ix << "-" << iy << "-" << iz << "-" << ixyz << "|ngrxyz|" << myBoxes[mn].ngrx << "-" << myBoxes[mn].ngry << "-" << myBoxes[mn].ngrz << "-" << myBoxes[mn].ngrxyz << endl;}
	}
	if ( dtilde.size() > 1 ) {
		int wheretoplace;

		for ( long inbound = 1; inbound < dtilde.size(); inbound++ ) { //runs to < size = 4
			wheretoplace = inbound - 1;
			dtraveled_inb = dtilde[(dtilde.size() - inbound - 1)] + ((dtilde[(dtilde.size() - inbound)] - dtilde[(dtilde.size() - inbound - 1)]) * 0.5); //###inbound < size so size - inbound = 1 so 1-1 = 0 no invalid access!
			ix = ( (pl->ox + (dtraveled_inb * uu) - myBoxes[mn].tx) / defgmedian_rd);
			iy = ( (pl->oy + (dtraveled_inb * vv) - myBoxes[mn].ty) / defgmedian_nd);
			iz = ( (pl->oz + (dtraveled_inb * ww) - myBoxes[mn].tz) / defgmedian_td);
			ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);
			//VERBOSE:
			//cout << "ixyz;out;" << ixyz << "\t";
			//
			_defms[wheretoplace] = this->myDefMS[mn][ixyz];

//###if ( _defms[wheretoplace] > this->ndefgpool ) {cout << "ERRTilde>1ForLongInbound::" << _defms[wheretoplace] << endl;}
		}
		//append last grain
		wheretoplace = dtilde.size() - 1;
		dtraveled_inb = dtilde[0] * 0.5; //MK::##MAYBE WRONG?
		ix = ( (pl->ox + (dtraveled_inb * uu) - myBoxes[mn].tx) / defgmedian_rd);
		iy = ( (pl->oy + (dtraveled_inb * vv) - myBoxes[mn].ty) / defgmedian_nd);
		iz = ( (pl->oz + (dtraveled_inb * ww) - myBoxes[mn].tz) / defgmedian_td);
		ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz); 
		//VERBOSE:
		//cout << "ixyz;out;" << ixyz << "\t";
		//
		_defms[wheretoplace] = this->myDefMS[mn][ixyz];

//###if ( _defms[wheretoplace] > this->ndefgpool ) {cout << "ERRTilde>1AppLast::" << _defms[wheretoplace] << endl;}
	}
	pl->grainregions_inbound = _defms;

	//now, construct the regionends for the forward iterate for OUTBOUND (simple the dtilde values) and INBOUND (totallength - dtilde[i]) growth
	double* pathsegend = NULL;
	pathsegend = new double[nseg];
	if (pathsegend == NULL) { cout << "ERROR::During memory allocation for regionend path in rank,nuc,nbor " << myRank << ";" << mn << ";" << nb << endl; return; }
	this->myMemGuard += nseg * sizeof(double);
	//OUTBOUND, trivial
	for ( long outb = 0; outb < dtilde.size(); outb++ ) {
		pathsegend[outb] = dtilde[outb];
	}
	pl->regionends_outbound = pathsegend;

	//INBOUND needs care!
	double* _pathsegend = NULL;
	_pathsegend = new double[nseg];
	if (_pathsegend == NULL) { cout << "ERROR::During memory allocation for regionend path in rank,nuc,nbor " << myRank << ";" << mn << ";" << nb << endl; return; }
	this->myMemGuard += nseg * sizeof(double);
	if ( dtilde.size() == 1 ) {
		_pathsegend[0] = pl->totallength;
	}
	if ( dtilde.size() > 1 ) {
		for ( long inb = 1; inb < dtilde.size(); inb++ ) {
			_pathsegend[inb-1] = pl->totallength - dtilde[(dtilde.size() - inb - 1)]; //should be equivalent to dtilde[dtilde.size()-1] - dtilde[dtilde.size()-1]
		}
		//append last to avoid access violation
		_pathsegend[(dtilde.size() - 1)] = pl->totallength;
 	}
	pl->regionends_inbound = _pathsegend;

	//##debug check content of the grains
	//cout << "grainregions_out::" << endl;
	//for ( unsigned k = 0; k < nseg; k++) cout << pl->grainregions_outbound[k] << ";";
	//cout << endl << "grainregions_in::" << endl;
	//for ( unsigned k = 0; k < nseg; k++) cout << pl->grainregions_inbound[k] << ";";
	//cout << endl << "regionend_out::" << endl;
	//for ( unsigned k = 0; k < nseg; k++) cout << pl->regionends_outbound[k] << ";";
	//cout << endl << "regionend_in::" << endl;
	//for ( unsigned k = 0; k < nseg; k++) cout << pl->regionends_inbound[k] << ";";
	//cout << endl;

	dtilde.clear();
}
*/


void print_polyhedron(vector<face> faces,vector<vertex> vertices, long n_faces, string filename){
	voro::voronoicell polyh_cell;
	char *fname = new char[filename.length() + 1];
	strcpy(fname, filename.c_str());
	polyh_cell.init(-5,5,-5,5,-5,5);
	double m_x,m_y,m_z;
	for (long i=0; i<n_faces; i++){
		m_x=(vertices[faces[i].A].x+vertices[faces[i].B].x+vertices[faces[i].C].x)*2*ONETHIRD;
		m_y=(vertices[faces[i].A].y+vertices[faces[i].B].y+vertices[faces[i].C].y)*2*ONETHIRD;
		m_z=(vertices[faces[i].A].z+vertices[faces[i].B].z+vertices[faces[i].C].z)*2*ONETHIRD;
		polyh_cell.plane(m_x,m_y,m_z);
	}
	polyh_cell.draw_pov(0.,0.,0.,fname);
	delete [] fname;
}


void polyxx::sim_ICO_HEDRONMESHING()
{
	//PH::generate icosphere Mesh (vertex + facial information)
	vector<face> icoh_faces;
	long curr_faces, curr_vertices;
	long n_vertices = 12 + 20 * pow(3, meshRecursion);
	ico_vertices.reserve( n_vertices ); //MK::one page automatically

	//PH::list longer than necessary but avoid such case decision which vertices are already existing

	double t = (1 + sqrt(5.0)) * 0.5; //golden ratio
	vertex v;
	face f;

	//Generate all 12 icosahedron vertices
	//vertices 1-4
	v.x = -1;		v.y = t;		v.z = 0;		ico_vertices.push_back(v);
	v.x = 1;		v.y = t;		v.z = 0;		ico_vertices.push_back(v);
	v.x = -1;		v.y = -t;		v.z = 0;		ico_vertices.push_back(v);
	v.x = 1;		v.y = -t;		v.z = 0;		ico_vertices.push_back(v);
	//vertices 5-8
	v.x = 0;		v.y = -1;		v.z = t;		ico_vertices.push_back(v);
	v.x = 0;		v.y = 1;		v.z = t;		ico_vertices.push_back(v);
	v.x = 0;		v.y = -1;		v.z = -t;		ico_vertices.push_back(v);
	v.x = 0;		v.y = 1;		v.z = -t;		ico_vertices.push_back(v);
	//vertices 9-12
	v.x = t;		v.y = 0;		v.z = -1;		ico_vertices.push_back(v);
	v.x = t;		v.y = 0;		v.z = 1;		ico_vertices.push_back(v);
	v.x = -t;		v.y = 0;		v.z = -1;		ico_vertices.push_back(v);
	v.x = -t;		v.y = 0;		v.z = 1;		ico_vertices.push_back(v);
	
	curr_vertices = 12;

	//Generate all 20 icosahedron faces, each having three vertices A,B,C
    //ico_faces.reserve(20); //MK::each vector automatically instantiates a page

	//##MK::where is this from?
	//PH::assign the edges to the vertices
	//Faces 1-5
	f.A = 0;		f.B = 11;		f.C = 5;		ico_faces.push_back(f);
	f.A = 0;		f.B = 5;		f.C = 1;		ico_faces.push_back(f);
	f.A = 0;		f.B = 1;		f.C = 7;		ico_faces.push_back(f);
	f.A = 0;		f.B = 7;		f.C = 10;		ico_faces.push_back(f);
	f.A = 0;		f.B = 10;		f.C = 11;		ico_faces.push_back(f);
	//Faces 6-10
	f.A = 1;		f.B = 5;		f.C = 9;		ico_faces.push_back(f);
	f.A = 5;		f.B = 11;		f.C = 4;		ico_faces.push_back(f);
	f.A = 11;		f.B = 10;		f.C = 2;		ico_faces.push_back(f);
	f.A = 10;		f.B = 7;		f.C = 6;		ico_faces.push_back(f);
	f.A = 7;		f.B = 1;		f.C = 8;		ico_faces.push_back(f);
	//Faces 11-15
    f.A = 3;		f.B = 9;		f.C = 4;		ico_faces.push_back(f);
	f.A = 3;		f.B = 4;		f.C = 2;		ico_faces.push_back(f);
	f.A = 3;		f.B = 2;		f.C = 6;		ico_faces.push_back(f);
	f.A = 3;		f.B = 6;		f.C = 8;		ico_faces.push_back(f);
	f.A = 3;		f.B = 8;		f.C = 9;		ico_faces.push_back(f);
	//Faces 16-20
	f.A = 4;		f.B = 9;		f.C = 5;		ico_faces.push_back(f);
	f.A = 2;		f.B = 4;		f.C = 11;		ico_faces.push_back(f);
	f.A = 6;		f.B = 2;		f.C = 10;		ico_faces.push_back(f);
	f.A = 8;		f.B = 6;		f.C = 7;		ico_faces.push_back(f);
	f.A = 9;		f.B = 8;		f.C = 1;		ico_faces.push_back(f);

	curr_faces = 20;

#ifdef TESTWITHVOROXX //
	print_polyhedron(ico_faces, ico_vertices, curr_faces, "icosahedron_data.pov" );
#endif
	//##MK:::certainly minor cache benefit also in restructuring the instantiation order of the vertices and faces

	long new_faces;
	//for each recursion add 3 new vertices and 4 new faces
	for( int i = 0; i < meshRecursion; i++ ) { 

		icoh_faces.reserve( 20 * pow(4, i+1) );
		new_faces = 0;

		for( long i_f = 0; i_f < curr_faces; i_f++){
			//create new vertices
			//vertice 1
			v.x = 0.5*(ico_vertices[ico_faces[i_f].A].x + ico_vertices[ico_faces[i_f].B].x);
			v.y = 0.5*(ico_vertices[ico_faces[i_f].A].y + ico_vertices[ico_faces[i_f].B].y);
			v.z = 0.5*(ico_vertices[ico_faces[i_f].A].z + ico_vertices[ico_faces[i_f].B].z);
			ico_vertices.push_back(v);
			//vertice 2
			v.x = 0.5*(ico_vertices[ico_faces[i_f].A].x + ico_vertices[ico_faces[i_f].C].x);
			v.y = 0.5*(ico_vertices[ico_faces[i_f].A].y + ico_vertices[ico_faces[i_f].C].y);
			v.z = 0.5*(ico_vertices[ico_faces[i_f].A].z + ico_vertices[ico_faces[i_f].C].z);
			ico_vertices.push_back(v);
			//vertice 3
			v.x = 0.5*(ico_vertices[ico_faces[i_f].B].x + ico_vertices[ico_faces[i_f].C].x);
			v.y = 0.5*(ico_vertices[ico_faces[i_f].B].y + ico_vertices[ico_faces[i_f].C].y);
			v.z = 0.5*(ico_vertices[ico_faces[i_f].B].z + ico_vertices[ico_faces[i_f].C].z);
			ico_vertices.push_back(v);

			//generate new faces and save them in help vector icoh
			//face I: A-2-1 counter clockwise roundtrip after clockwise definition of the triangular face
			f.A = ico_faces[i_f].A;
			f.B = curr_vertices + 1;
			f.C = curr_vertices;
			icoh_faces.push_back(f);
			//face II: 1-3-B
			f.A = curr_vertices;
			f.B = curr_vertices+2;
			f.C = ico_faces[i_f].B;
			icoh_faces.push_back(f);
			//face III: 2-C-3
			f.A = curr_vertices+1;
			f.B = ico_faces[i_f].C;
			f.C = curr_vertices+2;
			icoh_faces.push_back(f);
			//face IV:  1-2-3
			f.A = curr_vertices;
			f.B = curr_vertices+1;
			f.C = curr_vertices+2;
			icoh_faces.push_back(f);
	
			curr_vertices += 3;
			new_faces += 4;
		}
		curr_faces = new_faces;
		ico_faces = icoh_faces;
		icoh_faces.clear();
	}
//cout << "curr faces = " << curr_faces << " curr_vertices = " <<  curr_vertices<< endl;

	numFaces = curr_faces;
	//PH::normalizing each vertex maps each vertex touch the unit sphere
	double l;
	for( int i = 0; i < curr_vertices; i++ ) {
		l = pow( ( SQR(ico_vertices[i].x) + SQR(ico_vertices[i].y) + SQR(ico_vertices[i].z) ), 0.5 ); //l never zero because isocahedron located at the center
		//PH::thus making the unit sphere better and better approximated by scaling the icosphere from the inside to touch the unit sphere...
		//http//blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
		ico_vertices[i].x /= l;
		ico_vertices[i].y /= l;
		ico_vertices[i].z /= l;
	}

	//up to here the icosahedron is in unit coordinates but not in absolute coordinates...

	//calculate facial area
	vertex a,b;	
	for ( long i_f = 0; i_f < numFaces; i_f++ ) {
		//edge a: A-B
		ico_faces[i_f].area = 0.0;
		
		//##PH::1E6 - possible for make a quick and yes a dirty... scaling to avoid num flaws...
		a.x = 1E6 * (ico_vertices[ico_faces[i_f].B].x - ico_vertices[ico_faces[i_f].A].x );
		a.y = 1E6 * (ico_vertices[ico_faces[i_f].B].y - ico_vertices[ico_faces[i_f].A].y );
		a.z = 1E6 * (ico_vertices[ico_faces[i_f].B].z - ico_vertices[ico_faces[i_f].A].z );
		//edge b: A-C
		b.x = 1E6 * (ico_vertices[ico_faces[i_f].C].x - ico_vertices[ico_faces[i_f].A].x );
		b.y = 1E6 * (ico_vertices[ico_faces[i_f].C].y - ico_vertices[ico_faces[i_f].A].y );
		b.z = 1E6 * (ico_vertices[ico_faces[i_f].C].z - ico_vertices[ico_faces[i_f].A].z );
	
		//calculate area from cross product
		//A_triangle=0.5*|CP(a,b)|
		ico_faces[i_f].area += SQR(a.y*b.z - a.z*b.y);
		ico_faces[i_f].area += SQR(a.x*b.z - a.z*b.x); //signs in CP (-/+) does not matter only because magnitude is aimed for
		ico_faces[i_f].area += SQR(a.x*b.y - a.y*b.x);
                ico_faces[i_f].area = pow( ico_faces[i_f].area, 0.5 ); //length
		ico_faces[i_f].area = (0.5 * SQR(1E-6)) * ico_faces[i_f].area;
	}
	//if desired to compare for approximation of unit sphere 4/3_PI_CUBE(1.0)
	polyvol = init_CALCICOISOVOLUME();

#ifdef POVRAYIO
	print_polyhedron( ico_faces, ico_vertices, curr_faces, "polyhedron_data.pov" );
#endif
}


void polyxx::sim_ICO_NUCPATHCONSTRUCTION( void )
{
	//MK::generate facenormals and corresponding path
	double l, ltest;
	vertex_l fn;
	ico_facenormals.reserve(numFaces);

	//calculate normalized facevectors and length
	//#####but the normal vector fn has a direction, why not taking the cross product?
	for ( long i_f = 0; i_f < numFaces; i_f++ ) {
		//center of the triangle here the icosphere is still center in an own coordinate system therefore the direction of the unit vector is known
		fn.x = ( ico_vertices[ico_faces[i_f].A].x + ico_vertices[ico_faces[i_f].B].x + ico_vertices[ico_faces[i_f].C].x ) * ONETHIRD;
		fn.y = ( ico_vertices[ico_faces[i_f].A].y + ico_vertices[ico_faces[i_f].B].y + ico_vertices[ico_faces[i_f].C].y ) * ONETHIRD;
		fn.z = ( ico_vertices[ico_faces[i_f].A].z + ico_vertices[ico_faces[i_f].B].z + ico_vertices[ico_faces[i_f].C].z ) * ONETHIRD;
		//normalize as it is a direction vector
		l = pow( (SQR(fn.x) + SQR(fn.y) + SQR(fn.z)) , 0.5 ); 
		
		//meshing to unitsphere renders that no vertex is located in the center and therefore it is safe to assume l > DOUBLE_ACCURACY ie 0.0
		fn.x /= l;
		fn.y /= l;
		fn.z /= l;

		//wall-collision test
		//only because the the unitsphere is in the center of the domain, fabs includes +- cases
		l = fabs( 0.5 * boxx / fn.x ); //half box because nucleus is located at 0.5,0.5,0.5 relative coordinates
		ltest = fabs( 0.5 * boxy / fn.y );
		if ( ltest < l ) l = ltest;
		ltest = fabs( 0.5 * boxz / fn.z );
		if ( ltest < l ) l = ltest;

		//fn.l is the endpoint of the path along which the boundary migrates from the nucleation site ox,y,z to ox,y,z + fn.l * (u,v,w)
		fn.l = l - EPSILONBALL; //so fn.l is the shortest length, minus epsilonball to keep personal distance to the wall

		ico_facenormals.push_back(fn);
	}


	//construct and analyze the deformation microstructure along the path
	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {

		//grain boundary migration is inspected as a Turnbull rate migration travel problem along the tracetory of the mypathes
		vector<path> mypathes;

		double ox = myNuclei[mynuc].x;
		double oy = myNuclei[mynuc].y;
		double oz = myNuclei[mynuc].z;
		//loop over all faces
		for( long i_f = 0; i_f < numFaces; i_f++ ) {
			path pline;

			pline.ox = ox;
			pline.oy = oy;
			pline.oz = oz;
			pline.direcx = ico_facenormals[i_f].x;
			pline.direcy = ico_facenormals[i_f].y;
			pline.direcz = ico_facenormals[i_f].z;

			pline.distcontact = NOT_DETERMINED_YET;
			pline.timecontact = NOT_DETERMINED_YET;
			pline.nsegments = NOT_DETERMINED_YET;

			pline.totallength = ico_facenormals[i_f].l;

			pline.grainregions = NULL;

			this->sim_pathsegmentation_nucleus( &pline, mynuc );

			pline.analyzedYet = NO;

			mypathes.push_back(pline);

//			//##MK::DEBUG
//			for ( unsigned long p = 0; p < mypathes[mypathes.size()-1].nsegments; p++)
//			cout << "MYNUC---" << mynuc << "---i_f---" << i_f << "---" << mypathes[mypathes.size()-1].grainregions[p].dgid << "\t\t" << mypathes[mypathes.size()-1].grainregions[p].regionends << "---hit---" << (ox + (pline.direcx * mypathes[mypathes.size()-1].grainregions[p].regionends)) << "\t\t" << (oy + (pline.direcy * mypathes[mypathes.size()-1].grainregions[p].regionends)) << "\t\t" << (oz + (pline.direcz * mypathes[mypathes.size()-1].grainregions[p].regionends)) << endl;
//			cout << "MYNUC||" << ox << ";" << oy << ";" << oz << "---" << pline.direcx << ";" << pline.direcy << ";" << pline.direcz << endl;
//			##MK::DEBUG
		}

//		cout << "MYNUCBOX---" << myBoxes[mynuc].xrd << ";" <<myBoxes[mynuc].ynd << ";" << myBoxes[mynuc].ztd << "---" << myBoxes[mynuc].tx << ";" << myBoxes[mynuc].ty << ";" << myBoxes[mynuc].tz << "---" << myBoxes[mynuc].ngrx << ";" << myBoxes[mynuc].ngry << ";" << myBoxes[mynuc].ngrz << endl;

		myIcoPaths.push_back (mypathes);
	}
}


void polyxx::sim_ICO_NBORPATHCONSTRUCTION( void )
{
	//generate paths from nucleus to neighbors
	//for all my myIDs.size() generate the lines pointing outbound from the nucleus to its neighbors
	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
		vector<path> localpath;

		double onx, ony, onz; //the nuc
		double nx, ny, nz; //one its its nbors
		onx = myNuclei[mynuc].x;
		ony = myNuclei[mynuc].y;
		onz = myNuclei[mynuc].z;
		double l, ltest;
		double dist = 0.0;

		for (long nbor = 0; nbor < this->myNeighbors[mynuc].size(); nbor++) {
			struct path pline;
			//get path geometry
			nx = myNeighbors[mynuc][nbor].x; //##MK::this is already in MICRON global coordinates of the domain box
			ny = myNeighbors[mynuc][nbor].y;
			nz = myNeighbors[mynuc][nbor].z;

			double u, v, w;
			//MK::for direction on - n because neighbors grow towards the nucleus!!
			u = onx - nx;
			v = ony - ny;
			w = onz - nz;
			double len = SQR(u) + SQR(v) + SQR(w);
			len = pow(len, 0.5);
			//direcx,y,z is a direction vector pointing from ox,y,z to the nborx,y,z site
			if ( len < EPSILONBALL ) { cout << "ERROR during path construction node;mynuc;nbor " << myRank << ";" << mynuc << ";" << nbor << endl; return; }
			u /= len;
			v /= len;
			w /= len;

			//MK::wall-collision, constraints, hitting a wall at point px,py,pz = nxyz + (l*direcxyz), only accept l > 0.0 and in all these cases the minimum
			l = INFINITE;
			ltest = INFINITE;

			//x
			if ( fabs(u) > DOUBLE_ACCURACY ) {
				ltest = (0.0 - nx) / u;
				if ( ltest > 0.0 && ltest <= l ) l = ltest;

				ltest = (boxx - nx) / u;
				if ( ltest > 0.0 && ltest <= l ) l = ltest;
			}

			if ( fabs(v) > DOUBLE_ACCURACY ) {
				ltest = (0.0 - ny) / v;
				if ( ltest > 0.0 && ltest <= l ) l = ltest;

				ltest = (boxy - ny) / v;
				if ( ltest > 0.0 && ltest <= l ) l = ltest;
			}

			if ( fabs(w) > DOUBLE_ACCURACY ) {
				ltest = (0.0 - nz) / w;
				if ( ltest > 0.0 && ltest <= l ) l = ltest;

				ltest = (boxz - nz) / w;
				if ( ltest > 0.0 && ltest <= l ) l = ltest;
			}

			if ( l >= (INFINITE - DOUBLE_ACCURACY) ) {cout << "ERR::PATHDISTANCE IN NBOR PATHSEGMENTATION STILL INFINITE!" << endl; return; }

			//MKPH::good guys better come never too close to bad girls ... wall
			l -= EPSILONBALL;
			
			pline.totallength = l;
			pline.ox = nx; //##MK::wrong was that an additional  - u*l; was accounted for! //because l was necessary fabs before
			pline.oy = ny; // - v*l;
			pline.oz = nz; // - w*l;

			pline.direcx = u;
			pline.direcy = v;
			pline.direcz = w;

			pline.grainregions = NULL;


			pline.distcontact = NOT_DETERMINED_YET;
			pline.timecontact = NOT_DETERMINED_YET;
			pline.nsegments = NOT_DETERMINED_YET;
			pline.analyzedYet = NO;

			//analyze the deformed microstructure along each path here
			this->sim_pathsegmentation_nbors( &pline, mynuc );


			localpath.push_back ( pline ); //copy-constructor



//			//##MK::DEBUG
//			for ( unsigned long p = 0; p < localpath[localpath.size()-1].nsegments; p++) 
//				cout << "NBOR---" << nbor << "---" << localpath[localpath.size()-1].grainregions[p].dgid << "---" << localpath[localpath.size()-1].grainregions[p].regionends << "---hit---" << (nx + (pline.direcx * localpath[localpath.size()-1].grainregions[p].regionends)) << "\t\t" << (ny + (pline.direcy * localpath[localpath.size()-1].grainregions[p].regionends)) << "\t\t" << (nz + (pline.direcz * localpath[localpath.size()-1].grainregions[p].regionends)) << endl;
//			cout << "NBOR---" << nbor << "||" << nx << ";" << ny << ";" << nz << "---" << pline.direcx << ";" << pline.direcy << ";" << pline.direcz << endl;
//			//##MK::DEBUG

		} //for all neighbors

//		cout << "MYNUCBOX---" << myBoxes[mynuc].xrd << ";" <<myBoxes[mynuc].ynd << ";" << myBoxes[mynuc].ztd << "---" << myBoxes[mynuc].tx << ";" << myBoxes[mynuc].ty << ";" << myBoxes[mynuc].tz << "---" << myBoxes[mynuc].ngrx << ";" << myBoxes[mynuc].ngry << ";" << myBoxes[mynuc].ngrz << endl;

		this->myNBorPaths.push_back ( localpath );
	} //for all my nuclei
	//cache local disorientations along the path before analyzing the path
}



/*
void polyxx::sim_VORO_PATHCONSTRUCTION( void )
{
	//for all my myIDs.size() generate the lines pointing outbound from the nucleus to its neighbors
	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
		vector<path> localpath;

		double ox, oy, oz;
		double nx, ny, nz;
		ox = myNuclei[mynuc].x;
		oy = myNuclei[mynuc].y;
		oz = myNuclei[mynuc].z;

		double dist = 0.0;
		for (long nbor = 0; nbor < this->nNeighborsPerCA; nbor++) {
			struct path pline;

			pline.ox = ox;
			pline.oy = oy;
			pline.oz = oz;
			//get path geometry
			nx = myNeighbors[mynuc][nbor].x;
			ny = myNeighbors[mynuc][nbor].y;
			nz = myNeighbors[mynuc][nbor].z;

			//MK::for direction n - o is essential that the code gets the regions correctly settled out!!
			double u, v, w;
			u = nx - ox;
			v = ny - oy;
			w = nz - oz;
			double len = SQR(u) + SQR(v) + SQR(w);
			len = pow(len, 0.5);

			pline.totallength = len;
			//direcx,y,z is a direction vector pointing from ox,y,z to the nborx,y,z site
			if ( len < EPSILONBALL ) { cout << "ERROR during path construction node;mynuc;nbor " << myRank << ";" << mynuc << ";" << nbor << endl; return; }
			u /= len;
			v /= len;
			w /= len;

			pline.analyzedYet = NO;
			pline.direcx = u;
			pline.direcy = v;
			pline.direcz = w;

			//analyze the deformed microstructure along each path here
			this->sim_pathsegmentation( &pline, mynuc, nbor ); //VORO PATH CONSTRUCT

			pline.distcontact = NOT_DETERMINED_YET;
			pline.timecontact = NOT_DETERMINED_YET;
			
			localpath.push_back ( pline ); //copy-constructor

#ifdef DETAILED_PROMPTS
	cout << myRank << ";mynuc;nbor;myMemGuard" << mynuc << ";" << nbor << ";" << myMemGuard << endl;
#endif
		} //for all neighbors

		//push back the local path to the pathlist
		this->myNBorPaths.push_back ( localpath );

	} //for all my nuclei
	//cache local disorientations along the path before analyzing the path

	cout << myRank << " all growth paths constructed." << endl;
}
*/


void polyxx::init_CalculateDefTexVolumeGrainResolved( long mynuc,double *myNucTextureVolumeContribution)
{
	long ixyz;
	vector<double> xmin,xmax,ymin,ymax,zmin,zmax;

	for ( long i=0; i < oripool.size(); i++) {
			myNucTextureVolumeContribution[i] = 0.0;
	}

	xmin.reserve(myBoxes[mynuc].ngrx);
	xmax.reserve(myBoxes[mynuc].ngrx);
	
	ymin.reserve(myBoxes[mynuc].ngry);
	ymax.reserve(myBoxes[mynuc].ngry);
	
	zmin.reserve(myBoxes[mynuc].ngrz);
	zmax.reserve(myBoxes[mynuc].ngrz);

	//x
	//Calculate collision of deformed grain with box boundaries
	for( long ix = 0; ix < myBoxes[mynuc].ngrx; ix++ ) {
		xmin.push_back( ix * defgmedian_rd + myBoxes[mynuc].tx);
		xmax.push_back( xmin[ix] + defgmedian_rd);
		if( xmin[ix] < 0.0) xmin[ix] = 0.0;
		if( xmax[ix] > boxx ) xmax[ix] = boxx;
		//if(mynuc==0){cout<<scientific<<"tx "<<myBoxes[mynuc].tx<<" ix "<<ix<<" xmin "<<xmin[ix]<<" xmax "<<xmax[ix]<<endl;};
	}
	//y
	for( long iy=0; iy < myBoxes[mynuc].ngry; iy++ ) {
		ymin.push_back( iy * defgmedian_nd + myBoxes[mynuc].ty );
		ymax.push_back( ymin[iy] + defgmedian_nd );
		if( ymin[iy] < 0.0 ) ymin[iy] = 0.0;
		if( ymax[iy] > boxy ) ymax[iy] = boxy;
	}
	//z
	for( long iz=0; iz < myBoxes[mynuc].ngrz; iz++ ) {
		zmin.push_back( iz * defgmedian_td + myBoxes[mynuc].tz );
		zmax.push_back( zmin[iz] + defgmedian_td );
		if( zmin[iz] < 0.0 ) zmin[iz] = 0.0;
		if( zmax[iz] > boxz ) zmax[iz] = boxz;
	}
	//Now calculate deformed-grain-volume which lays inside the computation box
	//And add this Volume to Texture for corresponding ori id
	//grains in x layer stacked in y stacked in z
	
	for( long iz=0; iz < myBoxes[mynuc].ngrz; iz++) {
		long zoff = myBoxes[mynuc].ngrxy * iz;

		for( long iy=0; iy < myBoxes[mynuc].ngry; iy++ ) {
			long yzoff = (myBoxes[mynuc].ngrx * iy) + zoff;

			for(long ix=0; ix < myBoxes[mynuc].ngrx; ix++ ) {
				ixyz = ix + yzoff;

				myNucTextureVolumeContribution[defgpool[myDefMS[mynuc][ixyz]].ori] += ( xmax[ix] - xmin[ix] ) * ( ymax[iy] - ymin[iy] ) *( zmax[iz] - zmin[iz] ); //##MK::still room for improvement caches!
			}
		}
	}
		
	xmin.clear(); xmax.clear();
	ymin.clear(); ymax.clear();
	zmin.clear(); zmax.clear();
}


void polyxx::sim_ICO_NBORGROWTH_1tstep ( long &mynuc, double &t_total, double &deltat ) // vector<long> &nbor_segids, vector<double> &nbor_rads) 
{
	double v, m , P , s;
	double dt, dist_seg;
	long defgid, nbor;
	long t = 0;
	for( nbor = 0; nbor < myNeighbors[mynuc].size(); nbor++ ) {
		dt = t_total - myNeighbors[mynuc][nbor].tincub + deltat;

		if ( dt > deltat )
			dt = deltat;
		//otherwise while loop catches dt < 0.0
		//calculate neighbor-specific values
		long rxgid = myNeighbors[mynuc][nbor].rxgid;
		long rxgoid = rxgpool[rxgid].ori;

		//migrate step by step to potential boundaries between different deformed grains

/*//DEBUG
defgid = myNBorPaths[mynuc][nbor].grainregions_inbound[nbor_segids[nbor]];
if ( defgid > ndefgpool ) {
	cout << defgid << "--" << nbor_segids[nbor] << "--nsegments" << myNBorPaths[mynuc][nbor].nsegments << endl;
	for (long jj=0; jj<myNBorPaths[mynuc][nbor].nsegments; jj++) { 	cout << "\t\t--" << myNBorPaths[mynuc][nbor].grainregions_inbound[jj] << "--" << endl; }
	cout << endl << endl << myRank << ";ndegpool;" << ndefgpool << ";defgid;" << defgid << ";" << rxgoid << ";" << this->GSHact << ";" << this->myBoxes[mynuc].tx << "|" << this->myBoxes[mynuc].ty << "|" << this->myBoxes[mynuc].tz << "||" << this->myBoxes[mynuc].ngrxyz << endl; 
}*/


		while ( dt > SMALLTIME && nbor_segids[nbor] < myNBorPaths[mynuc][nbor].nsegments) {
			defgid = myNBorPaths[mynuc][nbor].grainregions[nbor_segids[nbor]].dgid; //grainregions_inbound[nbor_segids[nbor]];

			//how far am I away from the next but potentially different property region? MIND THAT regionends_inbound MARKS THE ENDPOINT OF THE PATH
			dist_seg = myNBorPaths[mynuc][nbor].grainregions[nbor_segids[nbor]].regionends - nbor_rads[nbor]; //regionends_inbound[nbor_segids[nbor]] - nbor_rads[nbor];


			P = get_lowtrig_mobilityweight( rxgoid, defgpool[defgid].ori );
                        if ( mobilitymodel_option == SEBALD_GOTTSTEIN ) {
                            m = P * mGS + (1 - P) * mHAGB;
                            if (P < 0.0) 
				m = mLAGB;
                        }
                        if ( mobilitymodel_option == ROLLETT_HUMPHREYS ) {
                            m = P * mRHHAGB;
                        }

			//Turnbull rate model
			v = m * (physConstants.get_halfG_b2() * get_rho( defgid ) - 0.0); //##MK::currently no ZenerDrag get_rho(defgid, t)
			//if(mynuc==0&&nbor==0){cout<<" vN= "<<v<< " Gbhalfsq=" <<Gbhalfsq<<" m= "<<m<<"P= "<<P<<endl;t=1;}
			//attempt to jump for s through this homogeneous deformation microstructure region
			s = v * dt;


			//we observe it is far more likely that we continue migrating in a region than crossing so lift load from the branch predictor
			if ( s < dist_seg ) { //##MK:: no <= slight FLP flaw maybe...
				nbor_rads[nbor] += s;
				dt = 0.0;
				continue;
			}

			//else s > dist_seg
			s = dist_seg;
			nbor_segids[nbor]++;
			dt -= s/v; //division by zero not possible because then s certainly almost zero and such already continued
			nbor_rads[nbor] += s;
		} 
		//timestep growth finished->next neighbor 
	}
}


void polyxx::sim_ICO_NUCGROWTH_1tstep(long &mynuc, double &t_total, double &deltat, long &numCeased, double *myCurrentNucTexContrib, double* myTotalDeformedTexture)
{
	double v, m , P, s;
	double dt_nuc, dt, dist_seg;
	double dx, dy, dz, dV;
	long defgid;
	
	long rxgid = myNuclei[mynuc].rxgid;
	long rxgoid = rxgpool[rxgid].ori;
	
	//how much time to integrate?
	dt_nuc = t_total - myNuclei[mynuc].tincub + deltat;

	if ( dt_nuc > deltat )	
		dt_nuc = deltat;
	//otherwise while loop catches dt < 0.0
	
	for ( long n_f = 0; n_f < numFaces; n_f++ ) {
		dt = dt_nuc;

		if( nuc_ceased[n_f] == UNCEASED ) { //calculate growth as long as collision has not occurred
//cout << "Nucgrowth--->mynuc/nf = " << mynuc << "--" << n_f << " unceased." << endl;

			while ( dt > SMALLTIME && nuc_segids[n_f] < myIcoPaths[mynuc][n_f].nsegments ) {
				//calculate nucleus-specific values
				defgid = myIcoPaths[mynuc][n_f].grainregions[nuc_segids[n_f]].dgid; //grainregions_outbound[nuc_segids[n_f]];
				dist_seg = myIcoPaths[mynuc][n_f].grainregions[nuc_segids[n_f]].regionends - nuc_rads[n_f]; //regionends_outbound[nuc_segids[n_f]] - nuc_rads[n_f];

				P = get_lowtrig_mobilityweight( rxgoid, defgpool[defgid].ori );
                                if ( mobilitymodel_option == SEBALD_GOTTSTEIN ) { 
                                    m = P * mGS + (1 - P) * mHAGB;
                                    if (P < 0.0) 
                                        m = mLAGB;
                                }
                                if ( mobilitymodel_option == ROLLETT_HUMPHREYS ) {
                                    m = P * mRHHAGB;
                                }

				v = m * ( physConstants.get_halfG_b2() * get_rho( defgid ) - 0.0);
				//if(mynuc==0){cout<<" vNUC= "<<v<< " Gbhalfsq=" <<Gbhalfsq<<" m= "<<m<<"P= "<<P<<endl;}		
				s = v * dt;
				//if(mynuc==0){cout<<" vNUC= "<<v<< " Gbhalfsq=" <<Gbhalfsq<<" m= "<<m<<"P= "<<P<<" s: "<<s<<endl;}		

				//not migrated sufficiently far during dt to reach the deformation grain boundary
				if ( s < dist_seg ) { 
					dt = 0.0;

					dV = ONETHIRD * ico_faces[n_f].area * ( CUBE(nuc_rads[n_f] + s) - CUBE(nuc_rads[n_f]) ); //PH::equal area    //cout<<"dTextVol="<<dTextVol<<endl;
					myCurrentNucTexContrib[defgpool[defgid].ori] -= dV; //all this dV is attributed exclusively growing into one texture component
					myCurrentNucTexContrib[rxgoid] += dV;
					myTotalDeformedTexture[defgpool[defgid].ori] += dV;
					//#####MK::mismatch as a function of nucleus size and deformation structure length scale
					nuc_rads[n_f] += s;
					
					continue;
				}
			
				//overshooting over the deformation grain boundary
				s = dist_seg;
			
				nuc_segids[n_f]++; //reset in new grain
				dt -= s/v;
								
				dV = ONETHIRD * ico_faces[n_f].area * ( CUBE(nuc_rads[n_f] + s) - CUBE(nuc_rads[n_f]) ); //PH::equal area    //cout<<"dTextVol="<<dTextVol<<endl;
				myCurrentNucTexContrib[defgpool[defgid].ori] -= dV;
				myCurrentNucTexContrib[rxgoid] += dV;
				myTotalDeformedTexture[defgpool[defgid].ori] += dV;
				nuc_rads[n_f] += s;
				if(nuc_segids[n_f] == myIcoPaths[mynuc][n_f].nsegments) mysumWallCollisions ++;
					/*
					if(s > dist_seg) {
						s = dist_seg;
						nuc_segids[n_f]++;
						dt -= s/v;
					} else { dt = 0.0; }*/
				
					//double dV = volfactor * ico_faces[n_f].area * ( CUBE(nuc_rads[n_f] + s) - CUBE(nuc_rads[n_f]) ); //PH::equal area    //cout<<"dTextVol="<<dTextVol<<endl;
					//myNucTextureVolumeContribution[defgpool[defgid].ori] -= dV;
					//myNucTextureVolumeContribution[rxgoid] += dV;
					//nuc_rads[n_f] += s;
			} //end while
			
			//collision test:: MK::THE KEY PHYSICAL SIMPLIFICATION OF THE IMPINGEMENT PROBLEM MADE IN THE MODEL IS THAT
			//GRAIN BOUNDARY MIGRATION ALONG THE TRACE n_f CEASES WHEN IT INTERSECTS A SPHERICALLY GROWING NUCLEUS, THIS COULD BE IMPROVED BY TESSELLATING IN MORE DETAIL
			//ALSO THE NEIGHBORS BUT THEN THE FORMAL COMPLEXITY INCREASES AND WE OBTAIN AGAIN THE CELLULAR AUTOMATON APPORACH
		
			for (long nbor = 0; nbor < myNeighbors[mynuc].size(); nbor++ ) {
				dx = myNeighbors[mynuc][nbor].x - (nuc_rads[n_f] * ico_facenormals[n_f].x) - myNuclei[mynuc].x;
				dy = myNeighbors[mynuc][nbor].y - (nuc_rads[n_f] * ico_facenormals[n_f].y) - myNuclei[mynuc].y;
				dz = myNeighbors[mynuc][nbor].z - (nuc_rads[n_f] * ico_facenormals[n_f].z) - myNuclei[mynuc].z;
			
				//assuming the neighbor growing ?
				if( ( SQR(dx) + SQR(dy) + SQR(dz) ) <= SQR(nbor_rads[nbor]) ) {
					//Nucleus collision with neighbor
					nuc_ceased[n_f] = nbor;
					nbor_ceased[nbor] = YES;
					numCeased++;
					break;
				}
			} //test against all neighbors
		} //for all UNCEASED paths
	} //for all faces
}



inline void polyxx::cp_Texture( double *TextSource, double *Text ) {
	long i;
	for( i=0; i < oripool.size(); i++ ){
		Text[i] = TextSource[i];
	}
}

inline void polyxx::add_Texture( double *TextSource, double *Text ) {
	long i;
	for( i=0; i < oripool.size(); i++ ){
		Text[i] += TextSource[i];
	}
}


void polyxx::sim_ICO_PATHGROWTH( void )
{
	double T, neg_kT, t_total, deltat;
	
	//reserve data for tracking local evolution of the nucleus kinetics and texture
	size_t texsize = oripool.size() * sizeof(double);


	//initialize Time-Dependent Values
	for(long myTimeStep = 0; myTimeStep < myTimeAlloc; myTimeStep++) {
		myTotalRXVolumeAllNuclei.push_back( 0.0 );

		double* tt = NULL;
		tt = (double*) calloc(oripool.size() , sizeof(double));
		if ( tt == NULL ) { cout << "ERROR:: allocation of Time-Dependent values myRank;myTimeStep " << myRank << ";" << myTimeStep << endl; return; } 

		myTimeOripoolTextureMatrix.push_back( tt );
	}
	
	maxTimeStepAtPointOfCeasingGrowth = 0;
	minTimeStepAtPointOfCeasingGrowth = LONG_MAX;
	mysumWallCollisions = 0;

        long nPreAllocatedTimeSteps = myTimeAlloc;
cout << "Probe time allocs = " << this->myTimeSteps.size() << ";" << this->myTemperatures.size() << ";" << this->myTotalRXVolumeAllNuclei.size() << ";" << this->myTimeOripoolTextureMatrix.size() << endl;

	//Empty Texture for initialization
	double *ZeroText = NULL;//initial no texture assumed
	ZeroText = (double*) calloc( oripool.size() , sizeof(double)); //Zero
	if ( ZeroText == NULL ) { cout << "ERROR::" << myRank << " failed allocation of ZeroText" << endl; return; }
	
	//Current (per-time-step) Texture change, will be reinitialized and updated each timestep
	double *myCurrentNucTexContrib = NULL;
	myCurrentNucTexContrib = (double*) malloc(texsize);
	if ( myCurrentNucTexContrib == NULL ) { cout << "ERROR::" << myRank << " failed allocation of myCurrentNucTexContrib" << endl; return; }
	
	//Total consumed Deformed Texture over all timesteps, initialize with zero values
	double *myTotalDeformedTexture = NULL;
	myTotalDeformedTexture = (double*) calloc( oripool.size() , sizeof(double)); //Zero
	if ( myTotalDeformedTexture == NULL ) { cout << "ERROR::" << myRank << " failed allocation of myTotalDeformedTexure" << endl; return; }

	long mynuc;

//RUN COMPUTATION OF ALL NUCS SEQUENTIAL ...
	for (mynuc = 0; mynuc < myIDs.size(); mynuc++ ) {

		long numCeased = 0;

		//initialize neighbor paths
		for (long nbor = 0; nbor < myNeighbors[mynuc].size(); nbor++ ) {
			nbor_rads.push_back( 0.0 );
			nbor_segids.push_back( 0 );
			nbor_ceased.push_back( UNCEASED );
		}

		//initialize nuc-paths
		for ( long n_f = 0; n_f < numFaces; n_f++ ) { //numFaces the same for all nuclei because same meshRecursion!
			nuc_rads.push_back( 0.0 );
			nuc_segids.push_back( 0 );
			nuc_ceased.push_back( UNCEASED );
		}

		//initialize texture by counting how much volume of deformed grains
		//init_CalculateDefTexVolumeGrainResolved(mynuc, myNucTextureVolumeContribution);

		long tstep = -1;

		//lets go until growth on all pathes ceased, either because of impingement or wall-collision
		while (numCeased < numFaces) {

			tstep++;

//cout << "Growing " << mynuc << " timestep " << tstep << " myTimeAlloc " << myTimeAlloc << " nPre" << nPreAllocatedTimeSteps << ";" << numCeased << endl;
			cp_Texture(ZeroText, myCurrentNucTexContrib);//MK::zeroText is a clearing mask for myCurrentNuc init Texture Contribution with zero values for every timestep

			if ( tstep >= (nPreAllocatedTimeSteps-1) ) { //extend the precomputed time scheme

				cout << "RANK: " << myRank << " : Out of Range-Error : Time-Reallocation with Alloc=" << myTimeSteps.size() << endl;

				//get new memory to store intermediate data
				init_precomputeMyTimeTemperature( myTimeSteps.back() , myTimeAlloc);
				for( long j = 0; j < myTimeAlloc; j++ ) { 
					myTotalRXVolumeAllNuclei.push_back( 0.0 ); 
				}
				for( long m = 0; m < myTimeAlloc; m++ ) {
					double* textureBucket = NULL;
					textureBucket = (double*) calloc(oripool.size() , sizeof(double));
					if ( textureBucket == NULL ) { cout << "ERROR::Allocation error in myrank;mynuc at TimeOripool 1 " << myRank << ";" << mynuc << ";" << tstep << endl; return; }
					myTimeOripoolTextureMatrix.push_back( textureBucket );
				}

				nPreAllocatedTimeSteps = nPreAllocatedTimeSteps + myTimeAlloc;
				myTimeAlloc = nPreAllocatedTimeSteps;
			}

			//set physical quantities at time = tstep
			t_total = myTimeSteps[tstep];
			deltat = myTimeSteps.at( tstep + 1 ) - t_total; //at necessary to throw potential out of range error
			T = myTemperatures[tstep];

			neg_kT = (NEG_KBOLTZMANN / T);
			mGS = GSm0 * exp(GSHact * neg_kT );
			mLAGB = LAGBm0 * exp( LAGBHact * neg_kT );
			mHAGB = HAGBm0 * exp( HAGBHact * neg_kT );
			mRHHAGB = RHHAGBm0 * exp( RHHAGBHact * neg_kT );
//cout << "mGS/HAGB/LAGB/RHHABG = " << mGS << ";" << mHAGB << ";" << mLAGB << ";" << mRHHAGB << "--" << mobilitymodel_option << endl;
			physConstants.update_Constants(T);

			//all neighbors grow
			sim_ICO_NBORGROWTH_1tstep(mynuc, t_total, deltat ); //, nbor_segids, nbor_rads); //##why nbor_segids and nbor_rads are within the same class?
			
			//all nuclei growth and test for collision in neighbors
			sim_ICO_NUCGROWTH_1tstep(mynuc, t_total, deltat, numCeased, myCurrentNucTexContrib, myTotalDeformedTexture);
			
			//all Nuclei calculated and tested -> next timestep
			//MK::-->export here the volume of the grain at timestep tstep
			myTotalRXVolumeAllNuclei.at(tstep) += sim_CALCICOVOLUME(); //nuc_rads is part of the class!

			//cout << "myTotalRXVolumeAllNuclei = " << myTotalRXVolumeAllNuclei.at(tstep) << endl;
			add_Texture(myCurrentNucTexContrib , myTimeOripoolTextureMatrix.at(tstep));


			//printing necessary?
			if( nuc_printing[mynuc] == 0 ) continue;
			if( tstep % nuc_printing[mynuc] != 0) continue;
			stringstream subDirName , NucPrintPrefix;
			NucPrintPrefix << "RxPAT.GrainPrint."<< simulationid << "." << myRank << "." << WhoHasHowManyNucleiCumulated[myRank] + mynuc;
			subDirName << "GrainPrint_"<< simulationid <<"_NucID_"<< WhoHasHowManyNucleiCumulated[myRank] + mynuc;
			grainPrinter NucleusPrinter ( NucPrintPrefix.str(), myNuclei[mynuc].x, myNuclei[mynuc].y, myNuclei[mynuc].z, nuc_rads[0], numFaces, ico_faces[0], ico_vertices[0], myNeighbors[mynuc].size(), myNeighbors[mynuc][0], nbor_rads[0], nbor_ceased[0], subDirName.str());
			NucleusPrinter.print( tstep );

		} //end of the integration scheme because all paths of the mynuc have ceased 

		mysumIntegrationSteps = mysumIntegrationSteps + tstep;
//cout << mynuc << "\t\t" << numCeased << "\t\t" << tstep << "\t\t" << mysumWallCollisions << endl;


		//print grainshape for last timestep
		if( nuc_printing[mynuc] != 0 ){
			stringstream subDirName , NucPrintPrefix;
			NucPrintPrefix << "RxPAT.GrainPrint."<< simulationid << "." << myRank << "." << WhoHasHowManyNucleiCumulated[myRank] + mynuc;
			subDirName << "GrainPrint_"<< simulationid <<"_NucID_"<< WhoHasHowManyNucleiCumulated[myRank] + mynuc;
			strucPrinter StrucPrinter ( NucPrintPrefix.str(),
					defgmedian_rd, defgmedian_nd, defgmedian_td, &(myBoxes[mynuc]),
					myNuclei[mynuc].x, myNuclei[mynuc].y, myNuclei[mynuc].z,
					nuc_rads[0], numFaces, ico_faces[0], ico_vertices[0],
					myNeighbors[mynuc].size(), myNeighbors[mynuc][0], nbor_rads[0], nbor_ceased[0],
					subDirName.str());
			/*
			for (long ixyz = 0; ixyz < myBoxes[mynuc].ngrxyz; ixyz ++){
			cout <<" i=" << ixyz << " ori= " << defgpool[myDefMS[mynuc][ixyz]].ori <<endl;
			}
			*/
			StrucPrinter.add_OriData(myDefMS[mynuc],defgpool);
			StrucPrinter.print_DefStructure();
			StrucPrinter.print( tstep );
			StrucPrinter.print_Box();
			StrucPrinter.print_nbor_time_evo();

			nucleationPrinter NucleationPrinter(NucPrintPrefix.str(), &(myBoxes[mynuc]), defgpool, myDefMS[mynuc],
				defgmedian_rd, defgmedian_nd, defgmedian_td,
				&(myNuclei[mynuc]),	&(myNeighbors[mynuc][0]), myNeighbors[mynuc].size(), rxgpool ,
				&(oripool[0]), subDirName.str());
			NucleationPrinter.print();
		}

		//empty containers
		nbor_rads.clear();
		nbor_segids.clear();
		nbor_ceased.clear();

		nuc_rads.clear();
		nuc_segids.clear();
		nuc_ceased.clear();

		//MK::the key idea is have all nuclei on one local node run in the time stepping of the node
		//clearly some start earlier than others, ###here we have to correct for tincub distributions, MK:: but not at the moment
		//thus we sample at regular intervals in time the growth of all nuclei and aggregate in the node
		//which can subsequently be aggregated on the master into global statistics

		//Save whole information after nuc-growth totally ceased	
		//save last timestep fur current nuc (for adjustment of kinetic/texture array
		myCeaseTimeSteps.push_back(tstep);
		//find out the minimum/maximum number of timesteps over all myNucs to know which part of the TimeOripoolTextureMatrix is used
		if (tstep > maxTimeStepAtPointOfCeasingGrowth) maxTimeStepAtPointOfCeasingGrowth = tstep;
		if (tstep < minTimeStepAtPointOfCeasingGrowth) minTimeStepAtPointOfCeasingGrowth = tstep;

		//dont forget to log the last timestep, therefore calculate volume of current nuc
		double vol = sim_CALCICOVOLUME();
		//save Volume in Volume List and Distribution, add for calculation of meanVol
		myNucRXFinalVolume.push_back( vol );
		BinMyNucRXFinalVolumeIntoDistr( vol );
		//##PH::fill up coordination number array with 0 value as a dummy currently DEBUG
		myCoordinationNumbers.push_back( 0 );

//cout<<"Recent nuc: "<<mynuc<<" vol = "<<vol<<endl;

cout << "RANK: "<< myRank << " Volume of MyNuc = " << mynuc << " is now = " << vol << " micron^3" << endl;

	} //do so until all mynuc have been calculated

	//MK::end of the simulation, begin node-local postprocessing

	double cumVol = 0.0;
	long tstep;
	//Adjusting RX-Volume for late timesteps, this takes care that the volume of all nuclei that stopped growing at PointOfCeasingGrowth is accounted for even if other nuclei take longer to impingement
	//MK::Therewith it is assumed that grains to not coarsen significantly after they have impinged!
	for (mynuc = 0; mynuc < myIDs.size(); mynuc ++){
		for (tstep = myCeaseTimeSteps[mynuc] + 1; tstep <= maxTimeStepAtPointOfCeasingGrowth; tstep ++){
			myTotalRXVolumeAllNuclei[tstep] += myNucRXFinalVolume[mynuc];
		}
	}
	//CUMULATION OF ALL TEXTURE CHANGES OVER TIME + CONSUMED DEFORMATION TEXTURE
	for(tstep = 0; tstep <= maxTimeStepAtPointOfCeasingGrowth; tstep++ ) { 
		//CUMULATE
		//TimeDependentTextureChange + DeformedTexture -> myTotalDeformedTexture
		add_Texture(myTimeOripoolTextureMatrix[tstep] , myTotalDeformedTexture);
		//SAVE
		//myTotalDeformedTexture -> TimeDependentTextureDevelopment
		cp_Texture(myTotalDeformedTexture, myTimeOripoolTextureMatrix[tstep]);
		//
	} //for all timesteps
	myMaxTime = myTimeSteps[maxTimeStepAtPointOfCeasingGrowth];
	//PH::FREE THE WORLD! FREE WHALES! FREE EUROPE! FREE BEER! (MK::this I cannot agree with - others though: a definite yes!)

	free(ZeroText);
	free(myTotalDeformedTexture);
	free(myCurrentNucTexContrib);

	nbor_segids.clear();
	nbor_rads.clear();	
	nbor_ceased.clear();

	nuc_rads.clear();
	nuc_segids.clear();
	nuc_ceased.clear();
}


void polyxx::init_myGrainSizeDistrBinning( void )
{
	//##MK::requires all domains of same size
	//meanVol = (boxx*boxy*boxz) / (nNeighborsPerCA+1);
        meanVol = (1.0 / ppdensity); //average volume allotted to each nucleus

	//linearized binning from 0 to 
	grvolhist.N_Bins = (long) grvolhist.VolDistMaxValue / grvolhist.VolDistBinSize;
	if ( grvolhist.N_Bins <= 1 ) { cout << "ERROR::Insufficient number of bins!" << endl; return; }

	myGrainVolHistoCounts = (long*) calloc( grvolhist.N_Bins, sizeof(long) );
	myBinEnds.reserve(  grvolhist.N_Bins );

	for(long nbin = 0; nbin < grvolhist.N_Bins; nbin++) {
		myBinEnds.push_back( grvolhist.VolDistBinSize * (nbin+1) );
	}
}


inline void polyxx::BinMyNucRXFinalVolumeIntoDistr( double vol ) {
	long bin = (long) ( vol / ( meanVol * grvolhist.VolDistBinSize) );

	//##MK::what happens when grain is categorized larger than Nbins?
	if(bin < grvolhist.N_Bins ){ 
		myGrainVolHistoCounts[bin]++;
	}
}


void polyxx::pp_GrainVolDistributionData( void )
{
	//malloc changed to calloc!
	ensembleGrainVolHistoCounts = (long*) calloc ( grvolhist.N_Bins, sizeof(long) );	
}


void polyxx::pp_mapMyTimeSteps( void )
{
	dt_rediscrEnsemble = ( AllMaxTime - 0.0 ) / (nRediscrEnsembleRealTimeIntervals - 1); //PH::as I start indexing at 0
	long mytstep_cur = 0;
	double t_global;
	myTimeStepMapping.reserve(nRediscrEnsembleRealTimeIntervals);

	//##what is the attempt here why while and ifelse nested?
	for ( long tstep_g = 0; tstep_g < nRediscrEnsembleRealTimeIntervals; tstep_g++ ) {

		t_global = dt_rediscrEnsemble * tstep_g; //time starts at 0.0

		while( mytstep_cur < maxTimeStepAtPointOfCeasingGrowth && myTimeSteps[mytstep_cur] < t_global ) {
			mytstep_cur ++;
		}
		
		if( mytstep_cur > 0 ) 
			myTimeStepMapping.push_back( mytstep_cur ); //right end of the time interval
		else 
			myTimeStepMapping.push_back( 1 );
	}
}


void polyxx::pp_interpKineticsAndTextureData( void )
{
	myTimeStepDataSize = (1 + oripool.size() ) * nRediscrEnsembleRealTimeIntervals * sizeof(double); //+1 because for summary kinetic as well
	myTimeStepData = NULL;
	myTimeStepData = (double*) calloc( (1 + oripool.size())*nRediscrEnsembleRealTimeIntervals, sizeof(double) );
	if ( myTimeStepData == NULL ) { cout << "ERROR::Allocation problem in pp_interp.." << endl; }

        //MK::arrangement of this array is first nRediscrEnsembleRealTimeIntervals are kinetic data | then all timesteps for the first texture component, then the second and so forth 
	ensembleTimeStepData = (double*) calloc( (1 + oripool.size())*nRediscrEnsembleRealTimeIntervals, sizeof(double) ); //PH::necessary pair wiese existence require MPI_Reduce
	long tpos = 0;
	long mytstep_upper;

	//Push Kinetic Data to array
	for ( long tstep_g = 0; tstep_g < nRediscrEnsembleRealTimeIntervals; tstep_g++ ){
		mytstep_upper = myTimeStepMapping[tstep_g]; //id of the upper time from node-local time data

		if( mytstep_upper > 1 && mytstep_upper < maxTimeStepAtPointOfCeasingGrowth ) { //linear interpolation most likely case
			double changeRate = ( myTotalRXVolumeAllNuclei[mytstep_upper] - myTotalRXVolumeAllNuclei[mytstep_upper-1] );
			changeRate /= ( myTimeSteps[mytstep_upper] - myTimeSteps[mytstep_upper-1] ); ///dt 0 mutually excluded because then nothing would have grown

			double dtt = (dt_rediscrEnsemble * tstep_g) - myTimeSteps[mytstep_upper-1];

			myTimeStepData[tpos] = myTotalRXVolumeAllNuclei[mytstep_upper - 1] + (dtt * changeRate);
                        tpos++;
			continue;
		} 
		//else {
		if ( mytstep_upper == 1 ) 
			myTimeStepData[tpos] = myTotalRXVolumeAllNuclei[0];
		if ( mytstep_upper == maxTimeStepAtPointOfCeasingGrowth ) 
			myTimeStepData[tpos] = myTotalRXVolumeAllNuclei[maxTimeStepAtPointOfCeasingGrowth];
                
                tpos++;
	}
	//first all kinetic data are written then all nRediscr * oripool.size()


	//Push Texture Data to Array
	for (long oid = 0; oid < oripool.size(); oid++ ) { 
		//##MK::is there a way to avoid this? OripoolTextureMatrix is vector<double*> so the array pointed too are spread in the world
		
		for( long tstep_g = 0; tstep_g < nRediscrEnsembleRealTimeIntervals; tstep_g++ ) {
			mytstep_upper = myTimeStepMapping[tstep_g];

			if ( mytstep_upper > 1 && mytstep_upper < maxTimeStepAtPointOfCeasingGrowth ) { //linear interpolation
				double changeRate = (myTimeOripoolTextureMatrix[mytstep_upper][oid] - myTimeOripoolTextureMatrix[mytstep_upper-1][oid]); //provoking conflict misses is like asking for trouble...
				changeRate /= ( myTimeSteps[mytstep_upper] - myTimeSteps[mytstep_upper-1] ); //if you are teasy -[u-1] + [u] is probably more efficient...

				double dtt = (dt_rediscrEnsemble*tstep_g) - myTimeSteps[mytstep_upper - 1];
				myTimeStepData[tpos] = myTimeOripoolTextureMatrix[mytstep_upper-1][oid] + (dtt * changeRate);
				
				tpos++;
				continue;			
			}
			//else { //if outside time-row, set as min/max
			if (mytstep_upper == 1) 
				myTimeStepData[tpos] = myTimeOripoolTextureMatrix[0][oid];
			if (mytstep_upper == maxTimeStepAtPointOfCeasingGrowth) 
				myTimeStepData[tpos] = myTimeOripoolTextureMatrix[maxTimeStepAtPointOfCeasingGrowth][oid];
			//}
			tpos++;
		}
	}
//cout << "Interpolate kinetics tpos/oripoolsize/nRed/(1+ops)*nRed) = " << tpos << ";" << oripool.size() << ";" << this->nRediscrEnsembleRealTimeIntervals << ";" << (1+oripool.size())*this->nRediscrEnsembleRealTimeIntervals << endl;
}


void polyxx::pp_VolList( void )
{
	AllVolumes.resize(nCAEnsemble);
}


void polyxx::out_VolList( void ) {
	stringstream vol_out_fname, vol_out_stream;
	FILE* vol_out_file;

	double TrueTotalVol = 0.;
	for( long n = 0; n < nCAEnsemble; n++ ) {
		TrueTotalVol += AllVolumes[n];
	}
	double TrueMeanVol = TrueTotalVol / nCAEnsemble;

	vol_out_fname << "RxPAT." << simulationid << ".VolList.csv";
	//Now stream dat data
	vol_out_stream << "VolID;Volume[um^3];Sorted Volume[um^3];ln(Vi/Vmean);CumVol[%]"<<";;;;PredictedMeanVol:;"<< meanVol << ";TrueMeanVol:;" << TrueMeanVol << endl;
		
	//Sort Volume List in ascending manner
	vector<double> SortedVolList = AllVolumes;
        
	std::sort(SortedVolList.begin(),SortedVolList.end(),SortDblAscending);
        
	//report cumulated distribution and unsorted list
	double CumVol = 0.0;
	for( long nuc = 0; nuc < nCAEnsemble; nuc++ ) {
		CumVol += SortedVolList[nuc];
		vol_out_stream << nuc << ";" << AllVolumes[nuc] << ";" << SortedVolList[nuc] << ";" << log(SortedVolList[nuc]/TrueMeanVol)<< ";" << (100*CumVol) / TrueTotalVol << endl;
	}

	vol_out_file = fopen( vol_out_fname.str().c_str(),"w" );
	fprintf( vol_out_file,vol_out_stream.str().c_str() );
	fclose( vol_out_file );
}



void polyxx::out_Kinetics( void )
{
	stringstream kinetics_out_fname, kinetics_out_stream;
	FILE* kinetics_out_file;
	kinetics_out_fname << "RxPAT." << simulationid << ".Kinetics.csv";
	kinetics_out_stream << "TimeStep;Time;X" << endl;
	double TotalRXVolume = ensembleTimeStepData[nRediscrEnsembleRealTimeIntervals-1]; //PHMK::ASSUMES THAT THE SYSTEM COMPLETELY RECRYSTALLIZED!

	for( long tstep = 0; tstep < nRediscrEnsembleRealTimeIntervals; tstep++ ) {
		kinetics_out_stream << tstep << ";" << scientific << setprecision(15) << tstep*dt_rediscrEnsemble << ";" << ensembleTimeStepData[tstep]/TotalRXVolume << endl;
	}
	kinetics_out_file = fopen( kinetics_out_fname.str().c_str(),"w" );
	fprintf( kinetics_out_file, kinetics_out_stream.str().c_str() );
	fclose( kinetics_out_file );
}


void polyxx::out_Texture( void )
{
	stringstream texture_out_fname, texture_out_stream;
	FILE* texture_out_file;
	texture_out_fname << "RxPAT." << simulationid << ".Texture.csv";
	texture_out_stream<<"TimeStep;Time";
	double TotalRXVolume = ensembleTimeStepData[nRediscrEnsembleRealTimeIntervals-1];
	long ori;
	/*for (ori = 0; ori < noripool; ori++) {
		texture_out_stream << ";Ori" << ori << setprecision(3) << ": Bunge["<<RADTODEG*oripool[ori].bunge1 << "," << RADTODEG*oripool[ori].bunge2 << "," << RADTODEG*oripool[ori].bunge3 << "]";
	}*/
	//MK::20160216HOTFIX
	texture_out_stream << ";Random[...]";
	for (long s = 0; s < nstandardlagen; s++ ) {
		texture_out_stream << ";Ori" << s << "--Bunge["<< RADTODEG*standardlagen[s].bunge1 << "," << RADTODEG*standardlagen[s].bunge2 << "," << RADTODEG*standardlagen[s].bunge3 << "]";	
	}
		
	
	texture_out_stream << endl;

	double* Texture = ensembleTimeStepData + nRediscrEnsembleRealTimeIntervals; //PH::texture information comes after kinetic information
	long pos = nRediscrEnsembleRealTimeIntervals;
       
	double* buffer = NULL; buffer = new double[nstandardlagen+1];
	long scurr;
	
	for(long tstep = 0; tstep < nRediscrEnsembleRealTimeIntervals; tstep++) {
		texture_out_stream << tstep << ";" << scientific << setprecision(15) << tstep*dt_rediscrEnsemble;
		
		/*//fill buffer
		for (ori = 0; ori < noripool; ori++) {
			pos = nRediscrEnsembleRealTimeIntervals*ori + tstep;
			texture_out_stream << ";" << scientific << setprecision(15) << Texture[pos]/TotalRXVolume;
		}
		texture_out_stream << endl;
		*/
		
		//clean buffer
		for (long st = -1; st < nstandardlagen; st++ ) { 
			buffer[st+1] = 0.0;
		}
		
		for (ori = 0; ori < noripool; ori++ ) {
			pos = nRediscrEnsembleRealTimeIntervals*ori + tstep;
			
			scurr = oripool[ori].closestideal + 1; //[-1, nstandardlagen-1] correct for RANDOM
			
			buffer[scurr] = buffer[scurr] + Texture[pos];
		}
		
		//write out
		for( long stt = -1; stt < nstandardlagen; stt++ ) {
			texture_out_stream << ";" << scientific << setprecision(15) << buffer[stt+1] / TotalRXVolume;
		}
		
		texture_out_stream << endl;
	}
	
	delete [] buffer; buffer = NULL;
	
	texture_out_file = fopen(texture_out_fname.str().c_str(),"w" );
	fprintf( texture_out_file,texture_out_stream.str().c_str() );
	fclose( texture_out_file );
}


void polyxx::out_GrainVolDistribution( void )
{
	stringstream dist_out_fname, dist_out_stream;
	FILE* dist_out_file;
	dist_out_fname << "RxPAT." << simulationid << ".GrainVolDist.csv";
	dist_out_stream << "VolBinEnd;NormBinEnd;NormBinCount;BinCount" << endl;
	for( long bin = 0; bin < grvolhist.N_Bins; bin++ ) {
		dist_out_stream << myBinEnds[bin]*meanVol << ";" << myBinEnds[bin] <<";" << setprecision(8) << scientific << ( (double) ensembleGrainVolHistoCounts[bin] ) / nCAEnsemble << ";" << ensembleGrainVolHistoCounts[bin] << endl;
	}
	dist_out_file = fopen(dist_out_fname.str().c_str(),"w" );
	fprintf( dist_out_file,dist_out_stream.str().c_str() );
	fclose( dist_out_file );
}

void polyxx::report_RealTimeComputingData( void )
{
	if (myRank == MASTER) {
		AllComputingTimes = (double*) malloc (nRanks * myComputeTimeLog.get_entries() * sizeof(double));
		if (AllComputingTimes == NULL){	cout << "ERROR::While trying to get memory for AllComputingTimes." << endl; }

		AllCollisions = (long*) malloc(nRanks * sizeof(long));
		if ( AllCollisions == NULL) { cout << "ERROR::While trying to get memory for AllCollision.s" << endl; }

		AllIntegrationSteps = (long*) malloc(nRanks * sizeof(long));
		if ( AllIntegrationSteps == NULL) { cout << "ERROR::While trying to get memory for AllIntegrationSteps." << endl; }
	}


	MPI_Gather( &myComputeTimeLog.times[0], myComputeTimeLog.get_entries(), MPI_DOUBLE, AllComputingTimes, myComputeTimeLog.get_entries(), MPI_DOUBLE, MASTER, MPI_COMM_WORLD );

	MPI_Gather( &mysumWallCollisions, 1, MPI_LONG, AllCollisions, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );

	MPI_Gather( &mysumIntegrationSteps, 1, MPI_LONG, AllIntegrationSteps, 1, MPI_LONG, MASTER, MPI_COMM_WORLD );


	if (myRank == MASTER) {
		//create file and master appends information via stringstream
		stringstream realtime_log_fname, realtime_log_stream;
		FILE *realtime_log_file = NULL;
		realtime_log_fname << "RxPAT." << simulationid << ".ComputingTimes.log.csv";
		realtime_log_stream << "Rank";

		for (long n = 0; n < myComputeTimeLog.get_entries(); n++) {
			realtime_log_stream << ";" << myComputeTimeLog.titles[n];
		}
		realtime_log_stream << ";mysumWallCollisions;mysumIntegrationSteps" << endl;

		long offset;
		long RTpos;
		for ( int r = 0; r < nRanks; r++) {
			offset = r * myComputeTimeLog.get_entries();
			realtime_log_stream << r << ";" ;
			for (RTpos = 0; RTpos < myComputeTimeLog.get_entries(); RTpos ++){
				realtime_log_stream << AllComputingTimes[offset+RTpos] - realStartTimeOfComputing << ";" ;
			}
			realtime_log_stream << AllCollisions[r] << ";" << AllIntegrationSteps[r] << endl;
		}
		realtime_log_file = fopen (realtime_log_fname.str().c_str(), "w");
		if (realtime_log_file == NULL) { cout <<"ERROR: While trying to write to: " << realtime_log_fname << endl; return; }

		fprintf( realtime_log_file, realtime_log_stream.str().c_str());
		fclose (realtime_log_file);
	}

	if ( myRank == MASTER ) {
		free(AllComputingTimes);
		free(AllCollisions);
		free(AllIntegrationSteps);
	}

}
 

/*
void polyxx::out_Kinetic( void){

}
void polyxx::out_Texture( void ){

}
*/
double polyxx::init_CALCICOISOVOLUME( void )
{
	double vol = 0.0;
	for( long n_f = 0; n_f < numFaces; n_f++ ) {
		vol += ico_faces[n_f].area;
	}
	vol *= ONETHIRD;
	return vol;
}

double polyxx::sim_CALCICOVOLUME ( void ) //vector<double> &nuc_rads )
{
	double vol = 0.0;
	//Volume is calculated from small pyramides
	//V_P=1/3*A*R
	//since facial areas A_0 are calculated on unit sphere
	//A can directly be calculated through eq. A=A_0*R^2 
	//->V_P=1/3*A_0*R^3
	for ( long n_f = 0; n_f < numFaces; n_f++ ) {
		vol += ico_faces[n_f].area * CUBE(nuc_rads[n_f]);
	}
	vol *= ONETHIRD;

	return vol;
}

/*
void polyxx::sim_VORO_CALCVOLUMES (void)
{
	double vol;
	long n_p;
	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++ ) {
		voro::voronoicell *polyhedron;
		polyhedron=new voro::voronoicell;
		polyhedron->init(-0.5E6*boxx,0.5E6*boxx,-0.5E6*boxy,0.5E6*boxy,-0.5E6*boxz,0.5E6*boxz);
		n_p = 0;
		//VORO++: Initialize voronoicell with box sizes in microns and particle at origin
		//boxx,boxy,boxz
		for (long nbor = 0; nbor < nNeighborsPerCA; nbor++){
			if(n_p==polyhedron->p){
				if (polyhedron->max_radius_squared()<(myContacts[mynuc][nbor].dist*myContacts[mynuc][nbor].dist*1E12)){
				//cout<<"nbor= "<<nbor<<endl;
				//cout<<"Polyhedron rad: "<<sqrt(polyhedron->max_radius_squared())<<endl;
				//cout<<"Dist in microns= "<<myContacts[mynuc][nbor].dist*1E6<<endl;
				//cout<<"posx posy posz "<<myContacts[mynuc][nbor].posx*1E6-0.5E6*boxx<<" "<<myContacts[mynuc][nbor].posy*1E6-0.5E6*boxy<<" "<<-myContacts[mynuc][nbor].posz*1E6+0.5E6*boxy<<endl;
				break;
				}
			}
			
			//cout<<"nuc "<<mynuc<<" nbor "<<nbor<<" posx of nbor contact "<<myContacts[mynuc][nbor].posx<<endl;
			//cout<<"nuc "<<mynuc<<" nbor "<<nbor<<" xpos of cut particle in microns "<<myContacts[mynuc][nbor].posx*2E6-1E6*boxx<<endl;
			polyhedron->plane(myContacts[mynuc][nbor].posx*2E6-1E6*boxx,myContacts[mynuc][nbor].posy*2E6-1E6*boxy,myContacts[mynuc][nbor].posz*2E6-1E6*boxz);
			n_p=polyhedron->p;
		}
		//VORO++: Calculate voronoicell volume
		//debug

		//if(mynuc==0){
		//	stringstream polyhedron_fname,polyhedron_m_fname;
		//	polyhedron_fname<<"RXPATHTRACER.Sim." << simulationid << ".Node." << myRank << ".Phedron.pov";
		//	polyhedron_m_fname<<"RXPATHTRACER.Sim." << simulationid << ".Node." << myRank << ".Phedron_Mesh.pov";
		//	polyhedron->draw_pov_mesh( 0.,0.,0., polyhedron_m_fname.str().c_str() );
		//	polyhedron->draw_pov(0.,0.,0.,polyhedron_fname.str().c_str() );
		//}

		vol=1E-18*polyhedron->volume();
		myNucRXFinalVolume.push_back (vol);
		BinMyNucRXFinalVolumeIntoDistr(vol);
		myCoordinationNumbers.push_back ( polyhedron->number_of_faces() );
		//meanVol+=vol;
		meanCoordNum+=polyhedron->number_of_faces();
		//cout<<"polyhedron volume="<<vol<<endl;
		delete polyhedron;
	}
	//meanVol/=myIDs.size();
	//mean_coord_num/=myIDs.size();
	//cout<<"mean volume ="<<mean_vol<<" um^3"<<endl;
}

void polyxx::sim_VORO_PATHGROWTH( void )
{
	/*
	//##the following work in each MPI node is totally independent numerical integration along the PI local paths
	for ( long mynuc = 0; mynuc < this->myIDs.size(); mynuc++ ){

#ifdef DETAILED_PROMPTS
	cout << myRank << ";myRank;startSolving=" << mynuc << endl;
#endif

		//information relevant for all paths
		double telapsed, tincubo, tincubn, deltat; //elapsed real time, and incubation times
		double Temperature; //temperature
		double disto, distn; //traveled from o, //traveled from n, while distn = len - disto holds...
		long segid_o, segid_n;
		double dist_segend_o, dist_segend_n; //distance before boundary changes region
		double pathlen;
		long rxgid, nrxgid; //ids dereferencing from the rxgpool for nucleus and neighbor
		long rxgoid, nrxgoid; //ids dereferencing from the oripool for nucleus and neighbor, can be utilized to get mobilityweight
		double mGS, mHAGB, mLAGB, Po, Pn, mo, mn; //the mobilities
		long defgoid, defgnid; //in which segments are both currently growing?
		double vo, vn, so, sn, to, tn; //how far traveled with v* during deltat?
		double _kT;

		for (long nbor = 0; nbor < this->nNeighborsPerCA; nbor++ ) {

#ifdef DETAILED_PROMPTS
	cout << "\t\t" << "nbor=" << nbor << endl;
#endif

			//initialize the simulation
			rxgid = this->myNuclei[mynuc].rxgid;
			nrxgid = this->myNeighbors[mynuc][nbor].rxgid;

			rxgoid = rxgpool[rxgid].ori;
			nrxgoid = rxgpool[nrxgid].ori;

			tincubo = this->myNuclei[mynuc].tincub;
			tincubn = this->myNeighbors[mynuc][nbor].tincub;

			//start simulating at the lowest incubation time
			telapsed = tincubo;
			if (tincubn < tincubo) { telapsed = tincubn; }

			disto = RCRIT; 
			distn = RCRIT;
			pathlen = myNBorPaths[mynuc][nbor].totallength;
			//##be careful here because growth front might already be somewhere else!
			segid_o = 0;
			defgoid = myNBorPaths[mynuc][nbor].grainregions_outbound[segid_o];
			segid_n = 0;
			defgnid = myNBorPaths[mynuc][nbor].grainregions_inbound[segid_n];

#ifdef DETAILED_PROMPTS
	cout << "\t\t\t" << myRank << ";" << mynuc << ";" << nbor << ";" << tincubo << ";" << tincubn << "myRank;mynuc;nbor;tincubo;tincubn;pathlen;defgoid,defgnid=" << pathlen << ";" << defgoid << ";" << defgnid << endl;
#endif
			//get an initial deltaT
			Temperature = get_temperature( telapsed );
			_kT = (1.0 / (kboltzman * Temperature));
			mGS = GSm0 * exp( -1.0 * GSHact * _kT );
			//###^time-integration scheme debug version currently
			deltat = (minDist / (mGS * Gbhalfsq * rhoMax)); //no drag

			//###mind incubation time, no growth prior the incubation time!

			//simulate now!
			long iter = 0;
			while ( ((disto + EPSILONBALL) < (pathlen - distn + EPSILONBALL)) && ((disto + EPSILONBALL) < pathlen) && ((distn + EPSILONBALL) < pathlen) ) {

#ifdef DETAILED_PROMPTS
	cout << "\t\t\titer;telapsed;Temperature;deltat;disto;distn=" << iter << "\t\t" << telapsed << "\t\t" << Temperature << "\t\t" << deltat << "\t\t" << disto << "\t\t" << distn << endl;
#endif

				//update mobilities
				double _kT = (1.0 / (kboltzman * Temperature));
				mLAGB = LAGBm0 * exp( -1.0 * LAGBHact * _kT );
				mHAGB = HAGBm0 * exp( -1.0 * HAGBHact * _kT );
				mGS = GSm0 * exp( -1.0 * GSHact * _kT );

				//boundary migration through potentially several regions at constant T, first outbound, boundary width of 
				//rediscretize further the time integration of deltat during which the boundary might ...

				//outbound growth
				dist_segend_o = myNBorPaths[mynuc][nbor].regionends_outbound[segid_o] - disto;
				to = deltat;
				while ( (telapsed >= tincubo) && (to > SMALLTIME) && (segid_o < myNBorPaths[mynuc][nbor].nsegments) ) { //dist_segend_o > (BOUNDARYWIDTH_IN_BURGERSV * b) || (to > 0.0) ) {
					defgoid = myNBorPaths[mynuc][nbor].grainregions_outbound[segid_o];
					Po = get_lowtrig_mobilityweight( rxgoid, defgpool[defgoid].ori );
					mo = Po * mGS + (1 - Po) * mHAGB;
					if (Po < 0.0) mo = mLAGB;
					vo = mo * Gbhalfsq * this->get_rho( defgoid ); //, telapsed );
					//attempt largest possible jump of the boundary
					so = vo * to;

					dist_segend_o = myNBorPaths[mynuc][nbor].regionends_outbound[segid_o] - disto;
					//if this would cause the boundary to pass the region, consider only distance to the locla regionend
					if (so > dist_segend_o) {
						so = dist_segend_o;
						segid_o++;
					}

					//otherwise nothing to do because boundary will not touch the regionend with vo during deltat  
					dist_segend_o = dist_segend_o - so;
					to = to - (so/vo); //if so>dis_segend_o then to > SMALLTIME
					disto = disto + so;
#ifdef DETAILED_PROMPTS
	cout << "\t\t\toutbound;Po;mo;vo;to;" << Po << ";" << mo << ";" << vo << ";" << to << endl;
#endif
				}

				//inbound growth
				dist_segend_n = myNBorPaths[mynuc][nbor].regionends_inbound[segid_n] - distn;
				tn = deltat;
				while ( (telapsed >= tincubn) && (tn > SMALLTIME) && (segid_n < myNBorPaths[mynuc][nbor].nsegments) ) { //dist_segend_n > (BOUNDARYWIDTH_IN_BURGERSV * b) || (tn > 0.0) ) {
					defgnid = myNBorPaths[mynuc][nbor].grainregions_inbound[segid_n];
					Pn = get_lowtrig_mobilityweight( nrxgoid, defgpool[defgnid].ori );
					mn = Pn * mGS + (1 - Pn) * mHAGB;
					if (Pn < 0.0) mn = mLAGB;
					vn = mn * Gbhalfsq * this->get_rho( defgnid ); //, telapsed );
					sn = vn * tn; //attempt largest possible jump of the boundary

					//if this so would cause the boundary to pass the region
					dist_segend_n = myNBorPaths[mynuc][nbor].regionends_inbound[segid_n] - distn;
					if (sn > dist_segend_n) {
						sn = dist_segend_n;
						segid_n++;
					}

					//otherwise nothing to do because boundary will not touch the regionend with vo during deltat  
					dist_segend_n = dist_segend_n - sn;
					tn = tn - (sn / vn);
					distn = distn + sn;
#ifdef DETAILED_PROMPTS
	cout << "\t\t\tinbound;Pn;mn;vn;tn;" << Pn << ";" << mn << ";" << vn << ";" << tn << endl;
#endif
				}

				//now that both nuclei have been modelled growing for a time deltat, update time and state variables get new deltat
				telapsed = telapsed + deltat;
				Temperature = get_temperature( telapsed );
				_kT = (1.0 / ( kboltzman * Temperature));
				mGS = GSm0 * exp( -1.0 * GSHact * _kT );
				//###debug version currently
				deltat = (minDist / (mGS * Gbhalfsq * rhoMax)); //no drag

				iter++;
			} //while neither the nucleus nor the neighbor have managed the pathdistance or met in between...
			//PATH HAS BEEN ANALYZED

			//impingement in between or either nucleus or neighbor faster than the other, check RCRIT
			myNBorPaths[mynuc][nbor].timecontact = telapsed;
			if ( (disto + EPSILONBALL) > pathlen ) {
				myNBorPaths[mynuc][nbor].analyzedYet = NUCLEUS_OVERGREW_NEIGHBOR;
				myNBorPaths[mynuc][nbor].distcontact = pathlen;
			}
			else if ( (distn + EPSILONBALL) > pathlen ) {
				myNBorPaths[mynuc][nbor].analyzedYet = NEIGHBOR_OVERGREW_NUCLEUS;
				myNBorPaths[mynuc][nbor].distcontact = RCRIT;
			}
			else {
				myNBorPaths[mynuc][nbor].analyzedYet = IMPINGEMENT_ALONG_GROWTHPATH;
				myNBorPaths[mynuc][nbor].distcontact = disto;
			}

//cout << myRank << ";nuc;path;" << mynuc << ";" << nbor << ";telapsed;contact=" << myNBorPaths[mynuc][nbor].timecontact << ";" << myNBorPaths[mynuc][nbor].distcontact << endl;

		} //for each path
	} //for each my nucleus

	cout << myRank << " all my growth paths simulated." << endl;
}

void polyxx::sim_VORO_PRECONDITIONING( void )
{
	cout << myRank << "-th current memory load is " << this->myMemGuard << " Bytes." << endl;

	//set up the vector that contains all results unsorted then sort the results
	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++ ) {
		//generate struct array of contactinfo by analyzing all my points
		vector<contactinfo> cinfo;
		cinfo.reserve ( this->nNeighborsPerCA );

		double nucx = myNuclei[mynuc].x;
		double nucy = myNuclei[mynuc].y;
		double nucz = myNuclei[mynuc].z;

		for (long nbor = 0; nbor < this->nNeighborsPerCA; nbor++) {
			contactinfo c;

			c.when = myNBorPaths[mynuc][nbor].timecontact;
			c.dist = myNBorPaths[mynuc][nbor].distcontact;
			//position where the boundary made contact; compiler will optimize this as struct my path will very likely cache-local
			//VERBOSE:
			//cout<<"nucx nucy nucz "<<nucx<<" "<<nucy<<" "<<nucz<<endl;
			//
			c.posx = nucx + (myNBorPaths[mynuc][nbor].distcontact * myNBorPaths[mynuc][nbor].direcx); //myNeighbors[mynuc][nbor].x;
			c.posy = nucy + (myNBorPaths[mynuc][nbor].distcontact * myNBorPaths[mynuc][nbor].direcy); //myNeighbors[mynuc][nbor].y;
			c.posz = nucz + (myNBorPaths[mynuc][nbor].distcontact * myNBorPaths[mynuc][nbor].direcz); //myNeighbors[mynuc][nbor].z;

			cinfo.push_back ( c ); //copy constructor
		}

		//sort this array for ascending distance or time
		std::sort( cinfo.begin(), cinfo.end(), SortByDistanceFromNucSite );
		//std::sort( cinfo.begin(), cinfo.end(), SortByContactTime );

		this->myContacts.push_back ( cinfo ); //copy-constructor
	}
	cout << myRank << " all simulation results preconditioned." << endl;


	//now get distribution of average contact distance (taking into account the k-th nearest contact points)
	kmin = 4;
	kmax = nNeighborsPerCA;
	if (kmin < 4) { cout << "kmin too small in PRECONDITIONING!" << endl; }

	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++) {
		vector<double> avr_knear;
		avr_knear.reserve( kmax - kmin + 1);

		for (int howmany = kmin; howmany <= kmax; howmany++) {
			double avr = 0.0;

			for (int k = 0; k < howmany; k++) {
				avr += myContacts[mynuc][k].dist;
			}
			avr = avr / howmany;
			avr_knear.push_back ( avr ); //copy-constructor
		}

		this->myAverageDistance_knearest.push_back ( avr_knear );
	}
}
*/

void polyxx::out_CONTACTPOINTS ( void )
{
	stringstream log_contact_fname;
	ofstream log_contact_file;

	log_contact_fname << "RXPATHTRACER.Sim." << this->simulationid << ".Node." << myRank << ".ContactPoints.csv";
	cout << "File " << log_contact_fname.str().c_str() << " is opened now" << endl;
	log_contact_file.open ( log_contact_fname.str().c_str() );

	//header
	log_contact_file << "MYNUC;MYNUCRXGID;MYNUCTINCUB;NBOR;NBORRXGID;NBORTINCUB;CONTACTPOSX;CONTACTPOSY;CONTACTPOSZ;CONTACTWHERE;CONTACTWHEN;" << endl;

	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++ ) {
		for (long nbor = 0; nbor < myNeighbors[mynuc].size(); nbor++) {
			log_contact_file << mynuc << ";" << myNuclei[mynuc].rxgid << ";" << myNuclei[mynuc].tincub << ";" << nbor << ";" << myNeighbors[mynuc][nbor].rxgid << ";" << myNeighbors[mynuc][nbor].tincub << ";";
			log_contact_file << myContacts[mynuc][nbor].posx << ";" << myContacts[mynuc][nbor].posy << ";" << myContacts[mynuc][nbor].posz << ";" << myContacts[mynuc][nbor].dist << ";" << myContacts[mynuc][nbor].when << ";" << endl;
		}
	}

	log_contact_file.flush();
	log_contact_file.close();


	stringstream log_avr_knear_fname;
	ofstream log_avr_knear_file;

	log_avr_knear_fname << "RXPATHTRACER.Sim." << this->simulationid << ".Node." << myRank << ".KNearestAvRadius.csv";
	cout << "File " << log_avr_knear_fname.str().c_str() << " is opened now" << endl;
	log_avr_knear_file.open ( log_avr_knear_fname.str().c_str() );

	//header
	log_avr_knear_file << "MYNUC;MYNUCRXGID;MYNUCTINCUB;";
	for (int k = kmin; k <= kmax; k++) { log_avr_knear_file << k << ";"; }
	log_avr_knear_file << endl;

	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++ ) {
		log_avr_knear_file << mynuc << ";" << myNuclei[mynuc].rxgid << ";" << myNuclei[mynuc].tincub << ";";
		for (int k = 0; k < (kmax - kmin + 1); k++) {
			log_avr_knear_file << this->myAverageDistance_knearest[mynuc][k] << ";";
		}
		log_avr_knear_file << endl;
	}

	log_avr_knear_file.flush();
	log_avr_knear_file.close();


	double fgeo = (4.0/3.0) * _PI_;
	stringstream log_eqspherevol_knear_fname;
	ofstream log_eqspherevol_knear_file;

	log_eqspherevol_knear_fname << "RXPATHTRACER.Sim." << this->simulationid << ".Node." << myRank << ".KNearestEqvSphereVol.csv";
	cout << "File " << log_eqspherevol_knear_fname.str().c_str() << " is opened now" << endl;
	log_eqspherevol_knear_file.open ( log_eqspherevol_knear_fname.str().c_str() );

	//header
	log_eqspherevol_knear_file << "MYNUC;";
	for (int k = kmin; k <= kmax; k++) { log_eqspherevol_knear_file << k << ";"; }
	log_eqspherevol_knear_file << endl;

	for (long mynuc = 0; mynuc < myIDs.size(); mynuc++ ) {
		log_eqspherevol_knear_file << mynuc << ";";
		for (int k = 0; k < (kmax - kmin + 1); k++) {
			log_eqspherevol_knear_file << (fgeo * pow( this->myAverageDistance_knearest[mynuc][k], 3)) << ";";
		}
		log_eqspherevol_knear_file << endl;
	}

	log_eqspherevol_knear_file.flush();
	log_eqspherevol_knear_file.close();

	if (myRank == MASTER) { cout << myRank << " final output written successfully." << endl; }
}






void polyxx::out_GrainSizeDistribution_voro(void){
	long n_p = 0;
	double mean_coord_num = 0.;
	vector<double> myBinEnds;
	vector<long> myGrainVolHistoCounts;
	
	//For GrainVolDistribution Output
	double myBinSize=0.02;
	cout<<"reserve "<<(long) (10/myBinSize+1)<<endl;
	myBinEnds.reserve ((long) (10/myBinSize+1));//reserve for maximum 10*meanVol
	myGrainVolHistoCounts.reserve((long) (10/myBinSize+1));//reserve for maximum 10*meanVol
	//Initialize GrainVolDistribution Data
	for(long i=0;i<(long) (10/myBinSize+1);i++){
	myBinEnds[i]=myBinSize*(i+1);
	myGrainVolHistoCounts[i]=0;
	}

	stringstream vol_out_fname,vol_out_stream;
	FILE* vol_out_file;

	vol_out_fname << "RXPATHTRACER.Sim." << this->simulationid << ".Node." << myRank << ".VolandCoordNum_List_Voro.csv";
	//print out data to file
	cout<<"Output Volume List File"<<endl;
	vol_out_stream<<"MEAN:;"<<meanVol<<";"<<mean_coord_num<<endl;
	vol_out_stream<<"ID;Volume[um^3];Coord. Number"<<endl;
	long maxEntry=0;
	long bin;
	//Output Vol List Data
	//and generate Distribution Data
	for(long mynuc = 0; mynuc < myIDs.size(); mynuc++){
		bin=(long) (myNucRXFinalVolume[mynuc]/(meanVol*myBinSize));
		myGrainVolHistoCounts[bin]++;
		if(bin>maxEntry) maxEntry=bin;
		//cout<<"mynuc: "<<mynuc<<" vol = "<<myNucRXFinalVolume[mynuc]<<endl;
		//cout<<"bin= "<<bin<<" vol = "<<myNucRXFinalVolume[mynuc]<<" maxEntry="<<maxEntry<<endl;
		if (myNucRXFinalVolume[mynuc]<0.){cout<<"ERROR: neg vol"<<endl;};
		//output Vol and Coord Num List Data
		vol_out_stream<<mynuc<<";"<<myNucRXFinalVolume[mynuc]<<";"<<myCoordinationNumbers[mynuc]<<endl;
	}
	vol_out_file=fopen(vol_out_fname.str().c_str(),"w");
	fprintf(vol_out_file,vol_out_stream.str().c_str());
	fclose(vol_out_file);
	if((maxEntry*myBinSize)<5) maxEntry=(long) 5/myBinSize;
	if(maxEntry>499) maxEntry=499;
	cout<<"Output Dist Data"<<endl;
	//Output Distribution Data
	stringstream dist_out_fname, dist_out_stream;
	FILE* dist_out_file;
	dist_out_fname << "RXPATHTRACER.Sim." << this->simulationid << ".Node." << myRank << ".GrainVolDist.csv";
	dist_out_stream<<"MEAN:| Vol[um^3] | Coord. Num | :;"<<meanVol<<";"<<mean_coord_num<<endl;
	dist_out_stream<<"VolBinEnd;NormBinEnd;NormBinCount;BinCount"<<endl;
	for(bin=0;bin<=maxEntry;bin++){
		dist_out_stream<<myBinEnds[bin]*meanVol<<";"<<myBinEnds[bin]<<";"<<setprecision(8)<<scientific<<((double)myGrainVolHistoCounts[bin])/myIDs.size()<<";"<<myGrainVolHistoCounts[bin]<<endl;
	}
	dist_out_file=fopen(dist_out_fname.str().c_str(),"w");
	fprintf(dist_out_file,dist_out_stream.str().c_str());
	fclose(dist_out_file);
	//Output Kinetic Data
	if(myTotalRXVolumeAllNuclei.size()>0){
		cout<<"Output Kinetic Data"<<endl;
		FILE *kinetic_out_file;
		stringstream kinetic_fname, kinetic_out_stream;
		kinetic_fname << "RXPATHTRACER.Sim." << this->simulationid << ".Node." << myRank << ".Kinetic.csv";
		kinetic_out_stream<<"TimeStep ; RX-Volume-Fraction ; RX-Volume"<<endl;
		for(long myTimeStep=0; myTimeStep <= maxTimeStepAtPointOfCeasingGrowth; myTimeStep++){
			kinetic_out_stream<<myTimeStep<<" ; "<<myTotalRXVolumeAllNuclei[myTimeStep]/myTotalRXVolumeAllNuclei[maxTimeStepAtPointOfCeasingGrowth]<<" ; "<<myTotalRXVolumeAllNuclei[myTimeStep]<<endl;
		}
		kinetic_out_file=fopen(kinetic_fname.str().c_str(),"w");
		fprintf(kinetic_out_file,kinetic_out_stream.str().c_str());
		fclose(kinetic_out_file);
	}
	//Output Texture Data
	long i;
	if(myTimeOripoolTextureMatrix.size()>0){
		cout<<"Output Texture Data"<<endl;
		FILE *texture_out_file;
		stringstream texture_fname, texture_out_stream;
		texture_fname << "RXPATHTRACER.Sim." << this->simulationid << ".Node." << myRank << ".Texture.csv";
		texture_out_stream<<"TimeStep ";
		for(i=0; i<oripool.size(); i++){
		texture_out_stream<<"; Ori-Vol["<<i<<"] [um]^3";
		}
		texture_out_stream<<endl<<"Bunge-Euler-Angles (phi1,PHI,phi2) [Deg]";
		for(i=0; i<oripool.size(); i++){
		texture_out_stream<<setprecision(2)<<fixed<<";("<<RADTODEG*oripool[i].bunge1<<", "<<RADTODEG*oripool[i].bunge2<<", "<<RADTODEG*oripool[i].bunge3<<")";
		}
		texture_out_stream<<endl;
		for(long myTimeStep=0; myTimeStep <= maxTimeStepAtPointOfCeasingGrowth; myTimeStep++){
			texture_out_stream<<myTimeStep;
			for(i=0; i<oripool.size(); i++){
				texture_out_stream<<setprecision(8)<<scientific<<" ; "<<1E18*myTimeOripoolTextureMatrix[myTimeStep][i];
			}
			texture_out_stream<<endl;
		}
		texture_out_file=fopen(texture_fname.str().c_str(),"w");
		fprintf(texture_out_file,texture_out_stream.str().c_str());
		fclose(texture_out_file);
	}
}
//auxiliaries


inline double polyxx::get_temperature( double t )
{
	//isothermal only?
	if ( temperature_isothermalonly == true ) { return Tiso; }

	//sight, more complex processing desired...
	double currT;
	long i = 0;

//cout << this->time[0] << ";" << this->temperature[0] << endl;
//cout << this->time[1] << ";" << this->temperature[1] << endl;
//cout << "i=" << i << "\ttime[i]=" << time[i] << "\tt=" << t << endl;

	//scan vector t and perform linear interpolation to get Ti(ti)
	int ntimetempstamps = time.size(); //these parameter are never changed in size!

	while ( time[i] < t && i < ntimetempstamps ) {
		i++;
	}

//cout << "afterwhile i is=" << i << endl;

	if (i == 0) {
//cout << "temperature[0]=" << temperature[0] << endl;
		return temperature[0];
	}

	if (i >= ntimetempstamps) {
//cout << "temperature[ntimetempstamps-1]=" << temperature[ntimetempstamps-1] << endl;
		return temperature[ntimetempstamps-1];
	}

	currT = temperature[i-1] + ( (temperature[i] - temperature[i-1]) / (time[i] - time[i-1]) ) * (t - time[i-1]);
//cout << "currT=" << currT << endl;

	return currT;
}


/*
inline void polyxx::update_intrinsicmobilities( void )
{
	mLAGB = LAGBm0 * exp (-1 * LAGBHact / (kboltzman * Temperature));
	mHAGB = HAGBm0 * exp (-1 * HAGBHact / (kboltzman * Temperature));
	mGS = GSm0 * exp (-1 * GSHact / (kboltzman * Temperature));
}
*/




inline double polyxx::get_myNuclei_mobilityweight ( long mynuc, long defgori_id )
{
	return myCA_MobWeightHash[(mynuc * oripool_firstdisjunct_rxg) + defgori_id];
}


inline double polyxx::get_lowtrig_mobilityweight ( long row , long col )
{
    
    long idx = (( 0.5 * row * (row + 1) ) + col);
    //would require row >= col but symmetric matrix so...
    
    if (col > row) {
       idx = (( 0.5 * col * (col + 1) ) + row);
    }
    
    return LowTrigMobWeightHash[idx];
}


/*MK::DEPRECATED
inline double polyxx::get_intrinsicmobility ( double mobWeight )
{
	if (mobWeight < 0.0) { return mLAGB; }
	//else
	return ( ( mobWeight * mGS ) + ( (1.0 - mobWeight) * mHAGB ) );
}
*/

inline double polyxx::get_rho ( long defgid ) //, double t )
{
	double rho = defgpool[defgid].rho0;

/*	double rho = defgpool[defgid].rho0;

	if ( kuhlmann_consider == true) {
		if ( t < defgpool[defgid].kuhlmann_tcrit ) { rho = defgpool[defgid].rho0 - (kuhlmann_alpha * log ( 1 + kuhlmann_beta * t ) ); }
		else { rho = KUHLMANNFINALRHO; }
	}
//cout << "get_rho::t\t\t" << t << "\t\t" << rho << endl; */

	return rho;
}


inline double polyxx::get_zener (double t) //currently done for each cell, very high formal complexity
{
	//###utilize ZenerForceConstant for a quick constant Zener drag
	if (zener_consider == true) {
		long i = 0;
		int ntimezdragstamps = zenertime.size();

		//scan vector t and perform linear interpolation to get pzi(ti)

		while (zenertime[i] < t && i < ntimezdragstamps) {
			i++;
		}
	//cout << "afterwhile i is=" << i << endl;

		if (i == 0) {
	//cout << "zenerforce[0]=" << zenerforce[0] << endl;
			return zenerforce[0];
		}

		if (i >= ntimezdragstamps) {
			//cout << "zenerforce[ntimezdragstamps-1]=" << zenerforce[ntimezdragstamps-1] << endl;
			return zenerforce[ntimezdragstamps-1];
		}

		double pzener = zenerforce[i-1] + ( (zenerforce[i] - zenerforce[i-1]) / (zenertime[i] - zenertime[i-1]) ) * (t - zenertime[i-1]);
	//cout << "pzener=" << pzener << endl;
		return pzener;
	}

	//no Zener drag
	return 0.0;
}


//obsolete because paths have integrator accuracy defined
/*
double polyxx::get_tmax_instantslope_cells ( double t )
{
	double tmin = INFINITE;

	//maximum intrinsic mobility
	//double m_max = GSm0 * exp( -1 * GSHact / (kboltzman * get_temperature(t)) );
	
	//current not absolute somewhen, as we can utilize explicit time-synchronization
	double m_max = GSm0 * exp( -1 * GSHact / (kboltzman * Temperature) ); 

	//maximum migration velocity
	//double v_max = xdefgpool_max * m_max * ((Gbhalfsq * rho_max) - 0.0); //##MK::201402PARA, even the fastest boundary will somewhen feel that the strongest component ceases
	double v_max = m_max * ((Gbhalfsq * rho_max) - 0.0); //##MK::201402PARA, even the fastest boundary will somewhen feel that the strongest component ceases


	//only a fraction of a cell in a time step
	tmin = ((MAXFILLPERSTEP * dcell) / v_max);

#ifdef DETAILED_PROMPTS
cout << "Temp " << Temperature << ";m_max " << m_max << ";v_max " << v_max << ";tmin " << tmin << endl;
#endif

	return tmin;
}
*/




double polyxx::get_tmax_instantslope_temp ( double t )
{
	if ( temperature_isothermalonly == true ) {
		return (INFINITE);
	}

	long i = 0; //scan vector t and get current linear heating rate
	long ntimetempstamps = time.size();

	//MK::reasonable disjunctness of time[i] is guaranteed to assure no occurence of division by zero error, as at one physical time, there can not be two temperatures!

	while (time[i] < t && i < ntimetempstamps) {
		i++;
	}
//cout << "afterwhile i is=" << i << "\tnstamps=" << ntimetempstamps << endl;

	if ( i == 0 ) 
		return (INITIAL_DELTAT);

	if (i < ntimetempstamps ) {
		double dT = ( fabs(temperature[i] - temperature[i-1]) / (time[i] - time[i-1]) ); 
		
		if (dT <= SMALL_HEATRATE) dT = SMALL_HEATRATE; //avoid division by zero error

		return (SMALL_HEAT / dT);
	}

	//if not already returned, then processing schedule ends, temperature stays constant with last value
	return (INFINITE);
}


double polyxx::get_tmax_instantslope_rho ( double t ) //unsigned short defgid )
{
	//if recovery is according to kuhlmann-wilsdorf approach rho(t) = rho(t=0) - alpha*ln(1+beta*t) and rho(t>t(rho=rhocrit)) = rhocrit
	//if t < tcrit instantaneous slope is known | drho/dt | = alpha*beta / (1 + beta * t)
	//t >= tcrit slope constant tmax = INFINITE
	//tcrit = (exp( (rhocrit - rho0)/-alpha ) - 1 ) / beta if alpha > 0.0 and beta > 0.0

	if ( kuhlmann_consider == true ) {
		//if ( t < defgpool[defgid].kuhlmann_tcrit ) {
		double dRHO = kuhlmann_alpha * kuhlmann_beta; //divisor >= 1.0
			
		if (dRHO <= SMALL_RECOVERYRATE) {dRHO = SMALL_RECOVERYRATE;} //avoid division by zero error

		return (SMALL_RHO / dRHO);
	} //else { return INFINITE; } //because no further change
	
	return (INFINITE);
}


double polyxx::get_tmax_instantslope_zener ( double t )
{
	//MK::reasonable disjunctness of time[i] is guaranteed to assure no occurence of division by zero error, 
	//as at one physical time, there can not be two different Zener drags, as bifuraction along breakaway is not modeled... :)!

	if (zener_consider == true) {
		long i;
		int ntimezdragstamps = zenertime.size();

		//scan vector t and get current linear rate of change
		while (zenertime[i] < t && i < ntimezdragstamps) {
			i++;
		}
		//cout << "afterwhile i is=" << i << "\tntimezdragstamps=" << ntimezdragstamps << endl;

		if ( i == 0 ) 
			return (INITIAL_DELTAT);

		if ( i < ntimezdragstamps ) {
			double dPZENER = ( fabs(zenerforce[i] - zenerforce[i-1]) / (time[i] - time[i-1]) );
			if (dPZENER <= SMALL_DRAGGINGRATE) dPZENER = SMALL_DRAGGINGRATE;

			return ( SMALL_ZENERFORCE / dPZENER );
		}
	} //else, when processing schedule ends, or Zener force is not even considered zenerforce stays constant with last value

	return (INFINITE);
}

//consider as obsolete and mount into PATH functions
/*
double polyxx::get_minimum_tmax ( double t ) //, unsigned short defgid )
{
	double min_tmax = INFINITE;
	double tmax_model = INFINITE; //analyze all physical mechanisms that need proper time discretization during the integration scheme

	tmax_model = get_tmax_instantslope_temp ( t );
	if ( tmax_model < min_tmax ) 
		min_tmax = tmax_model;

	tmax_model = get_tmax_instantslope_rho ( t );
	if ( tmax_model < min_tmax ) 
		min_tmax = tmax_model;

	tmax_model = get_tmax_instantslope_zener ( t );
	if ( tmax_model < min_tmax ) 
		min_tmax = tmax_model;
	
	//at maximum driving force with evaluate cells discretization
	//tmax_model = get_tmax_instantslope_cells ( t );
	//if ( tmax_model < min_tmax )
	//	min_tmax = tmax_model;


	return min_tmax;
}
*/
