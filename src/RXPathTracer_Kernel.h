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
 * Author: Markus Kuehbach, Paul Hoffrogge, markus.kuehbach at rwth-aachen.de
 * Revised code basis 20.05.2015
 */

#ifndef RXPATHTRACER_KERNEL_H
#define	RXPATHTRACER_KERNEL_H

#include "RXPathTracer_Defs.h"
#include "RXPathTracer_Types.h"
#include "RXPathTracer_TimeLogger.h"
#include "RXPathTracer_PhysicalConstants.h"
#include "RXPathTracer_GrainPrinter.h"
#include "RXPathTracer_NucleationPrinter.h"
//PH::vorosrc support added here
#include "./vorosrc/src/voro++.hh"

//...

class polyxx : public io, public randomClass, public mathMethods
{
public:
	polyxx();
	~polyxx();
        void specificallyDisorientednewOriFromReference( double *bunge, double sigma_rayl_dis2bunge, double *newOri );
        
	long check_disjunctness ( double * bunge );		//a new orientation that is not in the oripool already

	void init_adjustPRNGs( void );
	void init_MPIDatatypes( void );
	bool init_readparameter( const char * inputfname );
	bool read_parameter ( void );
	bool read_processing ( void );
	bool read_zener ( void );
	bool read_defgpool (void );
	bool read_rxgpool ( void );
	bool read_nucprintingschedule( void);

	void worker_interp_ideal( void );
	void worker_interp_defgpool ( void );
	void worker_interp_rxgpool ( void );
	void worker_interp_oripool( void );
	void worker_interp_processing( void );
	void worker_interp_zener( void );
	void worker_interp_printing( void );

	bool MasterParameterBcast( void );
	bool AllPrepareForInputBcast( void);
	bool WorkerInterpretInputBcast( long * elements );
	
	/*//PH::MPI DATA ROUTINES...
	bool rw_LongStorage ( long p_ID, long &value, bool init);
	bool rw_DoubleStorage ( long p_ID, double &value, bool init);
	
	void setup_DataArraysfor_MPI_Send( void );
	void setup_ParameterArraysfor_MPI_Send( void );
	void setup_ParameterArraysfor_MPI_Recv( void );
	void init_Parameters_from_MPI_data( void );
	void init_DynamicArrays_from_MPI_data( void );
	void rw_MPIDoubleParameter(bool init);
	void rw_MPILongParameter(bool init);

	void rw_MPIDoubleData( double &value, bool init);
	void rw_MPILongData(long &value, bool init);

	void init_rwMPIData(bool init);
	void init_rwStandardlagenMPIData(bool init);
	void init_rwTemperatureMPIData(bool init);
	void init_rwZenerMPIData(bool init);
	void init_rwDefgPoolMPIData(bool init);
	void init_rwRxgPoolMPIData(bool init);
	void init_rwOripoolMPIData(bool init);
	//...MPI DATA ROUTINES
	*/

	void init_someData( void );
	void init_myvectors( long wpieces );
	void init_workpartitioning_myNuclei ( void );
	void init_poissonpointprocess( void );
        double DisoriTwoGrains ( long o1, long o2 );
        long gbfacenucleationmodel( bool onlyat_hagb, long dg1, long dg2, double scale_dens2lamb, double scale_nuc2numb, double bndarea_micron );
        void init_gbface_nucleation_yz( bool onlyhagb );
        void init_gbface_nucleation_xz( bool onlyhagb );
        void init_gbface_nucleation_xy( bool onlyhagb );
        void init_large_defms_nucleation( void );
        void init_large_defms( void );
        
	void init_myBoxes ( void ); //##MK::take care that internal buffers are taken care of
        void init_myNuclei_BasedOnLargeDefMS( void );
        void init_myNuclei ( void );

	void init_myNeighbors ( long towhichnuc );
	void init_myNeighbors_BasedOnPointProcess ( long towhichnuc );
        void init_myNeighbors_BasedOnLargeDefMS( long towhichnuc, long rid ); //rid is largedefmsnuc reference id to get coordinates
	void init_myNuclei_nucleation ( void );
        void init_myNuclei_nucleation_BasedOnLargeDefMS( long towhichnuc, long refid ); //refid is a reference id that details which nucleus from largedefmsnuc was taken for a particular box

	void init_disori_lowertriangle ( void );
	void init_checkTimeSteppingOptimPotential( void );
	void init_MyTimeSteps( void );
	void init_precomputeMyTimeTemperature(double t_0, long nsteps);
	void init_myGrainSizeDistrBinning( void );
	void init_disori_fast ( void ); //uses myCA_MobWeightHash

	void sim_pathsegmentation ( pathP pl, long mn, long nb ); //performs the raytracing through the inhomogeneous 3D cuboid deformed grain aggregate from nuclation site mynuc to neighbor site nbor
	void sim_pathsegmentation_nbors( pathP pl, long mn ); //, long nb, double nxx, double nyy, double nzz );
	void sim_pathsegmentation_nucleus( pathP pl, long mn );
	void sim_ICO_HEDRONMESHING ( void );
	void sim_ICO_NUCPATHCONSTRUCTION ( void );
	void sim_ICO_NBORPATHCONSTRUCTION ( void );
	void sim_ICO_PATHGROWTH ( void );
	void sim_ICO_NBORGROWTH_1tstep( long &mynuc, double &t_total, double &deltat ); // vector<long> &nbor_segids, vector<double> &nbor_rads);
	//void sim_ICO_NUCGROWTH_1tstep( long &mynuc, double &t_total, double &deltat, long &numCeased, vector<double> &nuc_rads, vector<double> &nbor_rads, vector<long> &nuc_ceased, double *myNucTextureVolumeContribution, double volfactor);
	void sim_ICO_NUCGROWTH_1tstep( long &mynuc, double &t_total, double &deltat, long &numCeased, double *myNucTextureVolumeContribution, double* myTotalDeformedTexture);
	
	double sim_CALCICOVOLUME( void ); //vector<double> &nuc_rads);
	inline void BinMyNucRXFinalVolumeIntoDistr( double vol );

	//##MK::what is pp?

	//##MK
	void pp_mapMyTimeSteps( void );
	void pp_interpKineticsAndTextureData( void );
	void pp_GrainVolDistributionData( void );
	void pp_VolList( void );
	//ICO
	double init_CALCICOISOVOLUME( void );//calculates polyvol (V of icosphere for R=1)
	void init_Zero_Texture( double *Text);
	inline void cp_Texture(double *TextSource, double *Text);
	inline void add_Texture(double *TextSource, double *Text);
	void init_CalculateDefTexVolumeGrainResolved ( long mynuc , double *myNucTextureVolumeContribution);
	void print_pov_nucshape ( stringstream &fstream );
	void print_pov_nborshape (long nuc, stringstream &fstream );
	void output_grainshape( long nuc, long paint_id);
	//
	/*
	void sim_VORO_CALCVOLUMES ( void );
	void sim_VORO_PATHCONSTRUCTION ( void );
	void sim_VORO_PATHGROWTH ( void );
	void sim_VORO_PRECONDITIONING ( void ); //elaborated guess which convex hulls to analyze
	*/
	
	//IO of results
	void out_CONTACTPOINTS ( void );
	void out_GrainSizeDistribution_voro( void );
	void out_VolList( void ) ;
	void out_GrainVolDistribution( void );
	void out_Kinetics ( void );
	void out_Texture( void );
	//

	void report_RealTimeComputingData ( void );

	//auxiliaries
	double get_temperature ( double t );
	double get_lowtrig_mobilityweight ( long row , long col );
	double get_myNuclei_mobilityweight ( long mynuc, long defgori_id );
	//double get_intrinsicmobility ( double mobWeight );
	double get_rho ( long defgid ); //, double t );	//process modeling can not afford to always recover the grains completely
	double get_zener ( double t );
	//void update_intrinsicmobilities( void );
	//optimizing numerical integration scheme
	//double get_tmax_instantslope_cells( double t );
	double get_tmax_instantslope_temp ( double t );
	double get_tmax_instantslope_rho ( double t ); // unsigned short defgid );
	double get_tmax_instantslope_zener ( double t );
	//double get_minimum_tmax ( double t ); //, unsigned short defgid );
	long get_closest_standardlage ( double * quat );
	//

	//MPI datatypes
	MPI_Datatype MPI_IO_IdealOri_Type;
	MPI_Datatype MPI_IO_Defgrain_Type;
	MPI_Datatype MPI_IO_Rxgrain_Type;
	MPI_Datatype MPI_IO_Ori_Type;
	MPI_Datatype MPI_IO_Temperature_Type;
	MPI_Datatype MPI_IO_ZenerDrag_Type;
	MPI_Datatype MPI_IO_PrintSchedule_Type;

	MPI_IO_IdealOri* bcast_ideal;
	MPI_IO_Defgrain* bcast_dg;
	MPI_IO_Rxgrain* bcast_rxg;
	MPI_IO_Ori* bcast_ori;
	MPI_IO_Temperature* bcast_Tt;
	MPI_IO_ZenerDrag* bcast_zener;
	MPI_IO_PrintSchedule* bcast_nucprintingschedule;
	//
	timeLogger myComputeTimeLog;
	double realStartTimeOfComputing;
	double* AllComputingTimes;
	long* AllCollisions;
	long* AllIntegrationSteps;
	//##consider bundle in property struct of the class
	long modelType;
	long meshRecursion;
	long grainPrintingEnabled;

	idealP standardlagen;			//the world famous texture components
	long nstandardlagen;

	vector<ori> oripool;			//all polyxx know all orientations, but however this array is worked is worked on at runtime so better leave it up staying
	long noripool;

	//the grains and their properties
	defgP defgpool;					//deformed grains either from GIA or CPFEM
	long ndefgpool;
	rxgP rxgpool;					//rxTexture components
	long nrxgpool;

	double* LowTrigMobWeightHash;
	double* myCA_MobWeightHash;		//low memory overhead version of the lowertriangle matrix with myCAsomany * oripool_firstdisjunct_rxg elements, addressing via [(myCA* oripool_firstdisjunct_rxg)+oriid], oriid is always < oripool_firstdisjunct_rxg
	long oripool_firstdisjunct_rxg; /*marks the interval of items within the oripool array of disjunct ids 0, ..., oripool..._rxg - 1 for which myCA_MobWeightHash values should be cached*/
	
	long nTemp, nZener;

	/*
	//MPI DATA transfer from MASTER to the workers
	//Parameter Arrays:
	double *DoubleParameter;
	long *LongParameter;
	//Data Container:
	long nDoubleDataList, nLongDataList;
	vector <double> MPI_DoubleDatav;
	vector <long> MPI_LongDatav;

	double *MPI_DoubleP;
	long *MPI_LongP;
	//...MPI DATA;
	*/

	//work allocation should be dynamic can also change maybe if MPI job repartitioning is desired so better use a vector
	vector<long> myIDs;			//nuclei IDs! running from 0 to NumberNuclei - 1 global IDs, is only necessary on node level
	vector<int> WhoHasWhichNucleus;	//MPI rank IDs! so that all nodes know there are other CA executed on MPI nodes
	int* WhoHasHowManyNuclei;
	int* WhoHasHowManyNucleiCumulated;// For MPI_Gatherv Call
	
	long numberOfNucsToPrint;
	vector<int> NucIDsforPrinting;
	vector<int> NucPrintingPeriod;
	vector<int> nuc_printing;
	//maybe add how many nuclei each node got
	//ICOSPHERE DATA
	long numFaces;
	double polyvol;							//Icosphere-Volume for R=1, used for normalization on sphere volume
	vector<vertex> ico_vertices;
	vector<face> ico_faces;
	vector<vertex_l> ico_facenormals;
	vector< vector<path> > myIcoPaths;

	vector<double> nbor_rads;
	vector<long> nbor_segids, nuc_segids;
	//##MK::access to these is currently not thread safe!!!
	vector<double> nuc_rads;				//##Growth length for every face-path of active Nucleus is always overwritten!
	vector<long> nuc_ceased;				//##Identifier for every face-path which says whether it collided or not overwritten could also be bool!
	vector<long> nbor_ceased;				//Identifier for every Nbor which says whether it collided with the nucleus or not
	vector<long> myCeaseTimeSteps;			//Time step for every Nucleus, where growth ceased for all paths
	vector<double> myTimeSteps;
	vector<double> myTemperatures;
	vector<long> myTimeStepMapping;
	vector<double> myTotalRXVolumeAllNuclei;				//Total Recrystallized Volume (absolute), summed over all Nuclei for every timestep pyramid volume
	vector< double* > myTimeOripoolTextureMatrix;			//myTimeOripoolTextureMatrix[Time-Step][Tex-Comp] : Volume-Texture information for any timestep, sum over all nucs
	vector< double* > myLastTimeStepOripoolTexture;			//myLastTimeStepOripoolTexture[Nuc-ID][Tex-Comp] : Volume-Texture information for any nuc and last timestep
	
	
	long nRediscrEnsembleRealTimeIntervals;					//Time Output Accuracy
	double dt_rediscrEnsemble;
	size_t myTimeStepDataSize;
	double* myTimeStepData;
	double* ensembleTimeStepData;
	//
	//ICOSPHERE DATA
	double meanVol, meanCoordNum;

	//double VolDistBinSize, VolDistMaxValue;
	//long N_Bins;
	struct grvolbinning grvolhist;
	long* myGrainVolHistoCounts;
	long* ensembleGrainVolHistoCounts;
	vector <double> myBinEnds;

        vector<point> pointprocess;                     //a large point process to sample from
        vector<defglean> largedefms;                    //a large 3D cuboid aggregate of deformed grains to sample from
        struct defmsinfo largedefmsinfo;                //meta data regarding the structure
        vector<rxgo> largedefmsnuc;
        vector<nucsite> myNuclei;			//the nuclei which are growing in each domain
	vector<box> myBoxes;				//the local environment where growth is simulated
	vector<long*> myDefMS;				//implicit 1D array storing defgpool information
	vector< vector<nucsite> > myNeighbors; //the neighborhood
	vector< vector<path> > myNBorPaths;		//all paths, so the work that is done [myCA][path] pointer to an array of all the different pathes that are taken into account for the nucsite-th nucleus simulated
	vector< vector<contactinfo> > myContacts; //sorted list of all contact information for each nucleus
	vector< vector<double> > myAverageDistance_knearest;
	vector<rxgo*> largedefmsnuc_tree; //an index list guiding the search for neighbors
	vector<double> myNucRXFinalVolume;
	vector<double> AllVolumes;
	vector<int> myCoordinationNumbers;

	unsigned short* largedefmsnuc_treebins; //counts how many points happen to fall in one bin

	//double myPmax, myPmin;//Max/Min P-Value calculated for current Rank, utilized for adjusting timestep-size
	//double myRhomax;//Rhomax computed for current node
	
	//truct timings time_rediscr;
	long myTimeAlloc, maxTimeStepAtPointOfCeasingGrowth, minTimeStepAtPointOfCeasingGrowth;
	double myMaxTime, AllMaxTime;

	int kmin;
	int kmax;

	//simulation meta data
	long simulationid;
	long nCAEnsemble;				//so many nuclei simulated all together
	//long nNeighborsPerCA;			//MK::debug - so many neighbors placed in the vicinity of each, was like that! now each domain can have a different number of neighbors!
	//double Vref;					//macroscopic RVE volume to compare against
	double ppdensity;
        double ppscaling;                               //edge length of the RVE in micron on which a point process is initialized
        double boxx;					//debug in micron
	double boxy;
	double boxz;

	double minDist;					//minimum distance swept for numerical integration along the path
	double rhoMax; 

	//processing scheme
	vector<double> time;                //in (s)
	vector<double> temperature;         //in (K)
	double Tiso;

	//a priori coupled information on microchemistry
	vector<double> zenertime;			//in (s)
	vector<double> zenerforce;			//zenerfac * f/r in (J/m^2 * 1/m)
	double ZenerForceConstant;

	//deformation structure
	double defgmedian_rd;
	double defgmedian_nd;
	double defgmedian_td;

	//physical properties
	physicalConstants physConstants;
	double LAGBm0;
	double LAGBHact;
	double HAGBm0;
	double HAGBHact;
	double GSm0;
	double GSHact;
        
        double RHHAGBm0;
        double RHHAGBHact;
        double RHLAGBHAGBCut;
        double RHLAGBHAGBTrans;
        double RHLAGBHAGBExponent;

	//updated intrinsic mobilities
	double mLAGB;
	double mHAGB;
	double mGS;
        double mRHHAGB;
	struct optimtime optTimeStepping;


	double kuhlmann_alpha;		//dislocation density for isothermal recovery
	double kuhlmann_beta;		//time-dependent factor for isothermal recovery

	//physical caching variables
	double zenerfac;			//fgeo * gammaGB


	//further options
	long rvoption, zeneroption, quickisothermal_option, mobilitymodel_option; //PH: MPI_INT or MPI_LONG necessary for MPI Support
	bool temperature_isothermalonly;
	bool kuhlmann_consider;
	bool zener_consider;

	long binarymatrix_outputoption;

	//Parallel environment relevant information, once a simulation starts we do not further change the assignment of jobs
	int myRank;
	int nRanks;

	double myMemGuard;				//how much memory is currently consumed in bytes in this node
	
	long mysumWallCollisions;
	long mysumIntegrationSteps;
	long allsumWallCollisions;
	long allsumIntegrationSteps;

	bool mySuccess;
	randomClass rndlocal;
        randomClass rndlargedefms;
        
	FILE * rxpathtracer_input;
};
typedef class polyxx * polyxxP;


#endif

