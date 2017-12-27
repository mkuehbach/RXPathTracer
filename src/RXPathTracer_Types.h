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


#ifndef RXPATHT_TYPES_H
#define RXPATHT_TYPES_H

#include <iomanip>

struct point
{
    double x;
    double y;
    double z;
};
typedef point * pointP;


struct vertex
{
	double x;
	double y;
	double z;
};

struct face
{
	long A;				//three vertice -indices
	long B;
	long C;				
	double area;		//facial area of the pyramidal base triangle with height r = 1
};

struct nucsite
{
	double x;
	double y;
	double z;
	double tincub;		//there should be the possibility to have different incubation times for the same rxgOrientation
	//long defgid;
	long rxgid;
};
typedef struct nucsite * nucsiteP;


struct regionpair
{
	long dgid;
	double regionends;
};
typedef struct regionpair * regionpairP;


struct path
{
	double ox;					//origin of the path
	double oy;
	double oz;
	double direcx;				//normalized direction vector pointing from ox,oy,oz to the neighbor, could save here
	double direcy;
	double direcz;

	//heterogeneity of the path
	regionpairP grainregions;

	//long* grainregions_outbound; //grain ids! dereference from defgpool the deformed grains that define the heterogeneities through which the grain boundaries migrate along the path
	//double* regionends_outbound; //defines where along the trajectory ox,oy,ox to the neighbor along the path a new region is starting or the current ends
	//long* grainregions_inbound; //inbound are utilized for the neighbor growing towards the nucleus to avoid
	//double* regionends_inbound; //inintuitive indexing strategies during the growth path calculation

	//path solution
	double distcontact; //distance along the path where the boundaries would meet solving the vector equation: o + distcontact * direc
	double timecontact; //when would the boundaries meet

	double totallength;
	int nsegments;

	char analyzedYet; //path already solved?
	//padding bytes could be utilized
};
typedef struct path * pathP;

struct vertex_l
{
	double x;
	double y;
	double z;
	double l;
};

struct contactinfo
{
	double when;
	double dist;
	double posx;
	double posy;
	double posz;
};
typedef struct contactinfo * contactinfoP;

struct box
{
	double xrd;
	double ynd;
	double ztd;
	double xyz;

	double tx;		//translation vector for each simulation box that defines how the local grid of nuclei is displaced relative to the most front left bottom defgr corner
	double ty;
	double tz;
	
	int ngrx;			//number of deformed grains in this box
	int ngry;
	int ngrz;
	int ngrxy;
	int ngrxyz;
};
typedef struct box * boxP;


/*
struct neighbor
{
	short nx;
	short ny;
	short nz;
	long represents;
};
typedef struct neighbor * neighborP;


struct cell
{
	bool activity;
	char rhoFac;
	short shape; //at some point get rid of this because it causes 6 padding bytes!
	short ix;
	short iy;
	short iz;
	unsigned short defgid;
	double rxFrac;
};
typedef struct cell * cellP;
*/

struct rxg
{
	long ori;					//an ori-index associated to that rxg

	double tincub;				//an incubation time in seconds after which the nucleus starts growing
	//double Vhull_prior;		//postprocessing 1
	//double Vhull_now;			//postprocessing 2
};
typedef struct rxg * rxgP;


struct rxgo
{
    double cx;
    double cy;
    double cz;
    double bunge1;
    double bunge2;
    double bunge3;
    long oid; //the orientation from the orientation pool
    long rxgpoolid;
    //double tincub;
};
typedef struct rxgo * rxgoP;


struct defg
{
	long ori;
	double rho0;
	double rho;
	double kuhlmann_tcrit;
};
typedef struct defg * defgP;


struct defglean
{
    long defgpoolid;
    long ori;
    double rho0;
};
typedef struct defglean * defgleanP;


struct ori
{
	double bunge1;		//Bunge (3,1,3) convention Euler angle
	double bunge2;
	double bunge3;

	double q0;			//unit quaternion representing the orientation
	double q1;
	double q2;
	double q3;

	long closestideal;
};
typedef struct ori * oriP;


struct ideal
{
	double bunge1;		//Bunge (3,1,3) convention Euler angles
	double bunge2;
	double bunge3;

	double q0;			//unit quaternion representing the orientation
	double q1;
	double q2;
	double q3;

	double scatter;
	
	long id;			//self-reflecting
};
typedef struct ideal * idealP;


struct optimtime
{
	double myPmax;
	double myPmin;
	double myRhomax;
};
typedef struct optimtime otpimtimeP;


struct grvolbinning
{
	double VolDistBinSize;
	double VolDistMaxValue;
	long N_Bins;
}; 
typedef struct grvolbinning grvolbinningP;


struct defmsinfo
{
    long nx;
    long ny;
    long nz;
    long nxy;
    long nxyz;
    double dimx; //grain size in RD (x), ND(y), TD(z)
    double dimy;
    double dimz;
    double gbnucdens2num;
    double gbnucdrho2dens;
    double gbnucmaxscatter; 
    
};
typedef struct defmsinfo defmsinfoP;

//MPI IO
typedef struct 
{
	double bunge1;		//Bunge (3,1,3) convention Euler angles
	double bunge2;
	double bunge3;

	double q0;			//unit quaternion representing the orientation
	double q1;
	double q2;
	double q3;

	double scatter;
	
	long id;
} MPI_IO_IdealOri;


typedef struct 
{
	long ori;
	double rho0;
	double rho;
	double kuhlmann_tcrit;
} MPI_IO_Defgrain;


typedef struct
{
	long ori;					//an ori-index associated to that rxg
	double tincub;				//an incubation time in seconds after which the nucleus starts growing
} MPI_IO_Rxgrain;


typedef struct
{
	double bunge1;		//Bunge (3,1,3) convention Euler angle
	double bunge2;
	double bunge3;

	double q0;			//unit quaternion representing the orientation
	double q1;
	double q2;
	double q3;

	long closestideal;
} MPI_IO_Ori;


typedef struct
{
	double t;
	double T;
} MPI_IO_Temperature;


typedef struct
{
	double t;
	double pz;
} MPI_IO_ZenerDrag;


typedef struct {
	long id;
	long period;
} MPI_IO_PrintSchedule;


#endif