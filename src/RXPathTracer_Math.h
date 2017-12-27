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


#ifndef RXPATHTRACER_MATH_H
#define RXPATHTRACER_MATH_H

#include <stdlib.h>
#include "RXPathTracer_Random.h"

//numerics
#define DOUBLE_ACCURACY 		5e-16				//2^-51 as determined with MATLAB eps(2)
#define MYMATH_STANDARD_SEED	-4256
#define MINIMUM_ANGULAR_SPREAD          ((2.0)*(0.017453292))		//1.0/180.0_PI_
#define MAXIMUM_MISORI_FCC		(1.09606677)		//FCC 62.8/180.0*_PI_
#define SYMMETRIES_IN_FCC		24
#define _PI_					3.1415926535897932384626433832795
#define SQR(a)					((a)*(a))

//IPFZ Coloring
#define FAIL_MYMATH_NUMERIC		(-1.0)
#define EPS_PROXIMITY_IPFZ		(0.01)
#define IPF_COLOR_STRETCH_R		(0.5)
#define IPF_COLOR_STRETCH_G		(0.5)
#define IPF_COLOR_STRETCH_B		(0.5)
#define EPS_ENVIRONMENT			(1e-7)
#define RGBRANGE				255


using namespace std;


class mathMethods
{
public:
	mathMethods( void );
	~mathMethods( void );

	//helper LB
	void bubbleSort ( double arr [ ], int size );
	void sortInt( long arr [], long size );
	void swap ( double& x, double& y );
	void swapInt ( long& x, long& y );
	void sort(int n, double *ra);
        
        double fac( long x );
        double poissondistribution( double lambda, long k );

	//quaternion algebra
	void multiplyQuaternions( double *q, double* p, double* r );
	void QuatOnVector3D( double *q, double* v, double* r );

	//convert orientations among parametrizations, such as Bunge (3,1,3) rotation convention and quaternion
	void euler2quaternion( double * euler, double * q );
	void quaternion2Euler( double * quat, double * euler );


	//calculate disorientation among two orientations in various parametrizations
	double misorientationCubic( double pa1, double Pa, double pa2, double pb1, double Pb, double pb2 );
	void misorientationQuaternionCubic( double* p, double* q, double* quat  );
	double misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 );


	//generate new orientations with some scatter
	void randomOrientationShoemake( double * result );
	void randomMisorientationShoemake( double theta, double* qr );
	
	void rotateOrientation( double *oriOri, double angle, double u, double v, double w, double *newOri );
        void newOrientationFromReference( double *oriOri, double deviation, double *newOri );
	void newOrientationFromReferenceQuat( double *qbunge, double deviation, double *newquat );

	//identify orientation by RGB scheme
	void devtorefEuler2RGB ( double *bunge, double *ideal, double maxDev, unsigned char *rgb); //blue channel stretch from 0.0 to maxDev in radian, all other orientations white
	void project2fundamentalregion_ipfz ( double *qtest, double *xy );
	void bunge2ipfz( double phi1, double PHI, double phi2, unsigned char *rgb, double * pos);

	//an own PRNG
//private:
	randomClass r;
};


#endif
