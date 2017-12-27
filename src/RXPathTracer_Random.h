/*
	Copyright Markus Kühbach, Paul Wilhelm Hoffrogge, 2016

	RXPathTracer is an MPI-parallel program for simulating the growth stage during
	primary recrystallization. Its foundations are detailed in 
	
	  ''Efficient Recrystallization Microstructure Modeling by Utilizing Parallel Computation''
	  PhD thesis RWTH Aachen University
	  Successfully defended and accepted for publication August, 14th, 2017
	  Finally open source published in January, 2018
	  ´
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

#ifndef RXPATHTRACER_RANDOM_H
#define RXPATHTRACER_RANDOM_H

#include <math.h>
#include <time.h>
#include <stdint.h>
#include <float.h>

#define ISNAN(a) ((a) != (a))
#define _PI_ 3.1415926535897932384626433832795
#define EPSILON 0.01
#define EPS1 0.001
#define EPS2 1.0e-8
#define SQR(a) ((a)*(a))
#define CUBE(a) ((a)*(a)*(a))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define MYHASH(a,b) ( (((uint_fast64_t) max(a,b)) << 32) | ((uint_fast64_t) min(a,b)) )
#define SIGN(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)
#define XOR(a,b) ((!!a) != (!!b))


//defines constants relevant to the MLCG PRNG
//Park Miller Bays Durham parameter
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS DBL_EPSILON
#define RNMX (1.0-EPS)

//L'Ecuyer Combined MLCG ran2 Numerical Recipes C
#define IM1R2 2147483563
#define IM2R2 2147483399
#define AMR2 (1.0/IM1R2)
#define IMM1R2 (IM1R2-1)
#define IA1R2 40014
#define IA2R2 40692
#define IQ1R2 53668
#define IQ2R2 52774
#define IR1R2 12211
#define IR2R2 3791
#define NTABR2 32
#define NDIVR2 (1+IMM1R2/NTAB)
#define EPSR2 DBL_EPSILON
#define RNMXR2 (1.0-EPS)


#define DEFAULT_SEED	-46356

using namespace std;


class randomClass
{
public:
	randomClass( void ){ seed = DEFAULT_SEED; seedr2 = DEFAULT_SEED; }
	~randomClass( void ){};

	//setting the seeds
	void init( long sd ) { seed = sd; seedr2 = sd; }
	//void init(long sd) { long sd = (long) time(NULL); seed = sd; seedr2 = sd; }

	//type of implemented and tested generators
	double parkMiller( void );
	double leEcuyer( void );

	//interval functions
	double randomInterval( double rmin, double rmax );
	void generateShuffledLongArray( long lmin, long lmax, long * result );

	long randomIntervalLongInclusive_parkMiller( long lmin, long lmax );
	long randomIntervalLongInclusive_leEcuyer( long lmin, long lmax );


private:
	long seed;
	long seedr2;
};
typedef randomClass *randomClassP;


#endif