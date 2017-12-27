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
 * File:   Definitions.h
 * Author: Markus
 *
 * Created on 27. Dezember 2013, 13:34
 */

#ifndef RXPATHTRACER_DEFS_H
#define	RXPATHTRACER_DEFS_H

#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>

#include <string>

#include <string.h>
#include <algorithm>

#include <mpi.h>
#include <stdexcept>

#include <limits.h>

//einbindung von stdlib.h und stdint.h is nicht notwendig, da bereits durch mymath eingebunden
//#include <stdlib.h>

#include "RXPathTracer_Io.h"
#include "RXPathTracer_Math.h"


using namespace std;

//FOR MPI
#define NDOUBLEPARAMETER 25
#define NLONGPARAMETER 18

#define MASTER					0



//icosphere definitions
#define ONETHIRD				(3.333333333333333333333333333e-1)
#define FOURTHIRDPI				(4.188790204786390984616857844372670512262892532500141094)
#define TSTEPMINALLOC			8
#define RADTODEG				(57.295779513082320876798154814105170)

//icosphere POVRAY definitions
#define N_MSPAINT 4

//nametags
#define SPHEREMESH_MODEL				1
#define VORONOI_MODEL					2


//possible path solutions
#define NOT_DETERMINED_YET              (0x00)
#define NUCLEUS_OVERGREW_NEIGHBOR       (0x01)
#define NEIGHBOR_OVERGREW_NUCLEUS       (0x02)
#define IMPINGEMENT_ALONG_GROWTHPATH    (0x03)

//memory guard
#define BYTE2MEGABYTE					1048576

//unit conversion
#define MICROMETER_PER_METER			(1.0e6)
#define MILLIMETER_PER_MICROMETER		(1.0e-3)


//memory management
#define	VECTOR_REALLOCATION_FACTOR		2
//
//--NOTE:
//ALWAYS CHANGE WHEN INTRODUCING NEW PARAMETERS!!!
//PARAMETERLIST SEE polyxx::rw_LongStorage/rw_DoubleStorage
//--

//define offsets
#define LEFT                    0
#define RIGHT                   1
#define FRONT                   2
#define REAR                    3 //mind different order than in COREV3
#define BOTTOM                  4
#define TOP                     5

#define SIXDIRECTIONS           6

#define DEFAULT_NORI            100
#define DEFAULT_NDEFG           100
#define DEFAULT_NRXG            100

#define NRXG                    10
#define NSEGMENTS               1000
#define SEGLENGTH               (1e-6)
#define RCRIT                   (0.5e-6) //critical radius of a nucleus
#define MAXFILLPERSTEP			(0.1)
#define RHOMAX_DEFAULT          (1e-6) //1/micron^2 gives very small initial velocity response

#define NO                      false
#define YES                     true

#define NOTCEASED_YET			false
#define UNCEASED				-1
#define CEASED                                  1

//parameter consistency checks and bounds
#define MINDIST_DEFAULT					(1e-3)	//nanometer mind micrometer
#define MINNEIGHBORS_DEFAULT			(13)
#define MAXBOXSIZE_DEFAULT				(1000.0)  //micron
#define MINIMUM_REDISCR_TIMEINTERVALS	100

//local neighbors
#define EPSILONBALL						(1e-3) //to avoid overlapping nucleation sites micron
#define BOUNDARYWIDTH_IN_BURGERSV		(2.0)
#define SITE_SATURATED					(0.0)
#define SMALLTIME						(1e-12) //ps
#define	MINIMUM_DEFGRAIN_MEDIAN			(1.0) //micron

//integrator accuracy limits to assure reasonable numerical solutions during transient annealing
#define INITIAL_DELTAT					(1e-6) //micro s
#define INFINITE						(1e36) //s
#define NO_INFECTION					(0.0) //recrystallized fraction
#define ALMOST_FULLY_INFECTED			(1.0)
#define NO_CELL_USED_YET				(0)

#define EPS_INPLANE_ACCURACY			(1e-3) //nm micron
#define EPS_SEGMENTS_ELIMDOUBLE			(1e-3) //section below a nm are not considered


#define FULLY_RECRYSTALLIZED			(65534) //defgid run from 0 to ndefgpool - 1
#define CURRENTLY_INFECTED				(65533)	//##MK::201407 was (-1)
#define NO_GRAINASSIGNED				(65532)

#define USHORT_GRAINS_MAX				(65530) //maximum number of grains possible
#define INACTIVE						false
#define MASTERNUC						0

#define FGEOFACE						(1.0000) //1.0000 //1.0449		//1.0449 // 1 //1
#define FGEOEDGE						(0.7957) //1.1301 //0.8024 //1.1347 // 0.7957 //0.707106781186547
#define FGEODIAG						(0.6494) //1.1000 //0.5821 //1.0082 // 0.6494 //0.577350269189626 //###MK20121211, 26NN


//discretization to time integration scheme
#define SMALL_HEATRATE			(0.000277777777777777) //1 K/ 3600s
#define SMALL_RECOVERYRATE		(2.77e-2) // 10^14 1/m^2/3600 now micron^-2 s^-1
#define SMALL_DRAGGINGRATE		(0.027777777e-6) // 100Pa/3600s in kg/(s^3um)
#define SMALL_HEAT				(1.0) //K
#define SMALL_RHO				(1.0e-01) //1/m^2 ~ 0.01 rho min micron^-2
#define SMALL_ZENERFORCE		(100.0e-6) //Pa
//#define SMALL_DISTANCE		this is the dcell
#define SMALL_VELOCITY			(2.77e-7) //1nm/3600s


//physical constants
#define kboltzman               (1.3806488e-23)
#define TOFFSET                 (273.15) //degree Celsius into Kelvin
#define KUHLMANNFINALRHO		(1e-2) //micron^-2
#define NEG_KBOLTZMANN			(-7.242971565e22) 
#define ECHARGE					(1.602176565e-19)


//mobility models
#define SEBALD_GOTTSTEIN        1
#define ROLLETT_HUMPHREYS       2

//disorientation/orientation related issues
#define DISORI                  ((1.5)*(0.01745329251994330)) // 1 degree raster
#define DISORICACHE_MEMLIMIT    (10.0) //maximum size of lower-triangle matrix implemented as implicitly addressed array of types double

//gb nucleation
#define LAGB_TO_HAGB_TRANSITION (0.261799387)
#define MINIMUM_DRHO			1.0e1
#define NO_NUCLEI				0
#define MINIMUM_DRHO			1.0e1
#define SCALING_LAMBDA			10.0
#define POISSON_CUMSUM_CUTOFF	19
#define POISSON_CUMSUM_TABLE	20 //cutoff+1
#define UNKNOWN_CANDIDATE       -1     
#define UNKNOWN                 -1

//JMAK correction global defines
#define FGEOSPHERE				(4.188790204786390)
#define FGEOSPHERESURF			(12.56637061435920)


//definition of ideal texture component / ideal / standardlagen
#define RANDOM					(-1)
#define MAX_FCC					(1.099) //MacKenzie..., 62.8/180*_PI_
#define PLACE_FOR_RANDOM		(1)


//definition triaxial ellipsoids
#define THREE_AXES			3
#define AXIS_X				0
#define AXIS_Y				1
#define AXIS_Z				2


//definitions that control I/O for RXGVolume
#define POS_BUNGE1			0
#define POS_BUNGE2			1
#define POS_BUNGE3			2
#define THREE_EULER_ANGLES	3

#define	NUC_VOLUME_LONG		0 //only fully transformed cells
#define NUC_VOLUME_DOUBLE	1 //also partial infection
#define NUC_SURFACE_LONG	2 //discrete cells currently active in the discrete interfaces


//passed arguments
#define INPUTFILE				1
#define	ID						2
#define MODELTYPE				3
/*
void exitus (const char *s);
#define _STR(x) _VAL(x)
#define _VAL(x) #x
#define ASSERT(cond) if(!(cond)) exitus("failed assertion:"__FILE__"line"_STR(__LINE__)":"#cond)
#define TOLFP(x) ((1.0)-(0.99999/((double) x))) //tolerance
#define ERRTXT(text) (text" : file: "__FILE__" line:"_STR(__LINE__)"\n")
#define stringsNotEqual(a,b) strcmp(a,b)
#define stringsEqual(a,b) !stringsNotEqual(a,b)


#define DEBUG

#ifndef DEBUG
        #define QUICKASSERT(cond)       ((void)0)
#else
        #define QUICKASSERT(cond)       ASSERT(cond)
#endif
*/


//class randomClass;
//typedef randomClass * randomClassP;

/*
class mathMethods;
typedef class mathMethods * mathMethodsP;
 */


#endif	/* DEFINITIONS_H */
