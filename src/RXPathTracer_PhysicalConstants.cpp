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

#include "RXPathTracer_Defs.h"
#include "RXPathTracer_Types.h"
#include "RXPathTracer_PhysicalConstants.h"

/*
physicalConstants::physicalConstants(double inG0, double indGdT, double inTmelt, double inb0, double indbdT_C, double indbdT_a, double indbdT_b){
	G0 = inG0;
	dGdT = indGdT;
	Tmelt = inTmelt;
	b0 = inb0;
	dbdT_C = indbdT_C;
	dbdT_a = indbdT_a;
	dbdT_b = indbdT_b;
}

physicalConstants::physicalConstants(double * coefficients){
	G0 = coefficients[0];
	dGdT = coefficients[1];
	Tmelt = coefficients[2];
	b0 = coefficients[3];
	dbdT_C = coefficients[4];
	dbdT_a = coefficients[5];
	dbdT_b = coefficients[6];
}
*/


void physicalConstants::init(double inG0, double indGdT, double inTmelt, double inb0, double indbdT_C, double indbdT_a, double indbdT_b){
	G0 = inG0;
	dGdT = indGdT;
	Tmelt = inTmelt;
	b0 = inb0;
	dbdT_C = indbdT_C;
	dbdT_a = indbdT_a;
	dbdT_b = indbdT_b;
}


void physicalConstants::init(double * coefficients){
	G0 = coefficients[0];
	dGdT = coefficients[1];
	Tmelt = coefficients[2];
	b0 = coefficients[3];
	dbdT_C = coefficients[4];
	dbdT_a = coefficients[5];
	dbdT_b = coefficients[6];
}


void physicalConstants::get_Coefficients(double * coefficients){
	coefficients[0] = G0;
	coefficients[1] = dGdT;
	coefficients[2] = Tmelt;
	coefficients[3] = b0;
	coefficients[4] = dbdT_C;
	coefficients[5] = dbdT_a;
	coefficients[6] = dbdT_b;
}

void physicalConstants::update_Constants ( double T){
	update_G(T);
	update_b(T);
	update_halfG_b2();
}

inline void physicalConstants::update_G ( double T){
	//temperature dependent shear modulus, see Nadal, M-H and Le Poac, P in J. of Appl. Phys. 93 2472 2003 for further details what to do close to the melting point
	G = G0 - (dGdT * (T / Tmelt));
        
        
        //##MK::PhdMK
        //overwritten for iron in accordance with Ledstetter, Reed, 1973, Temperature is in K
	//MK::for bcc-iron
	G = G0 * (1.0 - ((9e-10 * CUBE(T) - 7e-7*SQR(T) + 0.0003*T - 0.00028)));
}

inline void physicalConstants::update_b ( double T){
	//temperature dependent Burgers vector as a result of lattice elongation/contraction
	double TinDegC = T - TOFFSET;
	b = b0 * (1. + (dbdT_C * (dbdT_a * TinDegC + dbdT_b * SQR(TinDegC) ) ) );
        
        
         //##MK::PhdMK
        //MK::for bcc iron up to 800deg celsius
        b = b0;
}

inline void physicalConstants::update_halfG_b2( void ){
	halfG_b2 = 0.5 * G * SQR(b);
}

double physicalConstants::get_G( void ){
	return G;
}

double physicalConstants::get_b( void ){
	return b;
}

double physicalConstants::get_halfG_b2( void ){
	return halfG_b2;
}

long physicalConstants::get_NumberOfCoefficients( void ){
	return 7;
}