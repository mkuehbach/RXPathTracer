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


#ifndef RXPATHT_PHYS_CONST
#define RXPATHT_PHYS_CONST

class physicalConstants {
public:
	//physicalConstants(double inG0, double indGdT, double inTmelt, double inb0, double indbdT_C, double indbdT_a, double indbdT_b);
	//physicalConstants(double * coefficients);
	//physicalConstants(){};
	//~physicalConstants(){};

	void init( double inG0, double indGdT, double inTmelt, double inb0, double indbdT_C, double indbdT_a, double indbdT_b );
	void init( double * coefficients );

	void get_Coefficients( double * coefficients);
	long get_NumberOfCoefficients( void );
	void update_Constants( double T);
	double get_G( void );
	double get_b( void );
	double get_halfG_b2 ( void );
private:
	inline void update_G( double T);
	inline void update_b( double T);
	inline void update_halfG_b2 ( void );
		
	double G;//Youngs-Modulus
	//parameter for Youngs-Modulus dependence on temperature
	double G0, dGdT;
	double Tmelt;//melting temperature

	double b;//Burgers-Vector
	//parameter for Burgers-Vector dependence on temperature
	double b0;
	double dbdT_C;
	double dbdT_a;
	double dbdT_b;

	double halfG_b2; //0.5Gb^2
};
typedef class physicalConstants * physicalConstantsP;

#endif