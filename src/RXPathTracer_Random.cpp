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


#include "RXPathTracer_Random.h"


double randomClass::parkMiller( void ) 
{
	//MK::ok, seeds from -1 to -2^16 tested and acceptable Monobit and spectral pass rate
	//Random Number Generator Park-Miller-Bays-Durham "Numerical Recipes in C Second Edition (ran1), seed value zero is not allowed!
	//Ij+1 = ( a*Ij) % m with M = 2^31 - 1, a = 7^5 giving a period length of [0, m-1] with spectral property for triplets m^(1/k=3) = 1290 at most!
	//period is about 2.1e9
	int32_t j;
	int32_t k;
	static int32_t iy=0;
	static int32_t iv[NTAB];
	double temp;

	if( seed <= 0 || !iy )
	{
		if( -seed < 1 ) seed=1;
		else seed = -seed;
		for( j=NTAB+7;j>=0;j-- )
		{
			k = (int32_t) (seed/IQ);
			seed=(long) (IA*((int32_t)seed-k*IQ)-IR*k);
			if( seed < 0 ) seed += IM;
			if( j < NTAB ) iv[j] = (int32_t) seed;
		}
		iy=iv[0];
	}
	k=(int32_t)(seed/IQ);
	seed=(long)(IA*(seed-k*IQ)-IR*k);
	if( seed < 0 ) seed += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = (int32_t) seed;
	temp = AM*iy;

	if( (temp=AM*iy) > RNMX ) return RNMX;
	else return temp;
}


double randomClass::leEcuyer(void) 
{
	// L'Ecuyer combined as defined ran2 in Numerical Recipes in C Second Edition
	// combined period is at most 2.3e18
	//MK::ok, seeds can be negative
	int32_t j=0;
	int32_t k=0;
	static int32_t idum2 = 123456789;
	static int32_t iy = 0;
	static int32_t iv[NTABR2];
	double temp;

	if( seedr2 <= 0 )
	{
		if( -seedr2 < 1 ) seedr2 = 1;
		else seedr2 = -seedr2;

		idum2 = seedr2;

		for( j=NTAB+7; j>=0; j-- ) {
			k = (int32_t) (seedr2 / IQ1R2);
			seedr2 = (long) (IA1R2 * ( seedr2 - k * IQ1R2 ) - k * IR1R2);

			if( seedr2 < 0 ) seedr2 += IM1R2;
			if( j < NTABR2 ) iv[j] = (int32_t) seedr2;
		}
		iy=iv[0];
	}

	k = (int32_t) (seedr2/IQ1R2);
	seedr2 = (long) (IA1R2 * ( seedr2 - k * IQ1R2 ) - k * IR1R2);
	if( seedr2 < 0 ) seedr2 += IM1R2;
	k = (int32_t) (idum2 / IQ1R2);
	idum2 = IA2R2 * (idum2-k*IQ2R2)-k*IR2R2;
	if( idum2 < 0 ) idum2 += IM2R2;
	j = iy/NDIVR2;
	iy=iv[j]-idum2;
	iv[j] = (int32_t) seedr2;

	if( iy < 1 ) iy += IMM1R2;

	if( (temp = AMR2 * iy) > RNMXR2 ) return RNMXR2;
	else return temp;
}


long randomClass::randomIntervalLongInclusive_parkMiller( long lmin, long lmax ) 
{
	//MK::ok
	if( lmax < lmin ) return 0;
	lmax++; //making lmax inclusive

	double r = parkMiller();

	return (long) ( lmin + r * ( lmax - lmin ) );
}


long randomClass::randomIntervalLongInclusive_leEcuyer( long lmin, long lmax )
{
	//MK::ok
	if( lmax < lmin ) return 0;
	lmax++; //making lmax inclusive

	double r = leEcuyer();

	return lmin + ( long ) (r * ( lmax - lmin  ) );
}


double randomClass::randomInterval( double rmin, double rmax ) 
{
	//MK::ok
	if( rmax < rmin ) return 0.0;
	double r = parkMiller();
	return rmin + r * ( rmax - rmin );
}


void randomClass::generateShuffledLongArray( long lmin, long lmax, long * result ) 
{
	//MK::ok
	if( lmax <= lmin ) { result = NULL; return; }

	for( long i = 0; i < ( lmax - lmin ); i++ ) result[i] = i + lmin;

	for( long i = 0; i < ( lmax - lmin ); i++ )
	{
		//swop element at a random pos with the one at i successively
		long pos = i + randomIntervalLongInclusive_parkMiller( 0,  lmax - lmin - i - 1 );
		long temp = result[i];
		result[i] = result[pos];
		result[pos] = temp;
	}
}