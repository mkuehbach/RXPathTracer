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

#ifndef RXPATHTRACER_IO_H
#define RXPATHTRACER_IO_H

#include "RXPathTracer_Defs.h"

typedef double Real;

#define BUFSIZE 1024


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
 


#define READ   0
#define WRITE  1
#define APPEND 2
#define FAILURE 0
#define SUCCESS 1


class dataLine;
class dataBlock; 
typedef dataLine * dataLineP;
typedef dataBlock * dataBlockP;


using namespace std;

void reportError( char * message); //int rank = 0);
void exitus (const char * s);


typedef struct
{
	union
	{
		float f;
		long i;
		char *s;
	} d;
	char type;

} univData;
typedef univData * univDataP;

class dataLine
{
public:
	dataLine( void );
	dataLine( int size );
	~dataLine( void );

	univDataP dat;
	long dataCount;

	dataLineP next;
	dataLineP prev;
};

class dataBlock
{
public:
	dataBlock( void );
	~dataBlock( void );
	long columnCount;
	long lineCount;
	
	dataLineP first;
	dataLineP last;

	char head[BUFSIZE];
	char name[128];
};

class io
{

public:
	io( void  );
	~io( void );

	short open( short rw, const char * filename, FILE ** file );
	short open( const char * filename, ofstream * file );
	short open( const char * filename, ifstream file );
	short write( char * message, FILE * file );
	short write( string message, ofstream * file );

	dataBlockP readDataBlock( const char *name, const char *fileName );
	dataBlockP readDataBlock( const char *name, FILE * file );

	char *newString(const char *s);
	Real getReal( dataLineP line, long column );
	Real getReal2( dataLineP line, long column );
	long getInt( dataLineP line, long column );
	char * getString( dataLineP line, long column );

	Real geTReal( const char *s, dataBlockP db );   
	long geTInt( const char *s, dataBlockP db );
	char * geTString( const char *s, dataBlockP db );

};


#endif	/* START1DCA_IO_H */

