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

#include "RXPathTracer_Io.h"


void reportError( char * message) //int rank )
{
	//char * head = "ERROR: ";
	char err[BUFSIZ];
	sprintf(err,"%s%s","ERROR: ",message);
//	if( logFile && myRank == rank ) ioHdl.write( err, logFile );
	exitus( err );
}

void exitus (const char *s)
{
        printf("%s\n",s);
	//MPI_Finalize();
        exit(1);
}


io::io( void )
{

}

io::~io( void )
{

}

dataBlock::dataBlock( void )
{
	columnCount = 0;
	lineCount = 0;
	
	first = NULL;
	last = NULL;
}

dataBlock::~dataBlock( void )
{
	dataLineP dlp = NULL;
	dataLineP dlpnext = NULL;

	for( dlp = this->first; dlp; dlp = dlpnext )
	{
		dlpnext = dlp->next;
		delete dlp;
	}
}

dataLine::dataLine( void )
{
	dat = NULL;
	dataCount = 0;
	next = NULL;
	prev = NULL;
}

dataLine::dataLine( int size )
{
	dat = NULL;
	dataCount = size;
	next = NULL;
	prev = NULL;

	dat = (univDataP) calloc( dataCount, sizeof(univData) );
}

dataLine::~dataLine( void )
{
	for( int i=0; i<dataCount; i++ )
	{
		if( dat[i].type == 's' && dat[i].d.s )
			free(dat[i].d.s);
	}
	free( dat );
}

short io::open( short rw, const char * filename, FILE ** file )
{
	char rightStr[][3] = { "r", "w", "a" };
	char *rstr = rightStr[rw];

	(*file) = fopen(filename,rstr);
	if( (*file) ) return 1;
	return 0;
}

short io::open( const char * filename, ofstream * file )
{
	file->open(filename);
	if( !file->is_open() ) return 0;
	return 1;
}

short io::write( char * message, FILE * file )
{
	if( !file ) return FAILURE;
	return SUCCESS;
}

short io::write( string message, ofstream * file )
{
	cout << message;

	if( !file->is_open() )	return FAILURE;

	message.append("\n");

	(*file) << message;
	
	file->flush();

	return SUCCESS;
}

dataBlockP io::readDataBlock( const char *wanted, const char *fileName )
{
	FILE *fh;
	char buf1[BUFSIZE];
	char buf2[BUFSIZE];
	char dataTypes[100];
	char formatString[] = "%q";
	dataBlockP db = new dataBlock;

	fh = fopen( fileName, "r" );
	if(!fh)
	{
		char err[256];
		sprintf( err, "file %s does not exist\n", fileName ); 
		exitus(err); 
	}
	
	long blockExists = 0;
	while( fgets( buf1, BUFSIZE, fh ) )
	{
		char c1, c2;
		long cnt = sscanf( buf1, "%c%c%s%*s%s\n", &c1, &c2, db->name, dataTypes  );
		if( !cnt ) 		exitus("file appears to be empty!");

		if( cnt < 3 ) 		continue;
		if( c1 != '|' )		continue;
		if( c2 != '|' )		continue;

		if( stringsNotEqual( db->name, wanted ) ) 	continue;
		blockExists = 1;

		if( dataTypes[0] != '*' )	exitus("block header defective 1");
		if( dataTypes[1] != '(' )	exitus("block header defective 2");

		long i;
		db->columnCount = 0;
		for( i=0; dataTypes[i+2] != ')'; i++ )		db->columnCount++;
		break;
	}
	if( ! blockExists )
	{
		printf("Block -->%s<-- does not exist in file -->%s<--\n",wanted,fileName);
		delete db;
		return NULL;
	}

	while( fgets( buf1, BUFSIZE, fh ) )
	{	
		dataLineP line = new dataLine(db->columnCount);

		long weiter = 0;
		long i;
		for( i=0; i<db->columnCount; i++ )
		{
			formatString[1] = dataTypes[i+2];
			void * there = &line->dat[i].d;
			line->dat[i].type = dataTypes[i+2];

			while( (buf1[weiter]==' ') || (buf1[weiter]=='\t') ) 	weiter++;

			buf2[0] = 0;
			long cnt = sscanf( &buf1[weiter], "%s", buf2 );

			if( (!cnt) || ((buf2[0]=='|') && (buf2[1]=='|')) || (strlen(buf2)==0) )
			{
				if( i>0 )
				{
					char err[256];
					sprintf( err, "incomplete line in Block -->%s<--", wanted );
					exitus(err);
				}

				delete line;
				line = NULL;
				break;
			}

			weiter += strlen(buf2) +1;

			if(dataTypes[i+2] == 's')
			{
				there = malloc( strlen(buf2) +1 );
				line->dat[i].d.s = (char *) there;
			}

			cnt = sscanf( buf2, formatString, there );
			if( !cnt ) 	exitus("file error 6");
		}

		if( line )
		{
			db->lineCount++;
			if(db->last)	db->last->next = line;
			line->prev = db->last;
			line->next = NULL;
			db->last = line;
			if( !db->first )	db->first = line;
		}
		else
		{
			if( (buf2[0]=='|') && (buf2[1]=='|') )	break;
		}
	}
	fclose(fh);

	return db;
}

dataBlockP io::readDataBlock( const char *wanted, FILE *file )
{
	FILE *fh = file;
	char buf1[BUFSIZE];
	char buf2[BUFSIZE];
	char dataTypes[100];
	char formatString[] = "%q";
	dataBlockP db = new dataBlock;

	
	long blockExists = 0;
	while( fgets( buf1, BUFSIZE, fh ) )
	{
		char c1, c2;
		long cnt = sscanf( buf1, "%c%c%s%*s%s\n", &c1, &c2, db->name, dataTypes  );
		if( !cnt ) 		exitus("file appears to be empty!");

		if( cnt < 3 ) 	continue;
		if( c1 != '|' )		continue;
		if( c2 != '|' )		continue;

		if( stringsNotEqual( db->name, wanted ) ) 	continue;
		blockExists = 1;

		if( dataTypes[0] != '*' )	exitus("block header defective 1");
		if( dataTypes[1] != '(' )	exitus("block header defective 2");

		long i;
		db->columnCount = 0;
		for( i=0; dataTypes[i+2] != ')'; i++ )		db->columnCount++;
		break;
	}
	if( ! blockExists )
	{
		printf("Block -->%s<-- does not exist in file\n",wanted);
		delete db;
		return NULL;
	}

	while( fgets( buf1, BUFSIZE, fh ) )
	{

		dataLineP line = new dataLine(db->columnCount);

		long weiter = 0;
		long i;
		for( i=0; i<db->columnCount; i++ )
		{
			formatString[1] = dataTypes[i+2];
			void * there = &line->dat[i].d;
			line->dat[i].type = dataTypes[i+2];

			while( (buf1[weiter]==' ') || (buf1[weiter]=='\t') ) 	weiter++;

			buf2[0] = 0;
			long cnt = sscanf( &buf1[weiter], "%s", buf2 );

			if( (!cnt) || ((buf2[0]=='|') && (buf2[1]=='|')) || (strlen(buf2)==0) )
			{
				if( i>0 )
				{
					char err[256];
					sprintf( err, "incomplete line in Block -->%s<--", wanted );
					exitus(err);
				}

				delete line;
				line = NULL;
				break;
			}

			weiter += strlen(buf2) +1;

			if(dataTypes[i+2] == 's')
			{
				there = malloc( strlen(buf2) +1 );
				line->dat[i].d.s = (char *) there;
			}

			cnt = sscanf( buf2, formatString, there );
			if( !cnt ) 	exitus("file error 6");
		}

		if( line )
		{
			db->lineCount++;
			if(db->last)	db->last->next = line;
			line->prev = db->last;
			line->next = NULL;
			db->last = line;
			if( !db->first )	db->first = line;
		}
		else
		{
			if( (buf2[0]=='|') && (buf2[1]=='|') )	break;
		}
	}
	rewind(file);
    //fclose(fh); //###MK20130224, 20140302
	return db;
}



char *io::geTString( const char *s, dataBlockP db )
{
	dataLineP line;
	char err[256];

	for( line=db->first; line; line=line->next )
	{
                long last;
                char * parmValueString;

		char * parmNameString = line->dat[0].d.s;
		if( parmNameString[0] == '\"' )	parmNameString += 1;
		last = strlen(parmNameString) -1;
		if( parmNameString[last] == '\"' )	parmNameString[last] = 0;
		
		if( stringsNotEqual( s, parmNameString ) )		continue;

		parmValueString = line->dat[1].d.s;
		if( parmValueString[0] == '\"' )	parmValueString += 1;
		last = strlen(parmValueString) -1;
		if( parmValueString[last] == '\"' )	parmValueString[last] = 0;
		return parmValueString;
	}
	sprintf( err, "-->%s<-- gibts nicht im Block -->%s<--", s, db->name);
	exitus(err);
	return 0;
}


Real io::geTReal( const char *s, dataBlockP db ) //read first? value from a uds datablockP keyword
{
	char * parmValueString;
	double out;
	char dummy;
	long cnt;

    parmValueString = geTString( s, db );
	cnt = sscanf( parmValueString, "%lf%c", &out, &dummy );
	if( cnt < 1 )		exitus("can't read a Real.");
	if( cnt > 1 )		exitus("something's behind the Real");
	return out;
}


long io::geTInt( const char *s, dataBlockP db )
{
	char * parmValueString;
	long out;
	char dummy;
	long cnt;

	parmValueString = geTString( s, db );
	cnt = sscanf( parmValueString, "%ld%c", &out, &dummy );
	if( cnt < 1 )		exitus("can't read a long");
	if( cnt > 1 )		exitus("something's behind the long");
	return out;
}


Real io::getReal( dataLineP line, long column ) //read specific real from an uds matrix
{
 	QUICKASSERT( line->dat[column-1].type == 'f' );
	return line->dat[column-1].d.f;
}

Real io::getReal2( dataLineP line, long column )              //noch rictig machen!!, f¸r 2.Zahl aus Spalte
{
 	QUICKASSERT( line->dat[column-1].type == 'f' );
	return line->dat[column-1].d.f;
}


long io::getInt( dataLineP line, long column )
{
	QUICKASSERT( line->dat[column-1].type == 'i' );
	return line->dat[column-1].d.i;        //-1 for counting from 1... instead of 0...
}

char * io::getString( dataLineP line, long column )
{
	QUICKASSERT( line->dat[column-1].type == 's' );
	return line->dat[column-1].d.s;
}