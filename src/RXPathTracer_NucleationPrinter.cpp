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


#include <sys/stat.h>
#include "RXPathTracer_Defs.h"
#include "RXPathTracer_Types.h"
#include "RXPathTracer_NucleationPrinter.h"
nucleationPrinter::nucleationPrinter(const string pre_filename, boxP inBox, defgP inDefgpool, long *inDefgIds,
	double inDefgSize_x, double inDefgSize_y, double inDefgSize_z,
	nucsite *inNucleus,	nucsite *inNborNuclei, long inNumofNuclei, rxgP inRxgpool,
	oriP inOripool, const string dir){
	nucleiBox = inBox;
	defgpool = inDefgpool;
	defgx = inDefgSize_x;
	defgy = inDefgSize_y;
	defgz = inDefgSize_z;
	defgIds = inDefgIds;
	//
	nucleus = inNucleus;
	nbor_nuclei = inNborNuclei;
	numberOfNuclei = inNumofNuclei;
	rxgpool = inRxgpool;
	oripool = inOripool;
	xWalls = new double[nucleiBox->ngrx+1];
	yWalls = new double[nucleiBox->ngry+1];
	zWalls = new double[nucleiBox->ngrz+1];
	prefix_FileName = pre_filename;
	directory_Name = dir;

	if( directory_Name != ""){
		prefix_FileName.insert(0, "./" + directory_Name + "/");
		mkdir(("./" + directory_Name + "/").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
}
nucleationPrinter::~nucleationPrinter( void ){
	delete [] xWalls;
	delete [] yWalls;
	delete [] zWalls;
}
void nucleationPrinter::print( void ){
	
	stringstream fname, fstream;
	fname << prefix_FileName << ".Nucleation.uds";
	fstream << scientific << setprecision(14);
	print_box_header(fstream);
	print_box_data(fstream);
	print_deformation_struc_header(fstream);
	print_deformation_struc(fstream);
	print_nucleation_header(fstream);
	print_nuclei(fstream);


	FILE* f = fopen(fname.str().c_str(),"w");
	fprintf(f, fstream.str().c_str());
	fclose(f);
	fstream.clear();
}

void nucleationPrinter::print_box_data( stringstream &sstream){
	print_double_parameter(sstream, "BoxXmin", -.5 * nucleiBox->xrd);
	print_double_parameter(sstream, "BoxXmax", .5 * nucleiBox->xrd);
	
	print_double_parameter(sstream, "BoxYmin", -.5 * nucleiBox->ynd);
	print_double_parameter(sstream, "BoxYmax", .5 * nucleiBox->ynd);
	
	print_double_parameter(sstream, "BoxZmin", -.5 * nucleiBox->ztd);
	print_double_parameter(sstream, "BoxZmax", .5 * nucleiBox->ztd);
}
void nucleationPrinter::print_deformation_struc (stringstream &sstream){
	make_walls();
	double xmin, xmax;
	double ymin, ymax;
	double zmin, zmax;
	long ix, iy, iz;
	long ixyz = 0;
	for( iz = 0; iz < nucleiBox->ngrz; iz ++){
		zmin = zWalls[iz]; zmax = zWalls[iz+1];
		for (iy = 0; iy < nucleiBox->ngry; iy ++){
			ymin = yWalls[iy]; ymax = yWalls[iy+1];
			for (ix = 0; ix < nucleiBox->ngrx; ix ++){
				xmin = xWalls[ix]; xmax = xWalls[ix+1];
				ori orient = oripool[defgpool[defgIds[ixyz]].ori];
				double rho = defgpool[defgIds[ixyz]].rho0;
				print_defgrain(sstream, xmin, xmax, ymin, ymax, zmin, zmax,
					orient.bunge1, orient.bunge2, orient.bunge3, rho);
				ixyz ++;
			}
		}
	}
}

void nucleationPrinter::print_nuclei( stringstream &sstream){
	print_nucleus(sstream, *nucleus);
	for (long nbor = 0; nbor < numberOfNuclei; nbor ++){
		print_nucleus(sstream, nbor_nuclei[nbor]);
	}
}

inline void nucleationPrinter::print_box_header( stringstream &sstream ){
	sstream << "|| BoxInfo      || *(ss) || Parameter, Value" << endl;
}

inline void nucleationPrinter::print_double_parameter(stringstream &sstream, const string parName, double value){
	sstream << "\"" << parName << "\"" << "\t\t\t" << "\"" << value << "\"" << endl;
}

inline void nucleationPrinter::print_long_parameter(stringstream &sstream, const string parName, long value){
	sstream << "\"" << parName << "\"" << "\t\t\t" << "\"" << value << "\"" << endl;
}
inline void nucleationPrinter::print_deformation_struc_header (stringstream &sstream){
	sstream << "|| Deformation Structure || *(ffffffffff) || xmin, xmax, ymin, ymax, zmin, zmax, Phi1, PHI, Phi2, rho" << endl;
}

inline void nucleationPrinter::print_defgrain(stringstream &sstream, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
	double phi1, double Phi, double phi2, double rho){
	sstream << xmin << "\t" << xmax << "\t"  << ymin << "\t" << ymax << "\t" << zmin << "\t" << zmax << "\t";
	sstream << phi1 << "\t" << Phi << "\t" << phi2 << "\t" << rho << endl;
}


inline void nucleationPrinter::print_nucleation_header( stringstream &sstream ){
	sstream << "|| Nuclei      || *(ffffff) || x,y,z Phi1, PHI, Phi2" << endl;
}

inline void nucleationPrinter::print_nucleus(stringstream &sstream, nucsite &nuc){
	double relPos[3] = {nuc.x - .5 * nucleiBox->xrd,
							nuc.y - .5 * nucleiBox->ynd,
							nuc.z - .5 * nucleiBox->ztd};
	sstream << relPos[0] << "\t" << relPos[1] << "\t" << relPos[2] << "\t";
	sstream << oripool[rxgpool[nuc.rxgid].ori].bunge1 << "\t";
	sstream << oripool[rxgpool[nuc.rxgid].ori].bunge2 << "\t";
	sstream << oripool[rxgpool[nuc.rxgid].ori].bunge3 << endl;
}

inline void nucleationPrinter::make_walls( void ){
	xWalls[0] = 0. - nucleus->x;
	yWalls[0] = 0. - nucleus->y;
	zWalls[0] = 0. - nucleus->z;
	xWalls[nucleiBox->ngrx] = nucleiBox->xrd - nucleus->x;
	yWalls[nucleiBox->ngry] = nucleiBox->ynd - nucleus->y;
	zWalls[nucleiBox->ngrz] = nucleiBox->ztd - nucleus->z;
	for ( long ix = 1;  ix < nucleiBox->ngrx ; ix++){
		xWalls[ix] = ix * defgx + nucleiBox->tx - nucleus->x;
	}
	for ( long iy = 1;  iy < nucleiBox->ngry ; iy++){
		yWalls[iy] = iy * defgy + nucleiBox->ty - nucleus->y;
	}
	for ( long iz = 1;  iz < nucleiBox->ngrz ; iz++){
		zWalls[iz] = iz * defgz + nucleiBox->tz - nucleus->z;
	}
}