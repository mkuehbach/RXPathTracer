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


class nucleationPrinter {
public:
	nucleationPrinter(const string pre_filename, boxP inBox, defgP inDefgpool, long *inDefgIds, double inDefgSize_x, double inDefgSize_y, double inDefgSize_z, nucsite *inNucleus,
		nucsite *inNborNuclei, long inNumofNuclei, rxgP inRxgpool, oriP inOripool, const string dir = "");
	~nucleationPrinter( void );
	void print( void );
private:
	void print_box_data(stringstream &sstream);
	void print_deformation_struc (stringstream &sstream);
	void print_nuclei(stringstream &sstream);
	inline void print_box_header(stringstream &sstream);
	inline void print_long_parameter(stringstream &sstream, const string parName, long value);
	inline void print_double_parameter(stringstream &sstream, const string parName, double value);
	inline void print_deformation_struc_header (stringstream &sstream);
	inline void print_defgrain(stringstream &sstream, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
					double phi1, double Phi, double phi2, double rho);
	inline void print_nucleation_header(stringstream &sstream);
	inline void print_nucleus(stringstream &sstream, nucsite &nuc);
	inline void make_walls ( void );
	string directory_Name;
	string prefix_FileName;

	boxP nucleiBox;
	defgP defgpool;
	long *defgIds;
	double *xWalls, *yWalls, *zWalls;
	double defgx, defgy, defgz;
	nucsite *nucleus;
	nucsite *nbor_nuclei;
	long numberOfNuclei;
	rxgP rxgpool;
	oriP oripool;
};

