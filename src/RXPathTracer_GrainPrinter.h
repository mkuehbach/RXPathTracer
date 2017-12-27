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
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkCellData.h>
#include <vtkLongArray.h>
#include <vtkDoubleArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#define PHI_RES 50
#define THETA_RES 50

class grainPrinter {
public:
	grainPrinter(const string pre_filename,double nucx, double nucy, double nucz, double &nuc_rads,
			long nFaces, face &nuc_faces, vertex &nuc_vertices,
			long nNeighbors, nucsite &nbors, double &nbors_radius, long &nbors_activated,
			const string dir = "");
			~grainPrinter();
	void print( long tstep );
	void print_nbor_time_evo( void );

protected:
	string directory_Name;
	string prefix_FileName;
	string raw_prefix_FileName;
	//Nuc-Data
	double nucleus_pos[3];
	long numberOfFaces;
	double *nucleus_radii;
	face *nucleus_faces;
	vertex *nucleus_vertices;
	//Neighbor-Data
	long numberOfNeighbors;
	nucsite *neighbors;
	double *neighbors_radius;
	long *neighbors_activated;
	vector<long> neighbors_end_tstep;
	//
	void print_vtp_nucshape( void );
	void print_vtp_nborshape(const string prefix_fname, long id);
	void print_vtp_nucshape_manual(stringstream &fstream);
	
	//
	void print_pov_nborshape(stringstream &fstream);
	void print_pov_nucshape(stringstream &fstream);
	//Auxiliaries
	inline vertex v_multiply(double factor, vertex &v);
	//pvd-format
	inline bool file_exists(const string fname);
	inline void print_header_pvd(stringstream &sstream);
	inline void print_footer_pvd(stringstream &sstream);
	//vtk-format
	inline void print_header_vtk(stringstream &sstream, const string type="PolyData");
	inline void print_footer_vtk(stringstream &sstream, const string type="PolyData");
	inline void print_piece_header_vtp(long numberOfPoints, long numberOfVerts, long NumberOfLines, long NumberOfStrips, long NumberOfPolys, stringstream &sstream);
	inline void print_piece_footer_vtp( stringstream &sstream );
	inline void print_points_vtp(double *points, long numberOfPoints, stringstream &sstream);
	inline void print_triang_faces_vtp(long *connectivity,long numberOfTriangFaces, stringstream &sstream);
	//pov-format
	inline void print_vertex_pov(vertex v, stringstream &sstream);
	inline void print_edge_pov(vertex v1, vertex v2, stringstream &sstream);
	inline void print_triangle_pov(vertex v1, vertex v2, vertex v3, stringstream &sstream);
};

class strucPrinter: public grainPrinter{
	public: 
	strucPrinter (const string pre_filename, double def_grain_x, double def_grain_y, double def_grain_z, boxP defboxes,
		double nucx, double nucy, double nucz,
		double &nuc_rads, long nFaces, face &nuc_faces, vertex &nuc_vertices,
		long nNeighbors, nucsite &nbors, double &nbors_radius, long &nbors_activated,
		const string dir = "");
	~strucPrinter(void);
	
	void print_Box( void );
	void print_DefStructure( void );
	void add_OriData(long* defgids, defgP defgpool);
	private:
	bool appendOriData;
	bool appendDefIdData;
	long *oriIds;
	long *defgIds;
	double *xWalls;
	double *yWalls;
	double *zWalls;
	boxP defStructure;
	double boxx, boxy, boxz;
	double defgr_x, defgr_y, defgr_z;
	void print_vtp_boxshape( void );
	void print_vtr_defStrucShape( void );
	inline void make_walls(void);
};
