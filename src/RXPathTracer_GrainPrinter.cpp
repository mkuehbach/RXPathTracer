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
#include "RXPathTracer_GrainPrinter.h"
grainPrinter::grainPrinter(const string prefix_filename, double nucx, double nucy, double nucz, double &nuc_rads, long nFaces, face &nuc_faces, vertex &nuc_vertices, long nNeighbors, nucsite &nbors, double &nbors_radius, long &nbors_activated, const string dir){
	raw_prefix_FileName = prefix_filename;
	prefix_FileName = prefix_filename;
	nucleus_radii = &nuc_rads;
	numberOfFaces = nFaces;
	nucleus_faces = &nuc_faces;
	nucleus_vertices = &nuc_vertices;
	numberOfNeighbors = nNeighbors;
	neighbors_activated = &nbors_activated;
	neighbors = &nbors;
	neighbors_radius = &nbors_radius;
	neighbors_end_tstep.resize(numberOfNeighbors);
	for (long i = 0; i< numberOfNeighbors; i++){
	neighbors_end_tstep[i] = 0;
	}
	nucleus_pos[0] = nucx;
	nucleus_pos[1] = nucy;
	nucleus_pos[2] = nucz;
	directory_Name = dir;

	if( directory_Name != ""){
		prefix_FileName.insert(0, "./" + directory_Name + "/");
		mkdir(("./" + directory_Name + "/").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
}

grainPrinter::~grainPrinter(void){}

void grainPrinter::print( long tstep ){
	stringstream nuc_pov_fname, nuc_pov_fstream;
	stringstream nuc_vtp_fname, nuc_vtp_fstream;
	stringstream nbor_pov_fname, nbor_pov_fstream;

	nuc_pov_fname << prefix_FileName << ".GrainShape_" << tstep <<".pov";
	nuc_vtp_fname << prefix_FileName << ".GrainShape_" << tstep <<".vtp";
	nbor_pov_fname << prefix_FileName << ".NborShape_" << tstep <<".pov";
	/*	
	} else {
	nuc_pov_fname<< "./" << directory_Name << "/" << prefix_FileName << ".GrainShape_" << id <<".pov";
	nuc_vtp_fname<< "./" << directory_Name << "/" << prefix_FileName << ".GrainShape_" << id <<".vtp";
	nbor_pov_fname<< "./" << directory_Name << "/" << prefix_FileName << ".NborShape_" << id <<".pov";
	}
	*/
	
	print_pov_nborshape(nbor_pov_fstream);
	FILE* nbor_pov_file = fopen(nbor_pov_fname.str().c_str(),"w");
	fprintf(nbor_pov_file,nbor_pov_fstream.str().c_str());
	fclose(nbor_pov_file);
	nbor_pov_fstream.clear();

	print_pov_nucshape(nuc_pov_fstream);
	FILE* nuc_pov_file = fopen(nuc_pov_fname.str().c_str(),"w");
	fprintf(nuc_pov_file,nuc_pov_fstream.str().c_str());
	fclose(nuc_pov_file);
	nuc_pov_fstream.clear();

	nuc_vtp_fstream << scientific << setprecision(14);
	print_vtp_nucshape_manual(nuc_vtp_fstream);
	FILE* nuc_vtp_file = fopen(nuc_vtp_fname.str().c_str(),"w");
	fprintf(nuc_vtp_file,nuc_vtp_fstream.str().c_str());
	fclose(nuc_vtp_file);
	nuc_vtp_fstream.clear();

	print_vtp_nborshape(prefix_FileName, tstep);
}

void grainPrinter::print_nbor_time_evo( void ){
	for (long nb = 0; nb < numberOfNeighbors; nb ++){
		if (neighbors_end_tstep[nb] > 0){
			stringstream nbor_pvd_fname, nbor_pvd_fstream;
			nbor_pvd_fname << prefix_FileName << ".NborEvolution"<< nb <<".pvd";
			print_header_pvd(nbor_pvd_fstream);
			for (long tstep = 0; tstep < neighbors_end_tstep[nb]; tstep++){
				stringstream vtp_fname;
				vtp_fname << prefix_FileName << ".NborShape" << nb << "_" << tstep <<".vtp";
				if( file_exists(vtp_fname.str()) ){
					nbor_pvd_fstream << "<DataSet timestep=\"" << tstep <<"\" group=\"\" part=\"0\" " <<
						"file=\""<< raw_prefix_FileName << ".NborShape" << nb << "_" << tstep << ".vtp" << "\"/>" << endl;
				}
			}
			print_footer_pvd(nbor_pvd_fstream);
			FILE* nbor_pvd_file = fopen(nbor_pvd_fname.str().c_str(),"w");
			fprintf(nbor_pvd_file,nbor_pvd_fstream.str().c_str());
			fclose(nbor_pvd_file);
			nbor_pvd_fstream.clear();
		}
	}
}

inline void grainPrinter::print_header_pvd(stringstream &sstream){
	sstream << "<?xml version=\"1.0\"?>" <<endl;
	sstream << "<VTKFile type=\"Collection\" "<<
				"version=\"0.1\" " <<
				"byte_order=\"LittleEndian\" " <<
				"compressor=\"vtkZLibDataCompressor\">" <<endl;
	sstream << "<Collection>" << endl;
}

inline void grainPrinter::print_footer_pvd(stringstream &sstream){
	sstream << "</Collection>" << endl;
	sstream << "</VTKFile>";
}

void grainPrinter::print_vtp_nborshape(const string prefix_fname, long tstep ){
//for each neighbor that is activated
	for (long nb = 0; nb < numberOfNeighbors; nb ++){
		stringstream fname;
		fname << prefix_fname << ".NborShape" << nb << "_" << tstep <<".vtp";
		if (neighbors_activated[nb] == YES){
			if( tstep > neighbors_end_tstep[nb]) neighbors_end_tstep[nb] = tstep;
			vtkSmartPointer<vtkSphereSource> nborSphere =
				vtkSmartPointer<vtkSphereSource>::New();
			nborSphere->SetCenter(neighbors[nb].x - nucleus_pos[0], neighbors[nb].y - nucleus_pos[1], neighbors[nb].z - nucleus_pos[2]);
			nborSphere->SetRadius(neighbors_radius[nb]);
			nborSphere->SetPhiResolution(PHI_RES);
			nborSphere->SetThetaResolution(THETA_RES);
			nborSphere->Update();

			vtkSmartPointer<vtkXMLPolyDataWriter> writer =
				vtkSmartPointer<vtkXMLPolyDataWriter>::New();
			writer->SetFileName(fname.str().c_str());
			writer->SetInputConnection(nborSphere->GetOutputPort());
			//writer->SetDataMode(vtkXMLWriter::Ascii);
			writer->Write();
		}
	}
  //
}

void grainPrinter::print_vtp_nucshape_manual(stringstream &fstream){
	//4 3D-point-vectors
	long np;
	vertex p;
	double points [12];//x1,y1,z1, x2, ...
	long connectivity [12];
	points[0] = 0.;
	points[1] = 0.;
	points[2] = 0.;
	//print vtp-file-header
	print_header_vtk(fstream);
	for (long nf = 0; nf < numberOfFaces; nf++){
	//Save points A, B and C
		p =  v_multiply( nucleus_radii[nf], nucleus_vertices[nucleus_faces[nf].A] );
		points[3] = points[0] + p.x ; points[4] = points[1] + p.y ; points[5] = points[2] + p.z ;
		
		p =  v_multiply( nucleus_radii[nf], nucleus_vertices[nucleus_faces[nf].B] );
		points[6] = points[0] + p.x ; points[7] = points[1] + p.y ; points[8] = points[2] + p.z ;
		
		p =  v_multiply( nucleus_radii[nf], nucleus_vertices[nucleus_faces[nf].C] );
		points[9] = points[0] + p.x ; points[10] = points[1] + p.y ; points[11] = points[2] + p.z ;
	//Save connectivity
	//face 0 - A - B
	connectivity[0] = 0; connectivity[1] = 1; connectivity[2] = 2;
	//face 0 - A - C
	connectivity[3] = 0; connectivity[4] = 1; connectivity[5] = 3;
	//face 0 - B - C
	connectivity[6] = 0; connectivity[7] = 2; connectivity[8] = 3;
	//face A - B - C (pyramidal base-face)
	connectivity[9] = 1; connectivity[10] = 2; connectivity[11] = 3;

	//write the header for a tetragonal pyramid (consists of 4 points and 4 faces)
	print_piece_header_vtp(4, 0, 0, 0, 4, fstream);
	//print points
	print_points_vtp(&points[0], 4, fstream);
	//print faces
	print_triang_faces_vtp(&connectivity[0], 4, fstream);
	//write footer
	print_piece_footer_vtp(fstream);
	}
	print_footer_vtk(fstream);
}

inline void grainPrinter::print_header_vtk(stringstream &sstream, const string type){
	sstream << "<?xml version=\"1.0\"?>" << endl;
	sstream << "<VTKFile type=\"" << type << "\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
	sstream << "\t<" << type << ">" << endl;
}

inline void grainPrinter::print_footer_vtk(stringstream &sstream, const string type){
	sstream << "\t</"<< type << ">" << endl;
	sstream << "</VTKFile>" ;
}

inline void grainPrinter::print_piece_header_vtp(long numberOfPoints, long numberOfVerts, long numberOfLines, long numberOfStrips, long numberOfPolys, stringstream &sstream){
	sstream << "\t\t" << "<Piece NumberOfPoints=\"" << numberOfPoints << "\" NumberOfVerts=\"" << numberOfVerts <<"\" "
		<< "NumberOfLines=\"" << numberOfLines << "\" NumberOfStrips=\"" << numberOfStrips << "\" NumberOfPolys=\"" << numberOfPolys << "\">" << endl;
}

inline void grainPrinter::print_piece_footer_vtp( stringstream &sstream ){
	sstream << "\t\t" << "</Piece>" << endl;
}

inline void grainPrinter::print_points_vtp(double *points, long numberOfPoints, stringstream &sstream){
	sstream << "\t\t\t" << "<Points>" << endl;
	sstream << "\t\t\t\t" << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	sstream << "\t\t\t\t" ;
	for ( long p = 0; p < numberOfPoints; p++ ){
	 sstream << scientific << setprecision(12) << points[3*p] << " " << points[3*p+1] << " " << points[3*p+2] << " ";
	}
	sstream << endl;
	sstream << "\t\t\t\t" << "</DataArray>" << endl;
	sstream << "\t\t\t" << "</Points>" << endl;
}

inline void grainPrinter::print_triang_faces_vtp(long *connectivity,long numberOfTriangFaces, stringstream &sstream){
	sstream << "\t\t\t" << "<Polys>" <<endl;
	sstream << "\t\t\t\t" << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
	sstream << "\t\t\t\t" ;
	for ( long f = 0; f < numberOfTriangFaces ; f++ ){
	sstream << connectivity[3*f] << " " << connectivity[3*f+1] << " " << connectivity[3*f+2] << " ";
	}
	sstream << endl;
	sstream << "\t\t\t\t" << "</DataArray>" << endl;
	sstream << "\t\t\t\t" << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
	sstream << "\t\t\t\t" ;
	for ( long f = 0; f < numberOfTriangFaces ; f++ ){
	sstream << (f+1)*3 << " ";
	}
	sstream << endl;
	sstream << "\t\t\t\t" << "</DataArray>" << endl;
	sstream << "\t\t\t " << "</Polys>" << endl;
}
//


inline vertex grainPrinter::v_multiply(double factor, vertex &v){
	vertex v_out;
	v_out.x = factor * v.x;
	v_out.y = factor * v.y;
	v_out.z = factor * v.z;
	return v_out;
}

inline void grainPrinter::print_edge_pov(vertex v1, vertex v2, stringstream &sstream){
	sstream << "cylinder{<"<< v1.x << "," << v1.y << "," << v1.z << ">,<"
						<< v2.x << "," << v2.y << "," << v2.z << ">,r_ve}" << endl;
}

inline void grainPrinter::print_vertex_pov(vertex v, stringstream &sstream){
	sstream << "sphere{<" << v.x << "," << v.y<<","<< v.z << ">,r_ve}" << endl;
}

inline void grainPrinter::print_triangle_pov(vertex v1, vertex v2, vertex v3, stringstream &sstream){
	sstream << "triangle{<" << v1.x << "," << v1.y << "," << v1.z << ">, <"
						 << v2.x << "," << v2.y << "," << v2.z << ">, <"
						 << v3.x << "," << v3.y << "," << v3.z << ">}" << endl;
}

void grainPrinter::print_pov_nucshape(stringstream &fstream){
	for(long i_f = 0; i_f < numberOfFaces; i_f++){
		fstream<<scientific;
		if(nucleus_radii[i_f] > 1E-8){
			//cone-version
			/*fstream<<"cone{<"<<0<<","<<0<<","<<0<<">,"<<0.001<<",";
			fstream<<"<"<<1E6*nuc_rads[i_f]*ico_facenormals[i_f].x<<","<<1E6*nuc_rads[i_f]*ico_facenormals[i_f].y<<","<<1E6*nuc_rads[i_f]*ico_facenormals[i_f].z<<">,";
			fstream<<2E6*nuc_rads[i_f]/sqrt(numFaces)<<"}"<<endl;
			*/
			//triangular version
			//vertices
			print_vertex_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].A] ) , fstream );
			print_vertex_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].B] ) , fstream );
			print_vertex_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].C] ) , fstream );
			//edges vertice-vertice
			print_edge_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].A] ),
							v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].B]) , fstream );
			print_edge_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].A] ),
							v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].C]) , fstream );
			print_edge_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].B] ),
							v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].C]) , fstream );
			//edges to centroid
			vertex v_centroid;
			v_centroid.x = 0.;
			v_centroid.y = 0.;
			v_centroid.z = 0.;
			print_edge_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].A] ) , v_centroid, fstream);
			print_edge_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].B] ) , v_centroid, fstream);
			print_edge_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].C] ) , v_centroid, fstream);
			//4 pyramidal faces
			//base face
			print_triangle_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].A] ),
								v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].B] ),
								v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].C] ), fstream );
			//side faces
			print_triangle_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].A] ),
								v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].B] ), v_centroid, fstream);
			print_triangle_pov( v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].A] ),
								v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].C] ), v_centroid, fstream);
			print_triangle_pov(	v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].B] ),
								v_multiply( nucleus_radii[i_f], nucleus_vertices[nucleus_faces[i_f].C] ), v_centroid, fstream);
		}
	}
}

void grainPrinter::print_pov_nborshape(stringstream &fstream){
	fstream<<scientific;
	for (long i_nbor = 0; i_nbor < numberOfNeighbors; i_nbor++){
		if( neighbors_activated[i_nbor] == YES){
			fstream<<"sphere{<"<< ( neighbors[i_nbor].x - nucleus_pos[0] ) << ","
			   << ( neighbors[i_nbor].y - nucleus_pos[1] ) << ","
			   << ( neighbors[i_nbor].z - nucleus_pos[2] ) <<">,"
			   << neighbors_radius[i_nbor] <<"}"<<endl;
		}
	}
}

inline bool grainPrinter::file_exists(const string fname){
	struct stat buf;
	return (stat (fname.c_str(), &buf) == 0); 
}

//strucPrinter
strucPrinter::strucPrinter(const string pre_filename, double def_grain_x, double def_grain_y, double def_grain_z, boxP defboxes, double nucx, double nucy, double nucz, double &nuc_rads, long nFaces, face &nuc_faces, vertex &nuc_vertices, long nNeighbors, nucsite &nbors, double &nbors_radius, long &nbors_activated, const string dir)
	:grainPrinter(pre_filename, nucx, nucy, nucz, nuc_rads, nFaces,
	nuc_faces, nuc_vertices, nNeighbors, nbors, nbors_radius, nbors_activated, dir)
{
	appendOriData = false;
	appendDefIdData = false;
	defgr_x = def_grain_x;
	defgr_y = def_grain_y;
	defgr_z = def_grain_z;
	xWalls = new double[defboxes->ngrx + 1];
	yWalls = new double[defboxes->ngry + 1];
	zWalls = new double[defboxes->ngrz + 1];
	boxx = defboxes->xrd;
	boxy = defboxes->ynd;
	boxz = defboxes->ztd;
	defStructure = defboxes;
}

strucPrinter::~strucPrinter( void ){
	delete [] xWalls;
	delete [] yWalls;
	delete [] zWalls;
	if(appendOriData) delete [] oriIds;
	if(appendDefIdData) delete [] defgIds;
	
}
void strucPrinter::print_DefStructure( void ){
	make_walls();
	print_vtr_defStrucShape();
}

void strucPrinter::print_Box(){
	print_vtp_boxshape();
}

void strucPrinter::add_OriData(long* defgids, defgP defgpool){
	appendOriData = true;
	appendDefIdData = true;
	oriIds = new long[defStructure->ngrxyz];
	defgIds = new long[defStructure->ngrxyz];
	for (long ixyz = 0; ixyz < defStructure->ngrxyz; ixyz ++){
		oriIds[ixyz] = defgpool[defgids[ixyz]].ori;
		defgIds[ixyz] = defgids[ixyz];
	}

}

void strucPrinter::print_vtp_boxshape( void ){
	stringstream fname;
	fname << prefix_FileName << ".BoxShape.vtp";
	vtkSmartPointer<vtkCubeSource> box = vtkSmartPointer<vtkCubeSource>::New();
	box->SetCenter(0., 0. , 0.);
	box->SetXLength(boxx);
	box->SetYLength(boxy);
	box->SetZLength(boxz);
	box->Update();
	//
	vtkSmartPointer<vtkXMLPolyDataWriter> writer =
				vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(fname.str().c_str());
	writer->SetInputConnection(box->GetOutputPort());
	//writer->SetDataMode(vtkXMLWriter::Ascii);
	writer->Write();
}

void strucPrinter::print_vtr_defStrucShape( void ){
	stringstream fname;
	fname << prefix_FileName << ".DefStrucShape.vtr";
	
	vtkSmartPointer<vtkRectilinearGrid> defStrucGrid =
		vtkSmartPointer<vtkRectilinearGrid>::New();

	vtkSmartPointer<vtkDoubleArray> xCoords =
		vtkSmartPointer<vtkDoubleArray>::New();
	xCoords->SetNumberOfComponents(1);

	vtkSmartPointer<vtkDoubleArray> yCoords =
		vtkSmartPointer<vtkDoubleArray>::New();
	yCoords->SetNumberOfComponents(1);
  
	vtkSmartPointer<vtkDoubleArray> zCoords =
		vtkSmartPointer<vtkDoubleArray>::New();
	zCoords->SetNumberOfComponents(1);

	for (long ix = 0; ix <= defStructure->ngrx; ix ++){
		xCoords->InsertNextValue(xWalls[ix]);
	}

	for (long iy = 0; iy <= defStructure->ngry; iy ++){
		yCoords->InsertNextValue(yWalls[iy]);
	}

	for (long iz = 0; iz <= defStructure->ngrz; iz ++){
		zCoords->InsertNextValue(zWalls[iz]);
	}
	defStrucGrid->SetExtent(0, defStructure->ngrx, 0, defStructure->ngry, 0, defStructure->ngrz);
	defStrucGrid->SetXCoordinates(xCoords);
	defStrucGrid->SetYCoordinates(yCoords);
	defStrucGrid->SetZCoordinates(zCoords);
	cout << "nGrains= " << defStrucGrid->GetNumberOfCells() << endl;
	//OriData:
	if(appendOriData){
		vtkSmartPointer<vtkLongArray> oids = 
		  vtkSmartPointer<vtkLongArray>::New();
		oids->SetNumberOfComponents(1);
		oids->SetName("OrientationIDs");
		for (long ixyz = 0; ixyz < defStructure->ngrxyz; ixyz ++){
			oids->InsertNextValue(oriIds[ixyz]);
		}
		defStrucGrid->GetCellData()->AddArray(oids);
	}
	//
	if (appendDefIdData){
		vtkSmartPointer<vtkLongArray> defgids = 
		  vtkSmartPointer<vtkLongArray>::New();
		defgids->SetNumberOfComponents(1);
		defgids->SetName("DeformedGrainIDs");
		for (long ixyz = 0; ixyz < defStructure->ngrxyz; ixyz ++){
			defgids->InsertNextValue(defgIds[ixyz]);
		}
		defStrucGrid->GetCellData()->AddArray(defgids);
	}
	vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
    vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
	writer->SetFileName(fname.str().c_str());
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(defStrucGrid->GetProducerPort());
	#else
	writer->SetInputData(defStrucGrid);
	#endif
	//writer->SetDataMode(vtkXMLWriter::Ascii);
	writer->Write();
}

inline void strucPrinter::make_walls( void ){
	xWalls[0] = 0. - nucleus_pos[0]; 
	yWalls[0] = 0. - nucleus_pos[1];
	zWalls[0] = 0. - nucleus_pos[2];
	xWalls[defStructure->ngrx] = boxx - nucleus_pos[0];
	yWalls[defStructure->ngry] = boxy - nucleus_pos[1];
	zWalls[defStructure->ngrz] = boxz - nucleus_pos[2];
	for ( long ix = 1;  ix < defStructure->ngrx ; ix++){
		xWalls[ix] = ix * defgr_x + defStructure->tx - nucleus_pos[0];
	}
	for ( long iy = 1;  iy < defStructure->ngry ; iy++){
		yWalls[iy] = iy * defgr_y + defStructure->ty - nucleus_pos[1];
	}
	for ( long iz = 1;  iz < defStructure->ngrz ; iz++){
		zWalls[iz] = iz * defgr_z + defStructure->tz - nucleus_pos[2];
	}
}

//defgids[ixyz]... with ixyz = ix + (myBoxes[mn].ngrx * iy) + (myBoxes[mn].ngrxy * iz);