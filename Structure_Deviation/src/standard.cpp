#include "trigonalplanar.h"
#include "tetragonalplanar.h"
#include "tetrahedral.h"
#include "trigonal_bipyramidal.h"
#include "tetragonal_pyramidal.h"
#include "octahedral.h"

	void triplMetalicStructure::form_standard_triangle_NA(){
		_metal_atom=atom("HETATM", 895, "NA", " ", "NA", "A", 1339, 5.059, 13.759, 13.331, 1.00, 28.94, 0.0, "NA", 0.0);
		_tripl_atoms[0]=atom("ATOM", 382, "O", " ", "TYR", "A", 1156, 4.137, 11.329, 14.623, 1.00, 17.64, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 394, "O", " ", "ILE", "A", 1157, 6.825, 12.291, 12.638, 1.00, 14.52, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("HETATM", 953, "O", " ", "HOH", "A", 1458, 3.773, 13.122, 11.342, 1.00, 25.84, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_CA(){
		_metal_atom=atom("HETATM", 6750, "CA", " ", "CA", "B", 2057, 81.101, 13.281, 8.855, 1.00, 67.95, 0.0, "CA", 0.0);
		_tripl_atoms[0]=atom("ATOM", 1488, "OD1", " ", "ASP", "B", 57, 82.105, 10.333, 9.650, 1.00, 92.35, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 1496, "OD1", " ", "ASP", "B", 58, 79.223, 11.188, 8.490, 1.00, 96.01, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("ATOM", 2432, "NH1", " ", "ARG", "C", 48, 79.196, 15.145, 7.164, 1.00, 90.39, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_K(){
		_metal_atom=atom("HETATM", 1042, "K", " ", "K", "A", 204, 22.118, -14.492, 5.832, 0.80, 17.63, 0.0, "K", 0.0);
		_tripl_atoms[0]=atom("ATOM", 179, "O", "A", "GLN", "A", 645, 22.664, -16.060, 3.601, 0.50, 16.67, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 221, "O", " ", "ASP", "A", 648, 21.281, -12.947, 3.830, 1.00, 15.58, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("ATOM", 243, "OD1", " ", "ASN", "A", 651, 21.439, -12.632, 7.461, 1.00, 22.46, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_MG(){
		_metal_atom=atom("HETATM", 5770, "MG", " ", "MG", "C", 302, -27.090, 30.865, -49.081, 1.00, 50.19, 0.0, "MG", 0.0);
		_tripl_atoms[0]=atom("ATOM", 3506, "O", " ", "ALA", "C", 96, -28.842, 29.139, -47.430, 1.00, 32.85, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 3657, "OE2", " ", "GLU", "C", 116, -29.574, 31.826, -47.961, 1.00, 35.19, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("HETATM", 5934, "O", " ", "HOH", "C", 432, -29.352, 31.651, -50.747, 1.00, 35.30, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_MN(){
		_metal_atom=atom("HETATM", 5059, "MN", " ", "MN", "C", 500, 116.728, -3.106, 35.047, 1.00, 9.05, 0.0, "MN", 0.0);
		_tripl_atoms[0]=atom("ATOM", 779, "OD2", " ", "ASP", "C", 271, 118.598, -2.316, 33.962, 1.00, 6.13, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 1044, "OD1", " ", "ASN", "C", 303, 117.413, -4.988, 34.938, 1.00, 8.69, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("HETATM", 5077, "O3", "A", "NHC", "C", 0, 115.855, -1.719, 33.655, 0.50, 8.82, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_FE(){
		_metal_atom=atom("HETATM", 2271, "FE", " ", "FE", "A", 1351, 9.707, 12.322, -2.616, 1.00, 60.23, 0.0, "FE", 0.0);
		_tripl_atoms[0]=atom("ATOM", 837, "OE1", " ", "GLU", "A", 169, 8.583, 10.420, -1.934, 1.00, 54.16, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 1578, "OE2", " ", "GLU", "A", 266, 11.340, 13.345, -1.879, 1.00, 51.09, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("ATOM", 1597, "ND1", " ", "HIS", "A", 269, 10.934, 10.824, -4.046, 1.00, 43.95, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_CU(){
		_metal_atom=atom("HETATM", 4776, "CU", "A", "CU", "F", 154, -0.616, 15.483, 9.471, 0.55, 5.31, 0.0, "CU", 0.0);
		_tripl_atoms[0]=atom("ATOM", 3090, "ND1", " ", "HIS", "F", 46, 0.690, 16.307, 8.129, 1.00, 7.90, 0.0, "N", 0.0);
		_tripl_atoms[1]=atom("ATOM", 3127, "NE2", " ", "HIS", "F", 48, -0.575, 14.577, 11.169, 1.00, 6.87, 0.0, "N", 0.0);
		_tripl_atoms[2]=atom("ATOM", 4244, "NE2", " ", "HIS", "F", 120, -2.410, 16.095, 8.723, 1.00, 7.10, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_ZN(){
		_metal_atom=atom("HETATM", 2332, "ZN", "A", "ZN", "B", 3, 1.053, 35.245, -0.862, 0.50, 8.65, 0.0, "ZN", 0.0);
		_tripl_atoms[0]=atom("ATOM", 1815, "OE2", " ", "GLU", "B", 113, 2.594, 35.155, -2.045, 1.00, 7.97, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 2224, "OE2", " ", "GLU", "B", 160, 2.147, 36.210, 0.789, 1.00, 7.53, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("HETATM", 2718, "O", " ", "HOH", "B", 408, 0.034, 33.525, 0.156, 0.50, 3.52, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_CO(){
		_metal_atom=atom("HETATM", 1420, "CO", " ", "CO", "A", 204, 0.000, -47.565, 0.000, 0.25, 17.90, 0.0, "CO", 0.0);
		_tripl_atoms[0]=atom("ATOM", 1396, "NE2", " ", "HIS", "A", 173, 2.124, -47.567, -0.796, 1.00, 11.16, 0.0, "N", 0.0);
		_tripl_atoms[1]=atom("HETATM", 1422, "CL", " ", "CL", "A", 206, 0.000, -45.061, 0.000, 0.25, 13.25, 0.0, "CL", 0.0);
		_tripl_atoms[2]=atom("HETATM", 1573, "O", " ", "HOH", "A", 448, 0.000, -49.572, 0.000, 0.25, 19.08, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_PB(){
		_metal_atom=atom("HETATM", 5845, "PB", "A", "PB", "A", 1101, 7.203, 17.128, 19.342, 0.25, 51.08, 0.0, "PB", 0.0);
		_tripl_atoms[0]=atom("HETATM", 5913, "NA", " ", "PP9", "A", 901, 5.417, 17.028, 17.803, 1.00, 49.79, 0.0, "N", 0.0);
		_tripl_atoms[1]=atom("HETATM", 5940, "ND", " ", "PP9", "A", 901, 5.732, 16.116, 20.863, 1.00, 49.00, 0.0, "N", 0.0);
		_tripl_atoms[2]=atom("HETATM", 5963, "OXT", " ", "ACY", "A", 803, 8.331, 16.515, 18.965, 0.75, 65.20, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_HG(){
		_metal_atom=atom("HETATM", 2041, "HG", " ", "HG", "A", 269, 36.267, 14.504, -1.501, 0.10, 13.42, 0.0, "HG", 0.0);
		_tripl_atoms[0]=atom("ATOM", 1141, "N", " ", "LEU", "A", 147, 38.242, 14.978, -3.627, 1.00, 9.24, 0.0, "N", 0.0);
		_tripl_atoms[1]=atom("ATOM", 1658, "O", " ", "CYS", "A", 212, 33.729, 15.292, -0.643, 1.00, 10.29, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("ATOM", 1673, "O", " ", "GLU", "A", 214, 36.480, 12.601, 0.363, 1.00, 11.08, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_AS(){
		_metal_atom=atom("HETATM", 3960, "AS", " ", "CAC", "A", 7, -9.304, 54.486, 22.055, 0.70, 13.57, 0.0, "AS", 0.0);
		_tripl_atoms[0]=atom("ATOM", 328, "SG", " ", "CYS", "A", 68, -7.590, 54.797, 23.518, 1.00, 10.93, 0.0, "S", 0.0);
		_tripl_atoms[1]=atom("HETATM", 3961, "C1", " ", "CAC", "A", 7, -9.675, 56.399, 21.846, 0.70, 11.32, 0.0, "C", 0.0);
		_tripl_atoms[2]=atom("HETATM", 3962, "C2", " ", "CAC", "A", 7, -8.234, 53.882, 20.524, 0.70, 14.80, 0.0, "C", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_CD(){
		_metal_atom=atom("HETATM", 1439, "CD", " ", "CD", "A", 201, 50.949, 37.868, -24.756, 0.50, 20.38, 0.0, "CD", 0.0);
		_tripl_atoms[0]=atom("ATOM", 89, "OD1", " ", "ASP", "A", 15, 49.607, 37.511, -22.928, 1.00, 19.83, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 90, "OD2", " ", "ASP", "A", 15, 48.707, 37.459, -24.900, 1.00, 28.73, 0.0, "O", 0.0);
		_tripl_atoms[2]=atom("HETATM", 1608, "O", " ", "HOH", "A", 445, 52.844, 37.855, -22.866, 0.50, 29.98, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	void triplMetalicStructure::form_standard_triangle_NI(){
		_metal_atom=atom("HETATM", 1215, "NI", " ", "NI", "A", 101, 5.239, -8.573, 6.275, 1.00, 10.61, 0.0, "NI", 0.0);
		_tripl_atoms[0]=atom("ATOM", 40, "OE1", " ", "GLU", "A", 8, 3.512, -8.069, 5.448, 1.00, 19.05, 0.0, "O", 0.0);
		_tripl_atoms[1]=atom("ATOM", 17, "ND1", " ", "HIS", "A", 5, 6.556, -9.198, 4.803, 1.00, 13.01, 0.0, "N", 0.0);
		_tripl_atoms[2]=atom("ATOM", 390, "NE2", " ", "HIS", "A", 47, 6.082, -6.787, 6.872, 1.00, 9.14, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tripl_vertices[i]=Point3D(_tripl_atoms[i].X(), _tripl_atoms[i].Y(), _tripl_atoms[i].Z());
		}

		_tripl=triangle(_tripl_vertices);
		_tripl.find_vertex_order();
		print(_metal_atom, _tripl_atoms);

		cout<<"\nVertices:";
		print(_tripl_vertices);
	}

	string triplMetalicStructure::assign_standard_site_type(){
		string site_type;
		int aa=0, dna=0, rna=0, other=0, water=0;
		for (int x=0; x<_no_of_binding_atoms; x++){
			if (_tripl_atoms[x].get_molecule_type()=="amino_acid"){
				aa++;
			}else if (_tripl_atoms[x].get_molecule_type()=="dna_molecule"){
				dna++;
			}else if (_tripl_atoms[x].get_molecule_type()=="rna_molecule"){
				rna++;
			}else if (_tripl_atoms[x].get_molecule_type()=="other_biom"){
				other++;
			}else if (_tripl_atoms[x].get_residue_name()=="HOH"){
				water++;
			}
		}

		if (water==_no_of_binding_atoms){
			site_type="all_water";
		}else if ((aa!=0) && (dna==0) && (rna==0)){
			site_type="protein";
		}else if ((aa==0) && (dna!=0) && (rna==0)){
			site_type="dna";
		}else if ((aa==0) && (dna==0) && (rna!=0)){
			site_type="rna";
		}else if ((aa!=0) && (dna!=0) && (rna==0)){
			site_type="protein_dna_complex";
		}else if ((aa!=0) && (dna==0) && (rna!=0)){
			site_type="protein_rna_complex";
		}else if ((aa==0) && (dna!=0) && (rna!=0)){
			site_type="dna_rna_complex";
		}else if ((aa!=0) && (dna!=0) && (rna!=0)){
			site_type="protein_dna_rna_complex";
		}else{
			site_type="other";
		}

		return site_type;
	}

	void tetraplMetalicStructure::form_standard_tetragonal_planar_MG(){
		_metal_atom=atom("HETATM", 3775, "MG", " ", "MG", "B", 2001, 37.085, -10.544, 3.367, 0.46, 18.33, 0.0, "MG", 0.0);
		_tetrapl_atoms[0]=atom("ATOM", 2009, "OD1", " ", "ASP", "B", 1781, 38.410, -9.690, 4.611, 1.00, 18.02, 0.0, "O", 0.0);
		_tetrapl_atoms[1]=atom("HETATM", 4299, "O", " ", "HOH", "B", 2213, 36.061, -11.332, 2.014, 0.42, 23.52, 0.0, "O", 0.0);
		_tetrapl_atoms[2]=atom("HETATM", 4300, "O", " ", "HOH", "B", 2214, 35.332, -10.119, 4.343, 0.46, 16.81, 0.0, "O", 0.0);
		_tetrapl_atoms[3]=atom("HETATM", 4410, "O", " ", "HOH", "B", 2321, 38.567, -10.761, 2.190, 0.47, 23.06, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetrapl_vertices[i]=Point3D(_tetrapl_atoms[i].X(), _tetrapl_atoms[i].Y(), _tetrapl_atoms[i].Z());
		}

		_tetrapl=tetragonalplanar(_tetrapl_vertices);
		print(_metal_atom, _tetrapl_atoms);

		cout<<"\nVertices:";
		print(_tetrapl_vertices);
	}

	void tetraplMetalicStructure::form_standard_tetragonal_planar_FE(){
		_metal_atom=atom("HETATM", 5406, "FE", " ", "HEM", "C", 1142, 22.296, 5.242, -2.842, 1.00, 8.33, 0.0, "FE", 0.0);
		_tetrapl_atoms[0]=atom("HETATM", 5402, "NA", " ", "HEM", "C", 1142, 23.471, 6.317, -4.181, 1.00, 9.25, 0.0, "N", 0.0);
		_tetrapl_atoms[1]=atom("HETATM", 5403, "NB", " ", "HEM", "C", 1142, 23.428, 6.022, -1.301, 1.00, 8.53, 0.0, "N", 0.0);
		_tetrapl_atoms[2]=atom("HETATM", 5404, "NC", " ", "HEM", "C", 1142, 21.803, 3.653, -1.575, 1.00, 7.90, 0.0, "N", 0.0);
		_tetrapl_atoms[3]=atom("HETATM", 5405, "ND", " ", "HEM", "C", 1142, 21.830, 3.997, -4.453, 1.00, 8.75, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetrapl_vertices[i]=Point3D(_tetrapl_atoms[i].X(), _tetrapl_atoms[i].Y(), _tetrapl_atoms[i].Z());
		}

		_tetrapl=tetragonalplanar(_tetrapl_vertices);
		print(_metal_atom, _tetrapl_atoms);

		cout<<"\nVertices:";
		print(_tetrapl_vertices);
	}

	void tetraplMetalicStructure::form_standard_tetragonal_planar_CD(){
		_metal_atom=atom("HETATM", 804, "CD", " ", "CD", "A", 262, -0.392, 3.557, 0.089, 0.50, 18.27, 0.0, "CD", 0.0);
		_tetrapl_atoms[0]=atom("ATOM", 781, "OE1", " ", "GLU", "A", 258, -0.076, 1.474, -1.738, 1.00, 17.63, 0.0, "O", 0.0);
		_tetrapl_atoms[1]=atom("ATOM", 782, "OE2", " ", "GLU", "A", 258, -1.718, 2.353, -0.673, 1.00, 20.31, 0.0, "O", 0.0);
		_tetrapl_atoms[2]=atom("HETATM", 816, "O", " ", "HOH", "A", 273, 0.402, 5.192, 1.755, 1.00, 30.05, 0.0, "O", 0.0);
		_tetrapl_atoms[3]=atom("HETATM", 901, "O", " ", "HOH", "A", 358, -2.277, 4.524, 1.514, 1.00, 40.63, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetrapl_vertices[i]=Point3D(_tetrapl_atoms[i].X(), _tetrapl_atoms[i].Y(), _tetrapl_atoms[i].Z());
		}

		_tetrapl=tetragonalplanar(_tetrapl_vertices);
		print(_metal_atom, _tetrapl_atoms);

		cout<<"\nVertices:";
		print(_tetrapl_vertices);
	}

	void tetraplMetalicStructure::form_standard_tetragonal_planar_NI(){
		_metal_atom=atom("HETATM", 667, "NI", " ", "NI", "A", 398, 13.849, 21.950, 30.990, 1.00, 21.68, 0.0, "NI", 0.0);
		_tetrapl_atoms[0]=atom("ATOM", 1, "N", " ", "GLY", "A", 288, 12.334, 23.080, 30.892, 1.00, 20.93, 0.0, "N", 0.0);
		_tetrapl_atoms[1]=atom("ATOM", 5, "N", " ", "SER", "A", 289, 12.992, 21.082, 29.673, 1.00, 28.62, 0.0, "N", 0.0);
		_tetrapl_atoms[2]=atom("ATOM", 11, "N", " ", "HIS", "A", 290, 15.185, 20.681, 30.963, 1.00, 22.82, 0.0, "N", 0.0);
		_tetrapl_atoms[3]=atom("ATOM", 17, "ND1", " ", "HIS", "A", 290, 14.740, 23.013, 32.337, 1.00, 19.58, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetrapl_vertices[i]=Point3D(_tetrapl_atoms[i].X(), _tetrapl_atoms[i].Y(), _tetrapl_atoms[i].Z());
		}

		_tetrapl=tetragonalplanar(_tetrapl_vertices);
		print(_metal_atom, _tetrapl_atoms);

		cout<<"\nVertices:";
		print(_tetrapl_vertices);
	}

	string tetraplMetalicStructure::assign_standard_site_type(){
		string site_type;
		int aa=0, dna=0, rna=0, other=0, water=0;
		for (int x=0; x<_no_of_binding_atoms; x++){
			if (_tetrapl_atoms[x].get_molecule_type()=="amino_acid"){
				aa++;
			}else if (_tetrapl_atoms[x].get_molecule_type()=="dna_molecule"){
				dna++;
			}else if (_tetrapl_atoms[x].get_molecule_type()=="rna_molecule"){
				rna++;
			}else if (_tetrapl_atoms[x].get_molecule_type()=="other_biom"){
				other++;
			}else if (_tetrapl_atoms[x].get_residue_name()=="HOH"){
				water++;
			}
		}

		if (water==_no_of_binding_atoms){
			site_type="all_water";
		}else if ((aa!=0) && (dna==0) && (rna==0)){
			site_type="protein";
		}else if ((aa==0) && (dna!=0) && (rna==0)){
			site_type="dna";
		}else if ((aa==0) && (dna==0) && (rna!=0)){
			site_type="rna";
		}else if ((aa!=0) && (dna!=0) && (rna==0)){
			site_type="protein_dna_complex";
		}else if ((aa!=0) && (dna==0) && (rna!=0)){
			site_type="protein_rna_complex";
		}else if ((aa==0) && (dna!=0) && (rna!=0)){
			site_type="dna_rna_complex";
		}else if ((aa!=0) && (dna!=0) && (rna!=0)){
			site_type="protein_dna_rna_complex";
		}else{
			site_type="other";
		}

		return site_type;
	}

	void tetraMetalicStructure::form_standard_tetrahedron_NA(){
		_metal_atom=atom("HETATM", 1649, "NA", " ", "NA", "A", 111, 19.287, 42.742, 64.785, 0.50, 10.50, 0.0, "NA", 0.0);
		_tetra_atoms[0]=atom("ATOM", 504, "O", " ", "THR", "A", 63, 20.933, 40.587, 64.525, 1.00, 55.15, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 511, "O", " ", "VAL", "A", 64, 18.312, 41.081, 62.856, 1.00, 55.70, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 1332, "O", " ", "THR", "B", 63, 20.964, 43.117, 62.669, 1.00, 61.76, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("ATOM", 1339, "O", " ", "VAL", "B", 64, 18.221, 44.710, 63.120, 1.00, 50.21, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_CA(){
		_metal_atom=atom("HETATM", 7170, "CA", " ", "CA", "F", 202, 4.120, 158.980, 194.134, 1.00, 86.45, 0.0, "CA", 0.0);
		_tetra_atoms[0]=atom("ATOM", 4693, "OD1", " ", "ASP", "F", 58, 2.351, 160.293, 195.274, 1.00, 88.49, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 4705, "OD1", " ", "ASN", "F", 60, 5.373, 158.883, 196.971, 1.00, 84.220, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 4759, "OE1", " ", "GLU", "F", 67, 4.283, 161.116, 192.198, 1.00, 60.04, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("HETATM", 7332, "O", " ", "HOH", "F", 301, 3.938, 161.798, 194.280, 1.00, 64.44, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_K(){
		_metal_atom=atom("HETATM", 2220, "K", " ", "K", "A", 3, -12.208, 39.862, 2.378, 1.00, 18.70, 0.0, "K", 0.0);
		_tetra_atoms[0]=atom("ATOM", 2108, "O", " ", "TYR", "A", 980, -12.164, 37.326, 3.317, 1.00, 13.95, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 2141, "O" , " ", "ARG", "A", 983, -9.532, 40.177, 2.775, 1.00, 17.94, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 2158, "OXT", " ", "LYS", "A", 984, -10.972, 41.427, 0.094, 1.00, 30.62, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("HETATM", 2269, "O", " ", "HOH", "A", 52, -12.105, 38.308, 0.040, 1.00, 26.87, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_MG(){
		_metal_atom=atom("HETATM", 5279, "MG", " ", "MG", "C", 201, -35.096, 28.802, -8.147, 1.00, 114.37, 0.0, "MG", 0.0);
		_tetra_atoms[0]=atom("ATOM", 2705, "O", " ", "ALA", "C", 43, -36.468, 31.011, -9.040, 1.00, 122.44, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 2726, "OD1", " ", "ASP", "C", 46, -35.736, 28.429, -10.645, 1.00, 187.540, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 2735, "O", " ", "GLN", "C", 48, -37.849, 28.410, -7.853, 1.00, 141.23, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("HETATM", 5284, "O", " ", "HOH", "C", 301, -33.086, 26.588, -8.119, 1.00, 100.360, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_MN(){
		_metal_atom=atom("HETATM", 4371, "MN", " ", "MN", "B", 503, 54.832, 3.667, 19.900, 1.00, 73.40, 0.0, "MN", 0.0);
		_tetra_atoms[0]=atom("ATOM", 3328, "OD1", " ", "ASN", "B", 354, 54.132, 2.467, 21.758, 1.00, 45.48, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 3427, "OD2", " ", "ASP", "B", 368, 54.887, 5.382, 21.371, 1.00, 57.37, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("HETATM", 4350, "O2A", " ", "ANP", "B", 502, 56.631, 5.424, 18.877, 1.00, 75.93, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("HETATM", 4458, "O", " ", "HOH", "B", 646, 52.677, 4.498, 18.678, 1.00, 53.66, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_FE(){
		_metal_atom=atom("HETATM", 8624, "FE", " ", "FE", "E", 201, 5.395, 38.352, -29.207, 1.00, 72.84, 0.0, "FE", 0.0);
		_tetra_atoms[0]=atom("ATOM", 5933, "OE1", " ", "GLU", "E", 27, 5.486, 37.570, -31.309, 1.00, 50.10, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 6239, "OE1", " ", "GLU", "E", 62, 7.474, 37.927, -29.520, 1.00, 65.31, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 6267, "ND1", " ", "HIS", "E", 65, 4.626, 35.993, -28.673, 1.00, 48.40, 0.0, "N", 0.0);
		_tetra_atoms[3]=atom("HETATM", 8697, "O", " ", "HOH", "E", 305, 3.284, 38.834, -29.839, 1.00, 63.89, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_CU(){
		_metal_atom=atom("HETATM", 2477, "CU", " ", "CU", "A", 1157, -1.317, -16.354, 16.473, 0.20, 9.45, 0.0, "CU", 0.0);
		_tetra_atoms[0]=atom("ATOM", 379, "ND1", " ", "HIS", "A", 46, -1.684, -16.116, 18.415, 1.00, 7.85, 0.0, "N", 0.0);
		_tetra_atoms[1]=atom("ATOM", 399, "NE2", " ", "HIS", "A", 48, -0.955, -14.399, 15.081, 1.00, 7.52, 0.0, "N", 0.0);
		_tetra_atoms[2]=atom("ATOM", 507, "NE2", " ", "HIS", "A", 63, 0.367, -17.755, 16.528, 1.00, 9.29, 0.0, "N", 0.0);
		_tetra_atoms[3]=atom("ATOM", 1008, "NE2", " ", "HIS", "A", 120, -3.547, -15.670, 16.025, 1.00, 7.52, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_ZN(){
		_metal_atom=atom("HETATM", 5244, "ZN", " ", "ZN", "C", 302, -4.393, 24.268, 71.872, 1.00, 12.70, 0.0, "ZN", 0.0);
		_tetra_atoms[0]=atom("ATOM", 3160, "NE2", " ", "HIS", "C", 168, -4.915, 24.963, 73.701, 1.00, 11.62, 0.0, "N", 0.0);
		_tetra_atoms[1]=atom("ATOM", 3172, "OD2", " ", "ASP", "C", 170, -5.350, 25.293, 70.535, 1.00, 11.70, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 3266, "NE2", " ", "HIS", "C", 183, -4.515, 22.315, 71.681, 1.00, 12.14, 0.0, "N", 0.0);
		_tetra_atoms[3]=atom("ATOM", 3343, "ND1", " ", "HIS", "C", 196, -2.557, 24.973, 71.331, 1.00, 13.35, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_CO(){
		_metal_atom=atom("HETATM", 3863, "CO", " ", "CO", "D", 201, -11.941, 42.346, 22.922, 1.00, 43.93, 0.0, "CO", 0.0);
		_tetra_atoms[0]=atom("ATOM", 2957, "NE2", " ", "HIS", "D", 50, -13.644, 40.620, 22.073, 1.00, 30.51, 0.0, "N", 0.0);
		_tetra_atoms[1]=atom("ATOM", 3474, "NE2", " ", "HIS", "D", 122, -12.844, 43.782, 21.344, 1.00, 23.69, 0.0, "N", 0.0);
		_tetra_atoms[2]=atom("ATOM", 3825, "OE1", " ", "GLU", "D", 172, -10.381, 42.060, 20.855, 1.00, 43.61, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("HETATM", 3938, "O", " ", "HOH", "D", 308, -13.130, 43.241, 24.742, 1.00, 27.24, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	/*void tetraMetalicStructure::form_standard_tetrahedron_MO(){
		_metal_atom=atom();
		_tetra_atoms[0]=atom();
		_tetra_atoms[1]=atom();
		_tetra_atoms[2]=atom();
		_tetra_atoms[3]=atom();

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_AG(){
		_metal_atom=atom();
		_tetra_atoms[0]=atom();
		_tetra_atoms[1]=atom();
		_tetra_atoms[2]=atom();
		_tetra_atoms[3]=atom();

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}*/

	void tetraMetalicStructure::form_standard_tetrahedron_PB(){
		_metal_atom=atom("HETATM", 2189, "PB", " ", "PB", "A", 201, 44.539, 66.520, 40.782, 0.50, 43.28, 0.0, "PB", 0.0);
		_tetra_atoms[0]=atom("ATOM", 334, "OD1", " ", "ASP", "A", 52, 43.017, 67.227, 43.212, 1.00, 39.20, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 335, "OD2", " ", "ASP", "A", 52, 43.470, 69.025, 42.027, 1.00, 48.03, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 1619, "OE1", " ", "GLU", "B", 75, 46.164, 68.564, 39.214, 1.00, 63.56, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("ATOM", 1620, "OE2", " ", "GLU", "B", 75, 46.678, 66.549, 38.450, 1.00, 68.85, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_HG(){
		_metal_atom=atom("HETATM", 2197, "HG", " ", "HG", "A", 901, 29.757, 12.850, 56.480, 0.80, 158.83, 0.0, "HG", 0.0);
		_tetra_atoms[0]=atom("ATOM", 527, "ND1", " ", "HIS", "A", 797, 30.455, 14.870, 54.971, 1.00, 43.36, 0.0, "N", 0.0);
		_tetra_atoms[1]=atom("ATOM", 895, "SG", " ", "CYS", "A", 848, 30.615, 11.003, 55.155, 1.00, 74.41, 0.0, "S", 0.0);
		_tetra_atoms[2]=atom("ATOM", 904, "O", " ", "SER", "A", 850, 30.139, 14.261, 59.302, 1.00, 49.97, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("ATOM", 912, "SG", " ", "CYS", "A", 851, 26.860, 12.237, 57.119, 1.00, 45.99, 0.0, "S", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_AS(){
		_metal_atom=atom("HETATM", 2458, "AS", " ", "ARS", "D", 502, 0.965, -34.313, 16.068, 1.00, 21.95, 0.0, "AS", 0.0);
		_tetra_atoms[0]=atom("ATOM", 2030, "OE1", " ", "GLU", "D", 24, 2.674, -33.575, 15.205, 1.00, 21.08, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 2290, "OD2", " ", "ASP", "D", 60, 0.043, -34.732, 14.318, 1.00, 21.47, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("HETATM", 2818, "O", " ", "HOH", "D", 591, -0.634, -33.109, 17.001, 1.00, 4.03, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("HETATM", 2819, "O", " ", "HOH", "D", 592, 1.525, -36.137, 17.274, 1.00, 3.22, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_CD(){
		_metal_atom=atom("HETATM", 2845, "CD", " ", "CD", "A", 1168, 26.348, 14.993, 4.272, 0.99, 7.56, 0.0, "CD", 0.0);
		_tetra_atoms[0]=atom("ATOM", 405, "OE1", " ", "GLU", "A", 26, 24.966, 13.332, 4.289, 1.00, 9.60, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 406, "OE2", " ", "GLU", "A", 26, 24.211, 14.323, 2.480, 1.00, 9.94, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 2089, "SG", " ", "CYS", "A", 133, 25.106, 16.905, 5.340, 1.00, 7.54, 0.0, "S", 0.0);
		_tetra_atoms[3]=atom("HETATM", 3085, "O", " ", "HOH", "A", 2239, 27.902, 13.972, 5.913, 1.00, 6.40, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	void tetraMetalicStructure::form_standard_tetrahedron_NI(){
		_metal_atom=atom("HETATM", 3349, "NI", " ", "NI", "A", 101, 55.743, 80.998, 54.823, 1.00, 103.11, 0.0, "NI", 0.0);
		_tetra_atoms[0]=atom("ATOM", 1, "N", " ", "GLY", "A", 1, 55.884, 80.846, 56.936, 1.00, 101.39, 0.0, "N", 0.0);
		_tetra_atoms[1]=atom("ATOM", 4, "O", " ", "GLY", "A", 1, 53.921, 79.889, 55.007, 1.00, 94.91, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 338, "NE2", " ", "HIS", "A", 43, 56.964, 79.202, 54.495, 1.00, 86.10, 0.0, "N", 0.0);
		_tetra_atoms[3]=atom("ATOM", 348, "NE2", " ", "HIS", "A", 44, 55.579, 81.175, 52.644, 1.00, 75.72, 0.0, "N", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	/*void tetraMetalicStructure::form_standard_tetrahedron_CR(){
		_metal_atom=atom();
		_tetra_atoms[0]=atom();
		_tetra_atoms[1]=atom();
		_tetra_atoms[2]=atom();
		_tetra_atoms[3]=atom();

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}*/

	void tetraMetalicStructure::form_standard_tetrahedron_CS(){
		_metal_atom=atom("HETATM", 12920, "CS", " ", "CS", "C", 402, 26.274, 9.205, -3.712, 1.00, 202.08, 0.0, "CS", 0.0);
		_tetra_atoms[0]=atom("HETATM", 6245, "O", " ", "MSE", "C", 163, 23.859, 7.212, -3.027, 1.00, 75.68, 0.0, "O", 0.0);
		_tetra_atoms[1]=atom("ATOM", 6267, "O", " ", "SER", "C", 165, 27.815, 6.466, -3.353, 1.00, 101.78, 0.0, "O", 0.0);
		_tetra_atoms[2]=atom("ATOM", 6282, "O", " ", "GLY", "C", 167, 28.083, 8.650, -6.566, 1.00, 71.71, 0.0, "O", 0.0);
		_tetra_atoms[3]=atom("HETATM", 13111, "O", " ", "HOH", "C", 543, 29.115, 10.652, -4.137, 1.00, 72.71, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_binding_atoms; i++){
			_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
		}

		_tetra=tetrahedron(_tetra_vertices);
		print(_metal_atom, _tetra_atoms);

		_tetra.measure_surface_areas_and_order_of_vertices();

		cout<<"\nVertices:";
		print(_tetra_vertices);
	}

	string tetraMetalicStructure::assign_standard_site_type(){
		string site_type;
		int aa=0, dna=0, rna=0, other=0, water=0;
		for (int x=0; x<_no_of_binding_atoms; x++){
			if (_tetra_atoms[x].get_molecule_type()=="amino_acid"){
				aa++;
			}else if (_tetra_atoms[x].get_molecule_type()=="dna_molecule"){
				dna++;
			}else if (_tetra_atoms[x].get_molecule_type()=="rna_molecule"){
				rna++;
			}else if (_tetra_atoms[x].get_molecule_type()=="other_biom"){
				other++;
			}else if (_tetra_atoms[x].get_residue_name()=="HOH"){
				water++;
			}
		}

		if (water==_no_of_binding_atoms){
			site_type="all_water";
		}else if ((aa!=0) && (dna==0) && (rna==0)){
			site_type="protein";
		}else if ((aa==0) && (dna!=0) && (rna==0)){
			site_type="dna";
		}else if ((aa==0) && (dna==0) && (rna!=0)){
			site_type="rna";
		}else if ((aa!=0) && (dna!=0) && (rna==0)){
			site_type="protein_dna_complex";
		}else if ((aa!=0) && (dna==0) && (rna!=0)){
			site_type="protein_rna_complex";
		}else if ((aa==0) && (dna!=0) && (rna!=0)){
			site_type="dna_rna_complex";
		}else if ((aa!=0) && (dna!=0) && (rna!=0)){
			site_type="protein_dna_rna_complex";
		}else{
			site_type="other";
		}

		return site_type;
	}


	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_NA(){
		_metal_atom=atom("HETATM", 10129, "NA", " ", "NA", "A", 3152, -24.117, 4.774, 4.369, 1.00, 32.14, 0.0, "NA", 0.0);
		_tribipyr_atoms_pl[0]=atom("HETATM", 10816, "O", " ", "HOH", "B", 4029, -25.396, 3.184, 4.888, 1.00, 27.39, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[1]=atom("HETATM", 10224, "O", " ", "HOH", "A", 3353, -23.876, 6.322, 2.800, 1.00, 26.53, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[2]=atom("HETATM", 10577, "O", " ", "HOH", "B", 3790, -22.136, 4.969, 5.897, 1.00, 20.42, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("ATOM", 1591, "O", " ", "ASN", "A", 190, -22.572, 3.010, 3.255, 1.00, 23.53, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[1]=atom("ATOM", 4125, "O", " ", "ASN", "B", 190, -24.921, 6.587, 5.935, 1.00, 22.63, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_CA(){
		_metal_atom=atom("HETATM", 5768, "CA", " ", "CA", "A", 302, -18.777, 49.103, 47.489, 1.00, 24.54, 0.0, "CA", 0.0);
		_tribipyr_atoms_pl[0]=atom("ATOM", 881, "OD1", " ", "ASN", "A", 114, -18.892, 46.918, 46.330, 1.00, 28.71, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[1]=atom("ATOM", 199, "O", " ", "GLY", "A", 25, -18.912, 49.750, 45.115, 1.00, 21.72, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[2]=atom("ATOM", 176, "O", " ", "PHE", "A", 23, -18.443, 51.207, 48.034, 1.00, 16.95, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("ATOM", 860, "O", " ", "TYR", "A", 112, -21.177, 49.359, 47.248, 1.00, 20.39, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[1]=atom("HETATM", 5989, "O", " ", "HOH", "A", 303, -19.939, 46.953, 49.666, 1.00, 6.00, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_MG(){
		_metal_atom=atom("HETATM", 5110, "MG", " ", "MG", "C", 2004, 54.733, 67.766, 5.798, 1.00, 26.17, 0.0, "MG", 0.0);
		_tribipyr_atoms_pl[0]=atom("HETATM", 6013, "O", " ", "HOH", "C", 2397, 54.686, 65.640, 5.507, 1.00, 33.03, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[1]=atom("HETATM", 5688, "O", " ", "HOH", "C", 2072, 52.640, 67.472, 5.910, 1.00, 27.15, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[2]=atom("HETATM", 5785, "O", " ", "HOH", "C", 2169, 54.499, 69.836, 5.774, 1.00, 20.75, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("HETATM", 6041, "O", " ", "HOH", "C", 2425, 53.373, 68.310, 4.151, 1.00, 27.30, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[1]=atom("HETATM", 5865, "O", " ", "HOH", "C", 2249, 54.646, 67.437, 7.801, 1.00, 25.22, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_MN(){
		_metal_atom=atom("HETATM", 2443, "MN", " ", "MN", "A", 2, -15.020, -16.928, 3.999, 0.50, 12.76, 0.0, "MN", 0.0);
		_tribipyr_atoms_pl[0]=atom("ATOM", 1357, "OD1", " ", "ASP", "A", 213, -17.205, -17.606, 3.407, 1.00, 13.93, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[1]=atom("ATOM", 1343, "OD2", " ", "ASP", "A", 211, -13.552, -18.556, 2.918, 1.00, 15.90, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[2]=atom("HETATM", 2468, "O3B", " ", "UDP", "A", 1, -13.672, -15.761, 5.130, 0.50, 15.92, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("HETATM", 2642, "O", " ", "HOH", "A", 524, -14.813, -15.377, 2.347, 1.00, 27.80, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[1]=atom("HETATM", 2462, "O1A", " ", "UDP", "A", 1, -15.395, -18.241, 5.739, 0.50, 13.23, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_FE(){
		_metal_atom=atom("HETATM", 4518, "FE", " ", "HEM", "D", 147, -19.888, -6.918, 94.231, 1.00, 16.87, 0.0, "FE", 0.0);
		_tribipyr_atoms_pl[0]=atom("ATOM", 3965, "NE2", " ", "HIS", "D", 92, -21.746, -6.234, 93.538, 1.00, 12.67, 0.0, "N", 0.0);
		_tribipyr_atoms_pl[1]=atom("HETATM", 4550, "ND", " ", "HEM", "D", 147, -19.113, -7.395, 92.467, 1.00, 16.63, 0.0, "N", 0.0);
		_tribipyr_atoms_pl[2]=atom("HETATM", 4534, "NB", " ", "HEM", "D", 147, -20.304, -6.489, 96.140, 1.00, 16.43, 0.0, "N", 0.0);
		_tribipyr_atoms_npl[0]=atom("HETATM", 4542, "NC", " ", "HEM", "D", 147, -20.084, -8.798, 94.645, 1.00, 17.26, 0.0, "N", 0.0);
		_tribipyr_atoms_npl[1]=atom("HETATM", 4523, "NA", " ", "HEM", "D", 147, -18.882, -5.183, 94.155, 1.00, 17.12, 0.0, "N", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_ZN(){
		_metal_atom=atom("HETATM", 1348, "ZN", " ", "ZN", "A", 999, 25.615, 56.463, 53.484, 1.00, 8.66, 0.0, "ZN", 0.0);
		_tribipyr_atoms_pl[0]=atom("ATOM", 957, "NE2", " ", "HIS", "A", 197, 24.914, 55.390, 51.933, 1.00, 7.19, 0.0, "N", 0.0);
		_tribipyr_atoms_pl[1]=atom("ATOM", 1036, "NE2", " ", "HIS", "A", 207, 27.004, 55.452, 54.523, 1.00, 8.24, 0.0, "N", 0.0);
		_tribipyr_atoms_pl[2]=atom("HETATM", 1341, "O1P", " ", "PAT", "B", 3, 24.565, 58.814, 52.829, 1.00, 32.21, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("ATOM", 991, "NE2", " ", "HIS", "A", 201, 24.096, 56.437, 54.939, 1.00, 7.40, 0.0, "N", 0.0);
		_tribipyr_atoms_npl[1]=atom("HETATM", 1342, "O2P", " ", "PAT", "B", 3, 26.751, 57.939, 53.331, 1.00, 29.06, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_CD(){
		_metal_atom=atom("HETATM", 1650, "CD", " ", "CD", "A", 901, 0.670, -107.745, -64.615, 1.00, 29.85, 0.0, "CD", 0.0);
		_tribipyr_atoms_pl[0]=atom("HETATM", 1656, "O", " ", "HOH", "A", 904, -1.735, -106.855, -64.211, 1.00, 19.11, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[1]=atom("HETATM", 1655, "O", " ", "HOH", "A", 903, 0.312, -110.131, -64.892, 1.00, 17.56, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[2]=atom("ATOM", 646, "OD2", " ", "ASP", "A", 85, 1.916, -105.802, -64.458, 1.00, 4.62, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("ATOM", 645, "OD1", " ", "ASP", "A", 85, 2.887, -107.453, -65.376, 1.00, 2.65, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[1]=atom("HETATM", 1657, "O", " ", "HOH", "A", 905, -0.195, -106.505, -67.061, 1.00, 29.26, 0.0, "O", 0.0 );


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_CO(){
		_metal_atom=atom("HETATM", 2788, "CO", " ", "CO", "A", 481, 18.366, 25.342, 16.442, 1.00, 39.49, 0.0, "CO", 0.0);
		_tribipyr_atoms_pl[0]=atom("ATOM", 1623, "NE2", " ", "HIS", "A", 331, 16.181, 25.895, 15.799, 1.00, 17.39, 0.0, "N", 0.0);
		_tribipyr_atoms_pl[1]=atom("ATOM", 1090, "OD2", " ", "ASP", "A", 262, 19.283, 25.775, 14.628, 1.00, 24.64, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[2]=atom("ATOM", 1874, "OE2", " ", "GLU", "A", 364, 17.382, 24.685, 18.458, 1.00, 26.64, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("ATOM", 2632, "OE2", " ", "GLU", "A", 459, 19.111, 23.142, 16.277, 1.00, 20.20, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[1]=atom("HETATM", 3029, "O", " ", "HOH", "A", 698, 20.137, 26.070, 17.764, 1.00, 33.49, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_HG(){
		_metal_atom=atom("HETATM", 2228, "HG", " ", "HG", "A", 401, -19.501, -67.999, 15.601, 1.00, 18.04, 0.0, "HG", 0.0);
		_tribipyr_atoms_pl[0]=atom("ATOM", 141, "SG", " ", "CYS", "A", 80, -21.387, -67.614, 16.891, 1.00, 15.48, 0.0, "S", 0.0);
		_tribipyr_atoms_pl[1]=atom("HETATM", 2315, "O", " ", "HOH", "A", 557, -20.879, -70.442, 13.587, 1.00, 17.73, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[2]=atom("HETATM", 2583, "O", " ", "HOH", "A", 825, -16.537, -66.580, 15.390, 1.00, 25.86, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("HETATM", 2524, "O", " ", "HOH", "A", 766, -17.757, -68.489, 14.494, 1.00, 22.55, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[1]=atom("ATOM", 288, "O", " ", "GLY", "A", 98, -17.921, -68.084, 18.093, 1.00, 14.01, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	void tribipyrMetalicStructure::form_standard_trigonal_bipyramid_NI(){
		_metal_atom=atom("HETATM", 4817, "NI", " ", "NI", "D", 203, 9.472, 21.852, 34.540, 1.00, 46.05, 0.0, "NI", 0.0);
		_tribipyr_atoms_pl[0]=atom("HETATM", 5396, "O", " ", "HOH", "D", 2093, 8.462, 23.907, 34.372, 1.00, 38.96, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[1]=atom("HETATM", 5446, "O", " ", "HOH", "D2", 143, 8.806, 20.627, 32.828, 1.00, 46.40, 0.0, "O", 0.0);
		_tribipyr_atoms_pl[2]=atom("HETATM", 5394, "O", " ", "HOH", "D2", 91, 9.988, 19.959, 34.894, 1.00, 44.75, 0.0, "O", 0.0);
		_tribipyr_atoms_npl[0]=atom("ATOM", 4178, "NE2", " ", "HIS", "D", 103, 10.595, 22.632, 36.088, 1.00, 40.19, 0.0, "N", 0.0);
		_tribipyr_atoms_npl[1]=atom("HETATM", 5395, "O", " ", "HOH", "D2", 92, 10.578, 22.839, 33.465, 1.00, 34.80, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_tribipyr_vertices_pl[i]=Point3D(_tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_tribipyr_vertices_npl[i]=Point3D(_tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_tribipyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tribipyr_atoms_pl[2].get_serial_no());

		_non_planar_nos.push_back(_tribipyr_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_tribipyr_atoms_npl[1].get_serial_no());

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);
	}

	string tribipyrMetalicStructure::assign_standard_site_type(){
		string site_type;
		int aa=0, dna=0, rna=0, other=0, water=0;
		for (int x=0; x<_no_of_planar_atoms; x++){
			if (_tribipyr_atoms_pl[x].get_molecule_type()=="amino_acid"){
				aa++;
			}else if (_tribipyr_atoms_pl[x].get_molecule_type()=="dna_molecule"){
				dna++;
			}else if (_tribipyr_atoms_pl[x].get_molecule_type()=="rna_molecule"){
				rna++;
			}else if (_tribipyr_atoms_pl[x].get_molecule_type()=="other_biom"){
				other++;
			}else if (_tribipyr_atoms_pl[x].get_residue_name()=="HOH"){
				water++;
			}
		}

		for (int x=0; x<_no_of_non_planar_atoms; x++){
			if (_tribipyr_atoms_npl[x].get_molecule_type()=="amino_acid"){
				aa++;
			}else if (_tribipyr_atoms_npl[x].get_molecule_type()=="dna_molecule"){
				dna++;
			}else if (_tribipyr_atoms_npl[x].get_molecule_type()=="rna_molecule"){
				rna++;
			}else if (_tribipyr_atoms_npl[x].get_molecule_type()=="other_biom"){
				other++;
			}else if (_tribipyr_atoms_npl[x].get_residue_name()=="HOH"){
				water++;
			}
		}

		if (water==_no_of_binding_atoms){
			site_type="all_water";
		}else if ((aa!=0) && (dna==0) && (rna==0)){
			site_type="protein";
		}else if ((aa==0) && (dna!=0) && (rna==0)){
			site_type="dna";
		}else if ((aa==0) && (dna==0) && (rna!=0)){
			site_type="rna";
		}else if ((aa!=0) && (dna!=0) && (rna==0)){
			site_type="protein_dna_complex";
		}else if ((aa!=0) && (dna==0) && (rna!=0)){
			site_type="protein_rna_complex";
		}else if ((aa==0) && (dna!=0) && (rna!=0)){
			site_type="dna_rna_complex";
		}else if ((aa!=0) && (dna!=0) && (rna!=0)){
			site_type="protein_dna_rna_complex";
		}else{
			site_type="other";
		}

		return site_type;
	}

	void tetrapyrMetalicStructure::form_standard_tetragonal_pyramid_NA(){
		cout<<"\nStandard Started";
		_metal_atom=atom("HETATM", 2828, "NA", "A", "NA", "A", 203, -5.473, -8.260, -11.029, 0.55, 16.96, 0.0, "NA", 0.0);
		_tetrapyr_atoms_pl[0]=atom("ATOM", 485, "OE2", "A", "GLU", "A", 56, -5.046, -9.322, -8.890, 0.25, 25.88, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[1]=atom("HETATM", 2801, "O1B", " ", "6U4", "A", 201, -5.923, -5.483, -9.601, 1.00, 17.67, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[2]=atom("HETATM", 2805, "O2A", " ", "6U4", "A", 201, -6.096, -5.918, -12.705, 1.00, 20.92, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[3]=atom("HETATM", 2944, "O", " ", "HOH", "A", 382, -5.315, -10.107, -12.368, 1.00, 32.01, 0.0, "O", 0.0);
		_tetrapyr_atom_npl=atom("ATOM", 311, "O", " ", "GLY", "A", 36, -2.955, -8.732, -10.973, 1.00, 12.54, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_planar_atoms; i++){
			_tetrapyr_vertices_pl[i]=Point3D(_tetrapyr_atoms_pl[i].X(), _tetrapyr_atoms_pl[i].Y(), _tetrapyr_atoms_pl[i].Z());
		}

		_tetrapyr_vertex_npl=Point3D(_tetrapyr_atom_npl.X(), _tetrapyr_atom_npl.Y(), _tetrapyr_atom_npl.Z());

		_planar_nos.push_back(_tetrapyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[3].get_serial_no());
		_non_planar_no.push_back(_tetrapyr_atom_npl.get_serial_no());

		cout<<"\nStandard Mid";
cout<<"\nForm Standard Struct";
		_tetrapyr=tetragonal_pyramidal(_tetrapyr_vertices_pl, _tetrapyr_vertex_npl);
		cout<<"\nEnd Forming Standard Struct";
		print(_metal_atom, _tetrapyr_atoms_pl, _no_of_planar_atoms);
		print(_tetrapyr_atom_npl, _no_of_non_planar_atom);

		cout<<"\nVertices:";
		print(_tetrapyr_vertices_pl, _no_of_planar_atoms);
		print(_tetrapyr_vertex_npl, _no_of_non_planar_atom);
		cout<<"\nStandard End";
	}

	void tetrapyrMetalicStructure::form_standard_tetragonal_pyramid_CA(){
		cout<<"\nStandard Started";
		_metal_atom=atom("HETATM", 4875, "CA", " ", "CA", "H", 521, 18.840, -16.752, 38.793, 0.43, 18.84, 0.0, "CA", 0.0);
		_tetrapyr_atoms_pl[0]=atom("ATOM", 3305, "O", " ", "LYS", "H", 169, 17.270, -14.971, 38.649, 1.00, 25.67, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[1]=atom("ATOM", 3350, "O", " ", "THR", "H", 172, 18.090, -17.476, 36.795, 1.00, 23.37, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[2]=atom("HETATM", 5299, "O", " ", "HOH", "H", 519, 20.300, -18.556, 39.302, 1.00, 34.66, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[3]=atom("HETATM", 5599, "O", " ", "HOH", "H", 747, 19.276, -15.941, 40.979, 1.00, 37.38, 0.0, "O", 0.0);
		_tetrapyr_atom_npl=atom("HETATM", 5305, "O", " ", "HOH", "H", 522, 20.147, -15.226, 37.426, 1.00, 34.48, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_planar_atoms; i++){
			_tetrapyr_vertices_pl[i]=Point3D(_tetrapyr_atoms_pl[i].X(), _tetrapyr_atoms_pl[i].Y(), _tetrapyr_atoms_pl[i].Z());
		}

		_tetrapyr_vertex_npl=Point3D(_tetrapyr_atom_npl.X(), _tetrapyr_atom_npl.Y(), _tetrapyr_atom_npl.Z());

		_planar_nos.push_back(_tetrapyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[3].get_serial_no());
		_non_planar_no.push_back(_tetrapyr_atom_npl.get_serial_no());

		cout<<"\nStandard Mid";
cout<<"\nForm Standard Struct";
		_tetrapyr=tetragonal_pyramidal(_tetrapyr_vertices_pl, _tetrapyr_vertex_npl);
		cout<<"\nEnd Forming Standard Struct";
		print(_metal_atom, _tetrapyr_atoms_pl, _no_of_planar_atoms);
		print(_tetrapyr_atom_npl, _no_of_non_planar_atom);

		cout<<"\nVertices:";
		print(_tetrapyr_vertices_pl, _no_of_planar_atoms);
		print(_tetrapyr_vertex_npl, _no_of_non_planar_atom);
		cout<<"\nStandard End";
	}

	void tetrapyrMetalicStructure::form_standard_tetragonal_pyramid_MG(){
		cout<<"\nStandard Started";
		_metal_atom=atom("HETATM", 1587, "MG", " ", "MG", "R", 171, 15.097, -25.469, -22.377, 0.50, 22.44, 0.0, "MG", 0.0);
		_tetrapyr_atoms_pl[0]=atom("ATOM", 929, "OD1", "A", "ASP", "R", 105, 16.375, -25.963, -23.822, 0.50, 22.16, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[1]=atom("HETATM", 1689, "O", " ", "HOH", "R", 273, 13.954, -27.445, -22.446, 0.50, 22.41, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[2]=atom("ATOM", 931, "OD2", "A", "ASP", "R", 105, 15.016, -27.655, -23.495, 0.50, 23.53, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[3]=atom("HETATM", 1687, "O", " ", "HOH", "R", 271, 15.939, -23.916, -22.499, 0.50, 25.92, 0.0, "O", 0.0);
		_tetrapyr_atom_npl=atom("ATOM", 888, "O", "A", "ARG", "R", 102, 13.446, -24.413, -23.301, 0.50, 20.97, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_planar_atoms; i++){
			_tetrapyr_vertices_pl[i]=Point3D(_tetrapyr_atoms_pl[i].X(), _tetrapyr_atoms_pl[i].Y(), _tetrapyr_atoms_pl[i].Z());
		}

		_tetrapyr_vertex_npl=Point3D(_tetrapyr_atom_npl.X(), _tetrapyr_atom_npl.Y(), _tetrapyr_atom_npl.Z());

		_planar_nos.push_back(_tetrapyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[3].get_serial_no());
		_non_planar_no.push_back(_tetrapyr_atom_npl.get_serial_no());

		cout<<"\nStandard Mid";
cout<<"\nForm Standard Struct";
		_tetrapyr=tetragonal_pyramidal(_tetrapyr_vertices_pl, _tetrapyr_vertex_npl);
		cout<<"\nEnd Forming Standard Struct";
		print(_metal_atom, _tetrapyr_atoms_pl, _no_of_planar_atoms);
		print(_tetrapyr_atom_npl, _no_of_non_planar_atom);

		cout<<"\nVertices:";
		print(_tetrapyr_vertices_pl, _no_of_planar_atoms);
		print(_tetrapyr_vertex_npl, _no_of_non_planar_atom);
		cout<<"\nStandard End";
	}

	void tetrapyrMetalicStructure::form_standard_tetragonal_pyramid_NI(){
		cout<<"\nStandard Started";
		_metal_atom=atom("HETATM", 1294, "NI", " ", "NI", "A", 320, 17.045, 2.879, 5.053, 1.00, 35.61, 0.0, "NI", 0.0);
		_tetrapyr_atoms_pl[0]=atom("HETATM", 1438, "O", " ", "HOH", "A", 464, 18.111, 3.138, 3.627, 1.00, 28.35, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[1]=atom("HETATM", 1378, "O", " ", "HOH", "A", 404, 15.829, 1.178, 4.229, 1.00, 41.16, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[2]=atom("ATOM", 4, "O", " ", "ALA", "A", 2, 18.107, 4.504, 5.864, 1.00, 30.78, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[3]=atom("ATOM", 1, "N", " ", "ALA", "A", 2, 15.878, 3.135, 6.828, 1.00, 39.35, 0.0, "N", 0.0);
		_tetrapyr_atom_npl=atom("HETATM", 1401, "O", " ", "HOH", "A", 427, 15.423, 4.168, 4.194, 1.00, 41.77, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_planar_atoms; i++){
			_tetrapyr_vertices_pl[i]=Point3D(_tetrapyr_atoms_pl[i].X(), _tetrapyr_atoms_pl[i].Y(), _tetrapyr_atoms_pl[i].Z());
		}

		_tetrapyr_vertex_npl=Point3D(_tetrapyr_atom_npl.X(), _tetrapyr_atom_npl.Y(), _tetrapyr_atom_npl.Z());

		_planar_nos.push_back(_tetrapyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[3].get_serial_no());
		_non_planar_no.push_back(_tetrapyr_atom_npl.get_serial_no());

		cout<<"\nStandard Mid";
cout<<"\nForm Standard Struct";
		_tetrapyr=tetragonal_pyramidal(_tetrapyr_vertices_pl, _tetrapyr_vertex_npl);
		cout<<"\nEnd Forming Standard Struct";
		print(_metal_atom, _tetrapyr_atoms_pl, _no_of_planar_atoms);
		print(_tetrapyr_atom_npl, _no_of_non_planar_atom);

		cout<<"\nVertices:";
		print(_tetrapyr_vertices_pl, _no_of_planar_atoms);
		print(_tetrapyr_vertex_npl, _no_of_non_planar_atom);
		cout<<"\nStandard End";
	}

	void tetrapyrMetalicStructure::form_standard_tetragonal_pyramid_CD(){
		cout<<"\nStandard Started";
		_metal_atom=atom("HETATM", 1651, "CD", " ", "CD", "B", 902, 13.019, -130.514, -47.301, 1.00, 31.40, 0.0, "CD", 0.0);
		_tetrapyr_atoms_pl[0]=atom("ATOM", 1200, "OE2", " ", "GLU", "B", 50, 11.105, -129.763, -46.370, 1.00, 3.86, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[1]=atom("HETATM", 1748, "O", " ", "HOH", "B", 950, 13.735, -132.006, -48.950, 1.00, 35.90, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[2]=atom("HETATM", 1716, "O", " ", "HOH", "B", 918, 15.193, -130.459, -46.743, 1.00, 21.33, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[3]=atom("ATOM", 1199, "OE1", "GLU", " ", "B", 50, 12.838, -129.110, -45.185, 1.00, 7.19, 0.0, "O", 0.0);
		_tetrapyr_atom_npl=atom("HETATM", 1709, "O", " ", "HOH", "B", 911, 11.312, -131.617, -47.901, 1.00, 28.28, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_planar_atoms; i++){
			_tetrapyr_vertices_pl[i]=Point3D(_tetrapyr_atoms_pl[i].X(), _tetrapyr_atoms_pl[i].Y(), _tetrapyr_atoms_pl[i].Z());
		}

		_tetrapyr_vertex_npl=Point3D(_tetrapyr_atom_npl.X(), _tetrapyr_atom_npl.Y(), _tetrapyr_atom_npl.Z());

		_planar_nos.push_back(_tetrapyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[3].get_serial_no());
		_non_planar_no.push_back(_tetrapyr_atom_npl.get_serial_no());

		cout<<"\nStandard Mid";
cout<<"\nForm Standard Struct";
		_tetrapyr=tetragonal_pyramidal(_tetrapyr_vertices_pl, _tetrapyr_vertex_npl);
		cout<<"\nEnd Forming Standard Struct";
		print(_metal_atom, _tetrapyr_atoms_pl, _no_of_planar_atoms);
		print(_tetrapyr_atom_npl, _no_of_non_planar_atom);

		cout<<"\nVertices:";
		print(_tetrapyr_vertices_pl, _no_of_planar_atoms);
		print(_tetrapyr_vertex_npl, _no_of_non_planar_atom);
		cout<<"\nStandard End";
	}

	void tetrapyrMetalicStructure::form_standard_tetragonal_pyramid_ZN(){
		cout<<"\nStandard Started";
		_metal_atom=atom("HETATM", 2416, "ZN", " ", "ZN", "A", 1157, 7.824, -14.869, 11.431, 0.15, 10.44, 0.0, "ZN", 0.0);
		_tetrapyr_atoms_pl[0]=atom("HETATM", 2444, "O4", " ", "SO4", "A", 1163, 9.179, -15.371, 10.262, 0.40, 13.67, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[1]=atom("ATOM", 1152, "N", " ", "ARG", "A", 143, 9.404, -15.516, 13.082, 1.00, 8.90, 0.0, "N", 0.0);
		_tetrapyr_atoms_pl[2]=atom("HETATM", 2706, "O", " ", "HOH", "A", 2245, 9.710, -16.182, 10.322, 0.60, 14.65, 0.0, "O", 0.0);
		_tetrapyr_atoms_pl[3]=atom("ATOM", 1155, "O", " ", "ARG", "A", 143, 8.214, -12.821, 12.522, 1.00, 12.17, 0.0, "O", 0.0);
		_tetrapyr_atom_npl=atom("HETATM", 2690, "O", " ", "HOH", "A", 2229, 6.836, -14.598, 9.856, 1.00, 30.58, 0.0, "O", 0.0);

		for (int i=0; i<_no_of_planar_atoms; i++){
			_tetrapyr_vertices_pl[i]=Point3D(_tetrapyr_atoms_pl[i].X(), _tetrapyr_atoms_pl[i].Y(), _tetrapyr_atoms_pl[i].Z());
		}

		_tetrapyr_vertex_npl=Point3D(_tetrapyr_atom_npl.X(), _tetrapyr_atom_npl.Y(), _tetrapyr_atom_npl.Z());

		_planar_nos.push_back(_tetrapyr_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_tetrapyr_atoms_pl[3].get_serial_no());
		_non_planar_no.push_back(_tetrapyr_atom_npl.get_serial_no());

		cout<<"\nStandard Mid";
cout<<"\nForm Standard Struct";
		_tetrapyr=tetragonal_pyramidal(_tetrapyr_vertices_pl, _tetrapyr_vertex_npl);
		cout<<"\nEnd Forming Standard Struct";
		print(_metal_atom, _tetrapyr_atoms_pl, _no_of_planar_atoms);
		print(_tetrapyr_atom_npl, _no_of_non_planar_atom);

		cout<<"\nVertices:";
		print(_tetrapyr_vertices_pl, _no_of_planar_atoms);
		print(_tetrapyr_vertex_npl, _no_of_non_planar_atom);
		cout<<"\nStandard End";
	}

	string tetrapyrMetalicStructure::assign_standard_site_type(){
		string site_type;
		int aa=0, dna=0, rna=0, other=0, water=0;
		for (int x=0; x<_no_of_planar_atoms; x++){
			if (_tetrapyr_atoms_pl[x].get_molecule_type()=="amino_acid"){
				aa++;
			}else if (_tetrapyr_atoms_pl[x].get_molecule_type()=="dna_molecule"){
				dna++;
			}else if (_tetrapyr_atoms_pl[x].get_molecule_type()=="rna_molecule"){
				rna++;
			}else if (_tetrapyr_atoms_pl[x].get_molecule_type()=="other_biom"){
				other++;
			}else if (_tetrapyr_atoms_pl[x].get_residue_name()=="HOH"){
				water++;
			}
		}

		if (_tetrapyr_atom_npl.get_molecule_type()=="amino_acid"){
			aa++;
		}else if (_tetrapyr_atom_npl.get_molecule_type()=="dna_molecule"){
			dna++;
		}else if (_tetrapyr_atom_npl.get_molecule_type()=="rna_molecule"){
			rna++;
		}else if (_tetrapyr_atom_npl.get_molecule_type()=="other_biom"){
			other++;
		}else if (_tetrapyr_atom_npl.get_residue_name()=="HOH"){
			water++;
		}

		if (water==_no_of_binding_atoms){
			site_type="all_water";
		}else if ((aa!=0) && (dna==0) && (rna==0)){
			site_type="protein";
		}else if ((aa==0) && (dna!=0) && (rna==0)){
			site_type="dna";
		}else if ((aa==0) && (dna==0) && (rna!=0)){
			site_type="rna";
		}else if ((aa!=0) && (dna!=0) && (rna==0)){
			site_type="protein_dna_complex";
		}else if ((aa!=0) && (dna==0) && (rna!=0)){
			site_type="protein_rna_complex";
		}else if ((aa==0) && (dna!=0) && (rna!=0)){
			site_type="dna_rna_complex";
		}else if ((aa!=0) && (dna!=0) && (rna!=0)){
			site_type="protein_dna_rna_complex";
		}else{
			site_type="other";
		}

		return site_type;
	}

	void octaMetalicStructure::form_standard_octa_NA(){
		_metal_atom=atom("HETATM", 2453, "NA", " ", "NA", "A", 202, 16.836, -0.925, 13.486, 1.00, 17.41, 0.0, "NA", 0.0);
		_octa_atoms_pl[0]=atom("HETATM", 2440, "O6", " ", "BCN", "A", 201, 15.263, 0.068, 12.444, 1.00, 17.83, 0.0, "O", 0.0);
		_octa_atoms_pl[1]=atom("HETATM", 2549, "O", " ", "HOH", "A", 396, 15.390, -1.929, 14.684, 1.00, 18.54, 0.0, "O", 0.0);
		_octa_atoms_pl[2]=atom("HETATM", 2585, "O", " ", "HOH", "A", 432, 18.259, 0.020, 12.264, 1.00, 18.95, 0.0, "O", 0.0);
		_octa_atoms_pl[3]=atom("HETATM", 2649, "O", " ", "HOH", "A", 496, 18.425, -1.928, 14.497, 1.00, 19.63, 0.0, "O", 0.0);
		_octa_atoms_npl[0]=atom("HETATM", 2484, "O", " ", "HOH", "A", 331, 16.844, 0.706, 14.755, 1.00, 17.65, 0.0, "O", 0.0);
		_octa_atoms_npl[1]=atom("HETATM", 2518, "O", " ", "HOH", "A", 365, 16.658, -2.505, 12.105, 1.00, 19.25, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	void octaMetalicStructure::form_standard_octa_CA(){
		_metal_atom=atom("HETATM", 2727, "CA", " ", "CA", "A", 303, -12.237, 26.114, 34.655, 1.00, 6.07, 0.0, "CA", 0.0);
		_octa_atoms_pl[0]=atom("ATOM", 648, "O", " ", "GLY", "A", 161, -13.659, 24.445, 33.496, 1.00, 5.78, 0.0, "O", 0.0);
		_octa_atoms_pl[1]=atom("ATOM", 632, "OD1", " ", "ASP", "A", 158, -12.677, 27.747, 32.850, 1.00, 2.00, 0.0, "O", 0.0);
		_octa_atoms_pl[2]=atom("ATOM", 792, "OD2", " ", "ASP", "A", 181, -10.864, 27.598, 35.732, 1.00, 6.98, 0.0, "O", 0.0);
		_octa_atoms_pl[3]=atom("ATOM", 817, "OE2", " ", "GLU", "A", 184, -11.831, 24.685, 36.366, 1.00, 8.95, 0.0, "O", 0.0);
		_octa_atoms_npl[0]=atom("ATOM", 637, "O", " ", "GLY", "A", 159, -13.925, 27.086, 35.858, 1.00, 8.17, 0.0, "O", 0.0);
		_octa_atoms_npl[1]=atom("ATOM", 660, "O", " ", "VAL", "A", 163, -10.582, 25.122, 33.219, 1.00, 2.12, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	void octaMetalicStructure::form_standard_octa_MG(){
		_metal_atom=atom("HETATM", 1377, "MG", " ", "MG", "A", 168, 4.626, 33.115, 18.949, 1.00, 7.82, 0.0, "MG", 0.0);
		_octa_atoms_pl[0]=atom("ATOM", 273, "OG1", " ", "THR", "A", 35, 2.905, 34.282, 18.998, 1.00, 10.25, 0.0, "O", 0.0);
		_octa_atoms_pl[1]=atom("HETATM", 1411, "O", " ", "HOH", "A", 1002, 3.633, 31.646, 17.915, 1.00, 9.42, 0.0, "O", 0.0);
		_octa_atoms_pl[2]=atom("HETATM", 1385, "O2B", " ", "GNP", "A", 167, 6.425, 32.127, 18.752, 1.00, 8.43, 0.0, "O", 0.0);
		_octa_atoms_pl[3]=atom("HETATM", 1410, "O", " ", "HOH", "A", 1001, 5.769, 34.698, 19.719, 1.00, 8.70, 0.0, "O", 0.0);
		_octa_atoms_npl[0]=atom("ATOM", 122, "OG", " ", "SER", "A", 17, 4.876, 34.000, 17.065, 1.00, 7.49, 0.0, "O", 0.0);
		_octa_atoms_npl[1]=atom("HETATM", 1380, "O2G", " ", "GNP", "A", 167, 4.129, 32.343, 20.710, 1.00, 9.01, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	void octaMetalicStructure::form_standard_octa_MN(){
		_metal_atom=atom("HETATM", 4506, "MN", " ", "MN", "A", 301, 0.055, -2.576, 7.483, 1.00, 9.02, 0.0, "MN", 0.0);
		_octa_atoms_pl[0]=atom("HETATM", 4509, "O2", " ", "AKG", "A", 302, 1.136, -0.906, 6.890, 1.00, 10.07, 0.0, "O", 0.0);
		_octa_atoms_pl[1]=atom("ATOM", 833, "NE2", " ", "HIS", "A", 121, -0.751, -2.951, 5.606, 1.00, 11.54, 0.0, "N", 0.0);
		_octa_atoms_pl[2]=atom("ATOM", 1259, "NE2", " ", "HIS", "A", 177, -1.225, -4.115, 8.265, 1.00, 8.95, 0.0, "N", 0.0);
		_octa_atoms_pl[3]=atom("HETATM", 4579, "O", " ", "HOH", "A", 437, 0.782, -2.110, 9.500, 1.00, 10.94, 0.0, "O", 0.0);
		_octa_atoms_npl[0]=atom("HETATM", 4511, "O5", " ", "AKG", "A", 302, -1.323, -1.025, 8.005, 1.00, 9.73, 0.0, "O", 0.0);
		_octa_atoms_npl[1]=atom("ATOM", 847, "OD1", " ", "ASP", "A", 123, 1.598, -4.034, 7.149, 1.00, 9.59, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	void octaMetalicStructure::form_standard_octa_FE(){
		_metal_atom=atom("HETATM", 4347, "FE", " ", "HEM", "A", 142, 48.284, 27.721, 30.645, 1.00, 11.81, 0.0, "FE", 0.0);
		_octa_atoms_pl[0]=atom("HETATM", 4379, "ND", " ", "HEM", "A", 142, 47.932, 26.229, 29.365, 1.00, 9.49, 0.0, "N", 0.0);
		_octa_atoms_pl[1]=atom("HETATM", 4352, "NA", " ", "HEM", "A", 142, 49.218, 26.419, 31.858, 1.00, 11.70, 0.0, "N", 0.0);
		_octa_atoms_pl[2]=atom("HETATM", 4363, "NB", " ", "HEM", "A", 142, 48.723, 29.217, 31.914, 1.00, 14.03, 0.0, "N", 0.0);
		_octa_atoms_pl[3]=atom("HETATM", 4371, "NC", " ", "HEM", "A", 142, 47.481, 29.046, 29.390, 1.00, 8.20, 0.0, "N", 0.0);
		_octa_atoms_npl[0]=atom("HETATM", 4390, "N", " ", "NO", "A", 143, 46.695, 27.446, 31.313, 1.00, 15.20, 0.0, "N", 0.0);
		_octa_atoms_npl[1]=atom("ATOM", 650, "NE2", " ", "HIS", "A", 87, 50.301, 27.991, 29.656, 1.00, 12.18, 0.0, "N", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	void octaMetalicStructure::form_standard_octa_ZN(){
		_metal_atom=atom("HETATM", 2538, "ZN", "A", "ZN", "A", 206, -26.103, 52.052, 2.761, 0.40, 16.93, 0.0, "ZN", 0.0);
		_octa_atoms_pl[0]=atom("ATOM", 1813, "OD1", " ", "ASP", "A", 147, -28.069, 51.656, 3.335, 1.00, 11.59, 0.0, "O", 0.0);
		_octa_atoms_pl[1]=atom("ATOM", 2371, "ND1", " ", "HIS", "A", 182, -25.094, 50.219, 3.603, 1.00, 11.93, 0.0, "N", 0.0);
		_octa_atoms_pl[2]=atom("HETATM", 2560, "O", " ", "HOH", "A", 319, -24.100, 52.627, 2.150, 1.00, 12.83, 0.0, "O", 0.0);
		_octa_atoms_pl[3]=atom("HETATM", 2554, "O", " ", "HOH", "A", 313, -26.670, 53.858, 1.906, 1.00, 13.02, 0.0, "O", 0.0);
		_octa_atoms_npl[0]=atom("ATOM", 1706, "NE2", " ", "HIS", "A", 140, -26.369, 50.968, 0.817, 1.00, 12.66, 0.0, "N", 0.0);
		_octa_atoms_npl[1]=atom("ATOM", 2278, "OE1", " ", "GLU", "A", 176, -25.682, 53.146, 4.584, 1.00, 12.33, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	void octaMetalicStructure::form_standard_octa_CO(){
		_metal_atom=atom("HETATM", 2339, "CO", " ", "NCO", "B", 201, 37.333, 84.916, 79.444, 0.80, 66.43, 0.0, "CO", 0.0);
		_octa_atoms_pl[0]=atom("HETATM", 2344, "N5", " ", "NCO", "B", 201, 35.472, 84.331, 79.366, 0.80, 67.59, 0.0, "N", 0.0);
		_octa_atoms_pl[1]=atom("HETATM", 2345, "N6", " ", "NCO", "B", 201, 37.866, 82.982, 79.402, 0.80, 68.59, 0.0, "N", 0.0);
		_octa_atoms_pl[2]=atom("HETATM", 2340, "N1", " ", "NCO", "B", 201, 39.240, 85.464, 79.436, 0.80, 69.14, 0.0, "N", 0.0);
		_octa_atoms_pl[3]=atom("HETATM", 2341, "N2", " ", "NCO", "B", 201, 36.805, 86.826, 79.391, 0.80, 68.05, 0.0, "N", 0.0);
		_octa_atoms_npl[0]=atom("HETATM", 2342, "N3", " ", "NCO", "B", 201, 37.428, 84.923, 77.439, 0.80, 68.00, 0.0, "N", 0.0);
		_octa_atoms_npl[1]=atom("HETATM", 2343, "N4", " ", "NCO", "B", 201, 37.114, 84.992, 81.366, 0.80, 68.47, 0.0, "N", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	void octaMetalicStructure::form_standard_octa_CD(){
		_metal_atom=atom("HETATM", 3319, "CD", " ", "CD", "L", 212, 46.051, 20.009, 8.855, 1.00, 39.64, 0.0, "CD", 0.0);
		_octa_atoms_pl[0]=atom("HETATM", 3348, "O", " ", "HOH", "L", 241, 45.687, 18.315, 7.132, 1.00, 31.63, 0.0, "O", 0.0);
		_octa_atoms_pl[1]=atom("HETATM", 3350, "O", " ", "HOH", "L", 243, 46.278, 18.291, 10.406, 1.00, 29.77, 0.0, "O", 0.0);
		_octa_atoms_pl[2]=atom("HETATM", 3349, "O", " ", "HOH", "L", 242, 46.777, 21.524, 10.573, 1.00, 26.08, 0.0, "O", 0.0);
		_octa_atoms_pl[3]=atom("ATOM", 1072, "ND2", " ", "ASN", "L", 138, 46.399, 21.864, 7.652, 1.00, 32.84, 0.0, "N", 0.0);
		_octa_atoms_npl[0]=atom("ATOM", 2942, "NE2", " ", "HIS", "H", 164, 43.630, 20.825, 10.406, 1.00, 30.38, 0.0, "N", 0.0);
		_octa_atoms_npl[1]=atom("HETATM", 3411, "O", " ", "HOH", "L", 304, 48.623, 19.584, 9.464, 1.00, 50.03, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	void octaMetalicStructure::form_standard_octa_NI(){
		_metal_atom=atom("HETATM", 3727, "NI", " ", "NI", "B", 970, -9.163, 35.799, 14.450, 0.50, 38.46, 0.0, "NI", 0.0);
		_octa_atoms_pl[0]=atom("HETATM", 4181, "O", " ", "HOH", "B", 2043, -11.023, 36.677, 14.878, 0.50, 31.54, 0.0, "O", 0.0);
		_octa_atoms_pl[1]=atom("ATOM", 2075, "NE2", " ", "HIS", "B", 735, -9.988, 34.548, 12.877, 1.00, 22.98, 0.0, "N", 0.0);
		_octa_atoms_pl[2]=atom("HETATM", 3831, "O", " ", "HOH", "A", 2096, -7.286, 34.893, 14.204, 0.50, 20.73, 0.0, "O", 0.0);
		_octa_atoms_pl[3]=atom("HETATM", 3832, "O", " ", "HOH", "A", 2097, -8.311, 36.900, 16.027, 0.50, 30.25, 0.0, "O", 0.0);
		_octa_atoms_npl[0]=atom("HETATM", 3883, "O", " ", "HOH", "A", 2148, -9.607, 34.284, 15.808, 0.50, 31.80, 0.0, "O", 0.0);
		_octa_atoms_npl[1]=atom("HETATM", 4187, "O", " ", "HOH", "B", 2049, -8.714, 37.322, 13.076, 0.50, 32.75, 0.0, "O", 0.0);


		for (int i=0; i<_no_of_planar_atoms; i++){
			_octa_vertices_pl[i]=Point3D(_octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			_octa_vertices_npl[i]=Point3D(_octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z());
		}

		_planar_nos.push_back(_octa_atoms_pl[0].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[1].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[2].get_serial_no());
		_planar_nos.push_back(_octa_atoms_pl[3].get_serial_no());

		_non_planar_nos.push_back(_octa_atoms_npl[0].get_serial_no());
		_non_planar_nos.push_back(_octa_atoms_npl[1].get_serial_no());

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		cout<<"\nVertices:";
		print(_octa_vertices_pl, _no_of_planar_atoms);
		print(_octa_vertices_npl, _no_of_non_planar_atoms);
	}

	string octaMetalicStructure::assign_standard_site_type(){
		string site_type;
		int aa=0, dna=0, rna=0, other=0, water=0;
		for (int x=0; x<_no_of_planar_atoms; x++){
			if (_octa_atoms_pl[x].get_molecule_type()=="amino_acid"){
				aa++;
			}else if (_octa_atoms_pl[x].get_molecule_type()=="dna_molecule"){
				dna++;
			}else if (_octa_atoms_pl[x].get_molecule_type()=="rna_molecule"){
				rna++;
			}else if (_octa_atoms_pl[x].get_molecule_type()=="other_biom"){
				other++;
			}else if (_octa_atoms_pl[x].get_residue_name()=="HOH"){
				water++;
			}
		}

		for (int x=0; x<_no_of_non_planar_atoms; x++){
			if (_octa_atoms_npl[x].get_molecule_type()=="amino_acid"){
				aa++;
			}else if (_octa_atoms_npl[x].get_molecule_type()=="dna_molecule"){
				dna++;
			}else if (_octa_atoms_npl[x].get_molecule_type()=="rna_molecule"){
				rna++;
			}else if (_octa_atoms_npl[x].get_molecule_type()=="other_biom"){
				other++;
			}else if (_octa_atoms_npl[x].get_residue_name()=="HOH"){
				water++;
			}
		}

		if (water==_no_of_binding_atoms){
			site_type="all_water";
		}else if ((aa!=0) && (dna==0) && (rna==0)){
			site_type="protein";
		}else if ((aa==0) && (dna!=0) && (rna==0)){
			site_type="dna";
		}else if ((aa==0) && (dna==0) && (rna!=0)){
			site_type="rna";
		}else if ((aa!=0) && (dna!=0) && (rna==0)){
			site_type="protein_dna_complex";
		}else if ((aa!=0) && (dna==0) && (rna!=0)){
			site_type="protein_rna_complex";
		}else if ((aa==0) && (dna!=0) && (rna!=0)){
			site_type="dna_rna_complex";
		}else if ((aa!=0) && (dna!=0) && (rna!=0)){
			site_type="protein_dna_rna_complex";
		}else{
			site_type="other";
		}

		return site_type;
	}

