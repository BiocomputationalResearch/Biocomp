
#ifndef OCTAHEDRAL_H_
#define OCTAHEDRAL_H_

#include <iostream>
#include <cmath>
#include <string>
#include <assert.h>
#include <bits/stdc++.h>
#include <sstream>
#include <iomanip>
#include "submolecule.h"
#include "ligand.h"
#include "tetrahedral.h"

using namespace geom;
#define PI 3.141592654
class octaMetalicStructure{
	int _no_of_binding_atoms;
	int _no_of_sides;
	int _no_of_planar_atoms;
	int _no_of_non_planar_atoms;
	int _no_of_vertex_vertex_angles;
	int _no_of_features_to_compare;	//To calculate mse
	atom _metal_atom;

	int _no_of_coords;
	atom _octa_atoms_pl[4];	//Four non-metal atoms which are nearest to metal atom that form the tetragonal plane
	atom _octa_atoms_npl[2];
	vector<long int> _planar_nos;
	vector<long int> _non_planar_nos;
	Point3D _octa_vertices_pl[4];	//Five vertices corresponding to Five non_metal atoms
	Point3D _octa_vertices_npl[2];
	octahedral _octa;

	tetraMetalicStructure _transformed_tetra_upper_left;
	tetraMetalicStructure _transformed_tetra_upper_right;
	tetraMetalicStructure _transformed_tetra_lower_left;
	tetraMetalicStructure _transformed_tetra_lower_right;
	tetrahedron _upper_left_tetra;
	tetrahedron _upper_right_tetra;
	tetrahedron _lower_left_tetra;
	tetrahedron _lower_right_tetra;
	atom _upper_left_tetra_atoms[4];
	atom _upper_right_tetra_atoms[4];
	atom _lower_left_tetra_atoms[4];
	atom _lower_right_tetra_atoms[4];

	int _pl_vertex_order[4];
	int _npl_vertex_order[2];
	int _upper_left_vertex_order[4];
	int _upper_right_vertex_order[4];
	int _lower_left_vertex_order[4];
	int _lower_right_vertex_order[4];

	double  _direction_cosines_pl_pl[4][3];	//Direction cosines of edges of standard tetragonal pyramid
	double  _direction_cosines_npl_pl[8][3];	//Direction cosines of edges of standard tetragonal pyramid
	double _vertex_vertex_distances_stand_pl_pl[4][3];
	double _vertex_vertex_distances_stand_npl_pl[8][3];
	double _vertex_vertex_angles_stand_pl_pl[4][4];
	double _vertex_vertex_angles_stand_npl_pl[24][4];
	double _vertex_vertex_distances_comp_pl_pl[4][3];
	double _vertex_vertex_distances_comp_npl_pl[8][3];
	double _vertex_vertex_angles_comp_pl_pl[4][4];
	double _vertex_vertex_angles_comp_npl_pl[24][4];
	double _dist_devt_from_cc_upper_left;
	double _dist_devt_from_cc_upper_right;
	double _dist_devt_from_cc_lower_left;
	double _dist_devt_from_cc_lower_right;
	double _devt_dist_vert_vert_upper_left[6];	//Deviation of edges of two tetrahedron
	double _devt_dist_vert_vert_upper_right[6];	//Deviation of edges of two tetrahedron
	double _devt_dist_vert_vert_lower_left[6];	//Deviation of edges of two tetrahedron
	double _devt_dist_vert_vert_lower_right[6];	//Deviation of edges of two tetrahedron
	double _devt_dist_vert_vert[12];	//Deviation of edges of two tetragonal pyramidal
	double _devt_angle_vertex_vertex_upper_left[12];
	double _devt_angle_vertex_vertex_upper_right[12];
	double _devt_angle_vertex_vertex_lower_left[12];
	double _devt_angle_vertex_vertex_lower_right[12];
	double _devt_angle_vertex_vertex[28];	//Deviation of angles of two tetragonal pyramidal
	double _mse[2];	//Mean Square Error due to distance deviations and angle deviations
	string _structure_size;

	void find_largest_side_and_update_vertex_order();
	octaMetalicStructure make_transformation();

	void compute_mse();
	void gen_report_vertex_vertex_distances(FILE*, atom*, double[][3], int);
	void gen_report_vertex_vertex_distances(FILE*, atom*, atom*, int*, int*, double[][3]);
	void gen_report_vertex_vertex_angles(FILE*, atom*, double[][4], int);
	void gen_report_vertex_vertex_angles(FILE*, atom*, atom*, int*, int*, double[][4]);
	void write_metal(FILE*);
	void write_residues(FILE*);
	void print(atom, atom*, int);
	void print(Point3D*, int);
	void print(Point3D, int);
public:
	octaMetalicStructure();
	octaMetalicStructure(tetrahedron upper_left, atom upper_left_atoms[], tetrahedron upper_right, atom upper_right_atoms[], tetrahedron lower_left, atom lower_left_atoms[], tetrahedron lower_right, atom lower_right_atoms[], atom metal_atom, vector<long int> platom, vector<long int> nonplatom);
	octaMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom, vector<long int> platom, vector<long int> nonplatom);
	octaMetalicStructure(tetraMetalicStructure upper_left_transformed, tetraMetalicStructure upper_right_transformed, tetraMetalicStructure lower_left_transformed, tetraMetalicStructure lower_right_transformed, atom metal_atom);
	void form_standard_octa_NA();
	void form_standard_octa_CA();
	void form_standard_octa_MG();
	void form_standard_octa_MN();
	void form_standard_octa_FE();
	void form_standard_octa_ZN();
	void form_standard_octa_CD();
	void form_standard_octa_CO();
	void form_standard_octa_NI();
	string assign_standard_site_type();

	octaMetalicStructure apply_transformations();
	atom get_metal();
	atom* get_octa_atoms_pl();
	atom* get_octa_atoms_npl();
	Point3D* get_octa_vertices_pl();
	Point3D* get_octa_vertices_npl();
	int* get_pl_vertex_order();
	int* get_npl_vertex_order();
	octahedral get_octa();
	tetrahedron get_upper_left_tetra();
	atom* get_upper_left_tetra_atoms();
	tetrahedron get_upper_right_tetra();
	atom* get_upper_right_tetra_atoms();
	tetrahedron get_lower_left_tetra();
	atom* get_lower_left_tetra_atoms();
	tetrahedron get_lower_right_tetra();
	atom* get_lower_right_tetra_atoms();
	int get_no_of_binding_atoms();
	double* get_side_deviations();
	double* get_angle_deviations_vertex_vertex();
	double* get_mse();
	string get_structure_size();

	void comparison(octaMetalicStructure &, string);
	void write_atoms_pdb_format(string fname);
	void write_atoms_combined_tex_format_upper_left(octaMetalicStructure & ctdrn, string fname);
	void write_atoms_combined_tex_format_upper_right(octaMetalicStructure & ctdrn, string fname);
	void write_atoms_combined_tex_format_lower_left(octaMetalicStructure & ctdrn, string fname);
	void write_atoms_combined_tex_format_lower_right(octaMetalicStructure & ctdrn, string fname);
	void measure_direction_cosines();
	void compare(octaMetalicStructure &);
	void gen_report_metal_binding_sites(octaMetalicStructure &, octaMetalicStructure &, octaMetalicStructure &, string, string, string, FILE*);
	//void free_memory();
};

#endif /* OCTAHEDRAL_H_ */
