
#ifndef TRIGONAL_BIPYRAMIDAL_H_
#define TRIGONAL_BIPYRAMIDAL_H_

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
class tribipyrMetalicStructure{
	int _no_of_binding_atoms;
	int _no_of_sides;
	int _no_of_planar_atoms;
	int _no_of_non_planar_atoms;
	int _no_of_vertex_vertex_angles;
	int _no_of_features_to_compare;	//To calculate mse
	atom _metal_atom;

	int _no_of_coords;
	atom _tribipyr_atoms_pl[3];	//Five non-metal atoms which are nearest to metal atom that form the triangle
	atom _tribipyr_atoms_npl[2];
	vector<long int> _planar_nos;
	vector<long int> _non_planar_nos;
	Point3D _tribipyr_vertices_pl[3];	//Five vertices corresponding to Five non_metal atoms
	Point3D _tribipyr_vertices_npl[2];
	trigonal_bipyramidal _tribipyr;
	tetrahedron _upper_tetra;
	tetrahedron _lower_tetra;
	atom _upper_tetra_atoms[4];
	atom _lower_tetra_atoms[4];

	int _pl_vertex_order[3];
	int _npl_vertex_order[2];

	double  _direction_cosines_pl[3][3];	//Direction cosines of edges of standard trigonal bipyramid
	double  _direction_cosines_npl[6][3];	//Direction cosines of edges of standard trigonal bipyramid
	double _vertex_vertex_distances_stand_pl_pl[3][3];
	double _vertex_vertex_distances_stand_npl_pl[6][3];
	double _vertex_vertex_angles_stand_pl_pl[3][4];
	double _vertex_vertex_angles_stand_npl_pl[18][4];
	double _vertex_vertex_distances_comp_pl_pl[3][3];
	double _vertex_vertex_distances_comp_npl_pl[6][3];
	double _vertex_vertex_angles_comp_pl_pl[3][4];
	double _vertex_vertex_angles_comp_npl_pl[18][4];
	double _dist_devt_from_cc_upper;
	double _dist_devt_from_cc_lower;
	double _devt_dist_vert_vert_upper[6];	//Deviation of edges of two tetrahedron
	double _devt_dist_vert_vert_lower[6];	//Deviation of edges of two tetrahedron
	double _devt_dist_vert_vert[9];	//Deviation of edges of two trigonal bipyramidal
	double _devt_angle_vertex_vertex_upper[12];
	double _devt_angle_vertex_vertex_lower[12];
	double _devt_angle_vertex_vertex[21];	//Deviation of angles of two trigonal bipyramidal
	double _mse[2];	//Mean Square Error due to distance deviations and angle deviations
	string _structure_size;

	//void found_planar_non_planar_standard(atom[], vector<int>, vector<int>);
	//bool check_for_coplanarity_standard(vector<Point3D>);
	//double compute_determinant(Point3D, Point3D, Point3D, Point3D);
	//void find_non_planar_atom_indices(atom[], vector<int>, vector<int>);
	void find_largest_side_and_update_vertex_order();
	tribipyrMetalicStructure make_transformation();
	//void make_planar_base_parallel_xy();
	//void shift_1st_vertex_to_origin();
	//void rotate_largest_side();
	//void rotate_3rd_vertex();
	//void check_and_change_orientation();
	//void reflect_and_shift_2nd_vertex_update_vertex_order(Point3D[], Point3D[], atom);
	//void reflect_4th_vertex_about_xy_plane(Point3D[], Point3D[], atom);
	//tribipyrMetalicStructure shift_circum_centre_to_origin();	//Shift circum centre to origin
	//bool check_distance(Point3D*, int*);
	void compute_mse();
	void gen_report_vertex_vertex_distances(FILE*, atom[], double[][3], int);
	void gen_report_vertex_vertex_distances(FILE*, atom[], atom*, double[][3]);
	void gen_report_vertex_vertex_angles(FILE*, atom[], double[][4], int);
	void gen_report_vertex_vertex_angles(FILE*, atom[], atom*, double[][4]);
	void write_metal(FILE*);
	void write_residues(FILE*);
	void print(atom, atom[], int);
	void print(Point3D[], int);
public:
	tribipyrMetalicStructure();
	tribipyrMetalicStructure(tetrahedron upper, atom upper_atoms[], tetrahedron lower, atom lower_atoms[], atom metal_atom);
	tribipyrMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom, vector<long int> platom, vector<long int> nonplatom);
	//tribipyrMetalicStructure(tetraMetalicStructure upper_transformed, tetraMetalicStructure lower_transformed, atom metal_atom);
	void form_standard_trigonal_bipyramid_NA();
	void form_standard_trigonal_bipyramid_CA();
	void form_standard_trigonal_bipyramid_MG();
	void form_standard_trigonal_bipyramid_MN();
	void form_standard_trigonal_bipyramid_FE();
	void form_standard_trigonal_bipyramid_ZN();
	void form_standard_trigonal_bipyramid_CD();
	void form_standard_trigonal_bipyramid_CO();
	void form_standard_trigonal_bipyramid_NI();
	void form_standard_trigonal_bipyramid_HG();
	string assign_standard_site_type();

	tribipyrMetalicStructure apply_transformations();
	atom get_metal();
	atom* get_tribipyr_atoms_pl();
	atom* get_tribipyr_atoms_npl();
	Point3D* get_tribipyr_vertices_pl();
	Point3D* get_tribipyr_vertices_npl();
	trigonal_bipyramidal get_tribipyr();
	tetrahedron get_upper_tetra();
	atom* get_upper_tetra_atoms();
	tetrahedron get_lower_tetra();
	atom* get_lower_tetra_atoms();
	int get_no_of_binding_atoms();
	double* get_side_deviations();
	double* get_angle_deviations_vertex_vertex();
	double* get_mse();
	string get_structure_size();

	void comparison(tribipyrMetalicStructure &, string);
	void write_atoms_pdb_format(string fname);
	void write_atoms_combined_tex_format_upper(tribipyrMetalicStructure & ctdrn, string fname);
	void write_atoms_combined_tex_format_lower(tribipyrMetalicStructure & ctdrn, string fname);
	void measure_direction_cosines();
	void compare(tribipyrMetalicStructure &);
	void gen_report_metal_binding_sites(tribipyrMetalicStructure &, tribipyrMetalicStructure &, tribipyrMetalicStructure &, string, string, string, FILE*);
	//void free_memory();
};

#endif /* TRIGONAL_BIPYRAMIDAL_H_ */
