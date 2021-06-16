
#ifndef TETRAGONAL_PYRAMIDAL_H_
#define TETRAGONAL_PYRAMIDAL_H_

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
class tetrapyrMetalicStructure{
	int _no_of_binding_atoms;
	int _no_of_sides;
	int _no_of_planar_atoms;
	int _no_of_non_planar_atom;
	int _no_of_vertex_vertex_angles;
	int _no_of_features_to_compare;	//To calculate mse
	atom _metal_atom;

	int _no_of_coords;
	atom _tetrapyr_atoms_pl[4];	//Four non-metal atoms which are nearest to metal atom that form the tetragonal plane
	atom _tetrapyr_atom_npl;
	vector<long int> _planar_nos;
	vector<long int> _non_planar_no;
	Point3D _tetrapyr_vertices_pl[4];	//Five vertices corresponding to Five non_metal atoms
	Point3D _tetrapyr_vertex_npl;
	tetragonal_pyramidal _tetrapyr;

	tetraMetalicStructure _transformed_tetra_left;
	tetraMetalicStructure _transformed_tetra_right;
	tetrahedron _left_tetra;
	tetrahedron _right_tetra;
	atom _left_tetra_atoms[4];
	atom _right_tetra_atoms[4];

	int _pl_vertex_order[4];
	int _npl_vertex_index;
	int _left_vertex_order[4];
	int _right_vertex_order[4];

	double  _direction_cosines_pl_pl[4][3];	//Direction cosines of edges of standard tetragonal pyramid
	double  _direction_cosines_npl_pl[4][3];	//Direction cosines of edges of standard tetragonal pyramid
	double _vertex_vertex_distances_stand_pl_pl[4][3];
	double _vertex_vertex_distances_stand_npl_pl[4][3];
	double _vertex_vertex_angles_stand_pl_pl[4][4];
	double _vertex_vertex_angles_stand_npl_pl[12][4];
	double _vertex_vertex_distances_comp_pl_pl[4][3];
	double _vertex_vertex_distances_comp_npl_pl[4][3];
	double _vertex_vertex_angles_comp_pl_pl[4][4];
	double _vertex_vertex_angles_comp_npl_pl[12][4];
	double _dist_devt_from_cc_left;
	double _dist_devt_from_cc_right;
	double _devt_dist_vert_vert_left[6];	//Deviation of edges of two tetrahedron
	double _devt_dist_vert_vert_right[6];	//Deviation of edges of two tetrahedron
	double _devt_dist_vert_vert[8];	//Deviation of edges of two tetragonal pyramidal
	double _devt_angle_vertex_vertex_left[12];
	double _devt_angle_vertex_vertex_right[12];
	double _devt_angle_vertex_vertex[16];	//Deviation of angles of two tetragonal pyramidal
	double _mse[2];	//Mean Square Error due to distance deviations and angle deviations
	string _structure_size;

	//void found_planar_non_planar_standard(atom[], vector<int>, vector<int>);
	//bool check_for_coplanarity_standard(vector<Point3D>);
	//double compute_determinant(Point3D, Point3D, Point3D, Point3D);
	//void find_non_planar_atom_indices(atom[], vector<int>, vector<int>);
	void find_largest_side_and_update_vertex_order();
	tetrapyrMetalicStructure make_transformation();

	//void make_planar_base_parallel_xy();
	//void shift_1st_vertex_to_origin();
	//void rotate_largest_side();
	//void rotate_3rd_vertex();
	//void check_and_change_orientation();
	//void reflect_and_shift_2nd_vertex_update_vertex_order(Point3D[], Point3D[], atom);
	//void reflect_4th_vertex_about_xy_plane(Point3D[], Point3D[], atom);
	//tetrapyrMetalicStructure shift_circum_centre_to_origin();	//Shift circum centre to origin
	//bool check_distance(Point3D*, int*);
	void compute_mse();
	void gen_report_vertex_vertex_distances(FILE*, atom*, double[][3], int);
	void gen_report_vertex_vertex_distances(FILE*, atom*, atom, double[][3]);
	void gen_report_vertex_vertex_angles(FILE*, atom*, double[][4], int);
	void gen_report_vertex_vertex_angles(FILE*, atom*, atom, double[][4]);
	void write_metal(FILE*);
	void write_residues(FILE*);
	void print(atom, atom*, int);
	void print(Point3D*, int);
	void print(Point3D, int);
public:
	tetrapyrMetalicStructure();
	tetrapyrMetalicStructure(tetrahedron left, atom left_atoms[], tetrahedron right, atom right_atoms[], atom metal_atom, vector<long int> platom, vector<long int> nonplatom);
	tetrapyrMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom, vector<long int> platom, vector<long int> nonplatom);
	tetrapyrMetalicStructure(tetraMetalicStructure left_transformed, tetraMetalicStructure right_transformed, atom metal_atom);
	void form_standard_tetragonal_pyramid_NA();
	void form_standard_tetragonal_pyramid_CA();
	void form_standard_tetragonal_pyramid_MG();
	void form_standard_tetragonal_pyramid_CD();
	void form_standard_tetragonal_pyramid_NI();
	void form_standard_tetragonal_pyramid_ZN();
	string assign_standard_site_type();

	tetrapyrMetalicStructure apply_transformations();
	atom get_metal();
	atom* get_tetrapyr_atoms_pl();
	atom get_tetrapyr_atom_npl();
	Point3D* get_tetrapyr_vertices_pl();
	Point3D get_tetrapyr_vertex_npl();
	tetragonal_pyramidal get_tetrapyr();
	tetrahedron get_left_tetra();
	atom* get_left_tetra_atoms();
	tetrahedron get_right_tetra();
	atom* get_right_tetra_atoms();
	int get_no_of_binding_atoms();
	double* get_side_deviations();
	double* get_angle_deviations_vertex_vertex();
	double* get_mse();
	string get_structure_size();

	void comparison(tetrapyrMetalicStructure &, string);
	void write_atoms_pdb_format(string fname);
	void write_atoms_combined_tex_format_left(tetrapyrMetalicStructure & ctdrn, string fname);
	void write_atoms_combined_tex_format_right(tetrapyrMetalicStructure & ctdrn, string fname);
	void measure_direction_cosines();
	void compare(tetrapyrMetalicStructure &);
	void gen_report_metal_binding_sites(tetrapyrMetalicStructure &, tetrapyrMetalicStructure &, tetrapyrMetalicStructure &, string, string, string, FILE*);
	//void free_memory();
};

#endif /* TETRAGONAL_PYRAMIDAL_H_ */
