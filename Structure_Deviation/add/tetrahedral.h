
#ifndef TETRAHEDRAL_H_
#define TETRAHEDRAL_H_

#include <iostream>
#include <cmath>
#include <string>
#include <assert.h>
#include <bits/stdc++.h>
#include <sstream>
#include <iomanip>
#include "submolecule.h"
#include "ligand.h"

using namespace geom;
#define PI 3.141592654
class tetraMetalicStructure{
	int _no_of_binding_atoms;
	int _no_of_sides;
	int _no_of_angles_with_circum_centre;
	int _no_of_vertex_vertex_angles;
	int _no_of_features_to_compare;	//To calculate mse
	atom _metal_atom;

	int _no_of_coords;
	atom _tetra_atoms[4];	//Four non-metal atoms which are nearest to metal atom that form the tetrahedron
	string _site_type_comparable;
	Point3D _tetra_vertices[4];	//Four vertices corresponding to Four non_metal atoms
	tetrahedron _tetra;
	atom _parallel_tetra_atoms[4];
	atom _parallel_metal_atom;
	atom _oriented_metal_atom;
	atom _oriented_tetra_atoms[4];
	atom _cc_shifted_tetra_atoms[4];	//Atoms of circumcentre shifted to origin
	atom _cc_shifted_metal_atom;

	Point3D _shifted_vertices[4];	//Shifted 1st vertex
	atom _shifted_metal;
	Point3D _rotated_vertices_largest_side[4];
	atom _rotated_metal_zaxis;
	atom _rotated_metal_yaxis;
	Point3D _parallel_vertices[4];
	int _vertex_order[4];
	Point3D _reflected_vertices_3rd_vertex[4];
	Point3D _reflected_vertices_4th_vertex[4];

	Point3D _oriented_vertices[4];
	Point3D _cc_shifted_vertices[4];
	//double _rotate_angle[2];

	double  _direction_cosines[6][3];	//Direction cosines of edges of Standard Tetrahedron
	double _cc_vertex_distances_stand[4][2];
	double _vertex_vertex_distances_stand[6][3];
	double _cc_vertex_angles_stand[6][3];
	double _vertex_vertex_angles_stand[12][4];
	double _cc_vertex_distances_comp[4][2];
	double _vertex_vertex_distances_comp[6][3];
	double _cc_vertex_angles_comp[6][3];
	double _vertex_vertex_angles_comp[12][4];
	double _dist_devt_from_cc;
	double _devt_dist_vert_vert[6];	//Deviation of edges of two tetrahedrons
	double _devt_angle_with_cc[6];	//Deviation of angles of two tetrahedrons
	double _devt_angle_vertex_vertex[12];
	double _mse[4];	//Mean Square Error due to two types of distance deviations and angle deviations
	bool _nearer;
	string _structure_size;

	void make_base_parallel_xy();
	void shift_1st_vertex_to_origin();
	void rotate_largest_side();
	void rotate_3rd_vertex();
	void check_and_change_orientation();
	void change_orientation();
	void reflect_and_shift_2nd_vertex_update_vertex_order(Point3D [], atom);
	atom reflect_4th_vertex_about_xy_plane(Point3D [], atom);
	tetraMetalicStructure shift_circum_centre_to_origin();	//Shift circum centre to origin
	tetraMetalicStructure shift_circum_centre_to_origin(int []);	//Shift circum centre to origin
	bool check_distance(Point3D[], int[]);
	void compute_mse();
	void gen_report_circumcentre_vertex_distances(FILE*, atom [], double [][2]);
	void gen_report_vertex_vertex_distances(FILE*, atom [], double [][3]);
	void gen_report_circumcentre_vertex_angles(FILE*, atom [], double [][3]);
	void gen_report_vertex_vertex_angles(FILE*, atom [], double [][4]);
	void write_metal(FILE*);
	void write_residues(FILE*);
	void print(atom, atom []);
	void print(Point3D []);
public:
	tetraMetalicStructure();
	tetraMetalicStructure(atom new_atm[], atom metal_atom);
	tetraMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom);
	tetraMetalicStructure(atom new_atm[], Point3D new_vert[], atom metal_atom, int v1, int v2, int v3, int v4);
	void form_standard_tetrahedron_NA();
	void form_standard_tetrahedron_CA();
	void form_standard_tetrahedron_K();
	void form_standard_tetrahedron_MG();
	void form_standard_tetrahedron_MN();
	void form_standard_tetrahedron_FE();
	void form_standard_tetrahedron_CU();
	void form_standard_tetrahedron_ZN();
	void form_standard_tetrahedron_CO();
	void form_standard_tetrahedron_MO();
	void form_standard_tetrahedron_AG();
	void form_standard_tetrahedron_PB();
	void form_standard_tetrahedron_HG();
	void form_standard_tetrahedron_AS();
	void form_standard_tetrahedron_CD();
	void form_standard_tetrahedron_NI();
	void form_standard_tetrahedron_CR();
	void form_standard_tetrahedron_CS();
	string assign_standard_site_type();

	tetraMetalicStructure apply_transformations();
	tetraMetalicStructure transform_tetrahedron();
	atom get_metal();
	atom* get_tetra_atoms();
	tetrahedron get_tetrahedron();
	int get_no_of_binding_atoms();
	double* get_side_deviations();
	double get_distance_deviations_from_cc();
	double* get_angle_deviations_with_cc();
	double* get_angle_deviations_vertex_vertex();
	double* get_mse();
	bool get_nearer();
	string get_structure_size();

	void comparison(tetraMetalicStructure &, string);
	void write_atoms_pdb_format(string fname);
	void write_atoms_single_tex_format(string fname);
	void write_atoms_combined_tex_format(tetraMetalicStructure & ctdrn, string fname);
	void measure_direction_cosines();
	void compare(tetraMetalicStructure &);
	void gen_report_metal_binding_sites(tetraMetalicStructure &, tetraMetalicStructure &, tetraMetalicStructure &, string, string, string, FILE*);
	void gen_pdb_model_file(tetraMetalicStructure &, tetraMetalicStructure &, tetraMetalicStructure &, FILE*);
};


#endif /* TETRAHEDRAL_H_ */
