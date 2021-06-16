
#ifndef TRIGONALPLANAR_H_
#define TRIGONALPLANAR_H_

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
class triplMetalicStructure{
	int _no_of_binding_atoms;
	int _no_of_sides;
	int _no_of_angles_with_circum_centre;
	int _no_of_vertex_vertex_angles;
	int _no_of_features_to_compare;	//To calculate mse
	atom _metal_atom;

	int _no_of_coords;
	atom _tripl_atoms[3];	//Three non-metal atoms which are nearest to metal atom that form the triangle
	string _site_type_comparable;
	Point3D _tripl_vertices[3];	//Three vertices corresponding to Three non_metal atoms
	triangle _tripl;
	atom _parallel_tripl_atoms[3];
	atom _parallel_metal_atom;
	atom _oriented_metal_atom;
	atom _oriented_tripl_atoms[3];
	atom _cc_shifted_tripl_atoms[3];	//Atoms of circumcentre shifted to origin
	atom _cc_shifted_metal_atom;

	Point3D _shifted_vertices[3];	//Shifted 1st vertex
	atom _shifted_metal;
	Point3D _rotated_vertices_largest_side[3];
	atom _rotated_metal_zaxis;
	atom _rotated_metal_yaxis;
	Point3D _parallel_vertices[3];
	int _vertex_order[3];
	Point3D _reflected_vertices_3rd_vertex[3];

	Point3D _oriented_vertices[3];
	Point3D _cc_shifted_vertices[3];

	double  _direction_cosines[][3];	//Direction cosines of edges of Standard Triangle
	double _cc_vertex_distances_stand[3][2];
	double _vertex_vertex_distances_stand[3][3];
	double _cc_vertex_angles_stand[3][3];
	double _vertex_vertex_angles_stand[3][4];
	double _cc_vertex_distances_comp[3][2];
	double _vertex_vertex_distances_comp[3][3];
	double _cc_vertex_angles_comp[3][3];
	double _vertex_vertex_angles_comp[3][4];
	double _dist_devt_from_cc;
	double _devt_dist_vert_vert[3];	//Deviation of edges of two triangles
	double _devt_angle_with_cc[3];	//Deviation of angles of two triangles
	double _devt_angle_vertex_vertex[3];
	double _mse[4];	//Mean Square Error due to distance deviations and angle deviations
	bool _nearer;
	string _structure_size;

	void make_triangle_parallel_xy();
	void shift_1st_vertex_to_origin();
	void rotate_largest_side();
	void rotate_3rd_vertex();
	void check_and_change_orientation();
	void change_orientation();
	void reflect_and_shift_2nd_vertex_update_vertex_order(Point3D*, atom);
	triplMetalicStructure shift_circum_centre_to_origin();	//Shift circum centre to origin
	bool check_distance(Point3D[]);
	void compute_mse();
	void gen_report_circumcentre_vertex_distances(FILE*, atom[], double[][2]);
	void gen_report_vertex_vertex_distances(FILE*, atom[], double[][3]);
	void gen_report_circumcentre_vertex_angles(FILE*, atom[], double[][3]);
	void gen_report_vertex_vertex_angles(FILE*, atom[], double[][4]);
	void write_metal(FILE*);
	void write_residues(FILE*);
	void print(atom, atom*);
	void print(Point3D*);
public:
	triplMetalicStructure();
	triplMetalicStructure(atom* new_atm, atom metal_atom);
	triplMetalicStructure(atom* new_atm, atom metal_atom, int v1, int v2, int v3);
	triplMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom);
	triplMetalicStructure(atom tri_atom[], Point3D* tri_vert, atom metal, int v1, int v2, int v3);
	void form_standard_triangle_NA();
	void form_standard_triangle_CA();
	void form_standard_triangle_K();
	void form_standard_triangle_MG();
	void form_standard_triangle_MN();
	void form_standard_triangle_FE();
	void form_standard_triangle_CU();
	void form_standard_triangle_ZN();
	void form_standard_triangle_CO();
	void form_standard_triangle_PB();
	void form_standard_triangle_HG();
	void form_standard_triangle_AS();
	void form_standard_triangle_CD();
	void form_standard_triangle_NI();
	void form_standard_triangle_CS();
	string assign_standard_site_type();

	triplMetalicStructure apply_transformations();
	triplMetalicStructure transform_triangle();
	atom get_metal();
	atom* get_tripl_atoms();
	atom* get_cc_shifted_tripl_atoms();
	Point3D* get_cc_shifted_tripl_vertices();
	triangle get_triangle();
	int get_no_of_binding_atoms();
	double* get_side_deviations();
	double get_distance_deviations_from_cc();
	double* get_angle_deviations_with_cc();
	double* get_angle_deviations_vertex_vertex();
	double* get_mse();
	bool get_nearer();
	bool get_updated_vertex_order_2nd_vertex();
	bool get_updated_vertex_order_stand();
	string get_structure_size();

	void comparison(triplMetalicStructure &, string);
	void write_atoms_pdb_format(string fname);
	void write_atoms_single_tex_format(string fname);
	void write_atoms_combined_tex_format(triplMetalicStructure & ctdrn, string fname);
	void measure_direction_cosines();
	void compare(triplMetalicStructure &);
	void gen_report_metal_binding_sites(triplMetalicStructure &, triplMetalicStructure &, triplMetalicStructure &, string, string, string, FILE*);
};

#endif /* TRIGONALPLANAR_H_ */
