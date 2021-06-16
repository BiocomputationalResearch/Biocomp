
#ifndef TETRAGONALPLANAR_H_
#define TETRAGONALPLANAR_H_

#include <iostream>
#include <cmath>
#include <string>
#include <assert.h>
#include <bits/stdc++.h>
#include <sstream>
#include <iomanip>
#include "submolecule.h"
#include "ligand.h"
#include "trigonalplanar.h"

using namespace geom;
#define PI 3.141592654
class tetraplMetalicStructure{
	int _no_of_binding_atoms;
	int _no_of_sides;
	int _no_of_vertex_vertex_angles;
	int _no_of_features_to_compare;	//To calculate mse
	int _no_of_coords;
	int _no_of_triangles;
	atom _metal_atom;
	tetragonalplanar _tetrapl;
	string _site_type_comparable;
	atom _tetrapl_atoms[4];
	Point3D _tetrapl_vertices[4];	//Four vertices corresponding to Four non_metal atoms

	atom _cc_shifted_left_triangle_atoms[3];	//Atoms of circumcentre shifted to origin Left Triangle
	atom _cc_shifted_right_triangle_atoms[3];	//Atoms of circumcentre shifted to origin Right Triangle

	int _left_vertex_order[3];
	int _right_vertex_order[3];
	int _vertex_order[4];

	Point3D _cc_shifted_vertices_left[3];
	Point3D _cc_shifted_vertices_right[3];

	double  _direction_cosines[4][3];	//Direction cosines of edges of Standard Triangle
	double _vertex_vertex_distances_stand[4][3];
	double _vertex_vertex_angles_stand[4][4];
	double _vertex_vertex_distances_comp[4][3];
	double _vertex_vertex_angles_comp[4][4];
	double _dist_devt_from_cc_left;
	double _dist_devt_from_cc_right;
	double _devt_dist_vert_vert_left[3];	//Deviation of edges of two triangles
	double _devt_dist_vert_vert_right[3];	//Deviation of edges of two triangles
	double _devt_angle_vertex_vertex_left[3];
	double _devt_angle_vertex_vertex_right[3];
	double _devt_dist_vert_vert[4];	//Deviation of edges of two trigonal bipyramidal
	double _devt_angle_vertex_vertex[4];	//Deviation of angles of two trigonal bipyramidal

	triangle _left_tri;
	triangle _right_tri;
	triplMetalicStructure _transformed_tri_left;
	triplMetalicStructure _transformed_tri_right;
	double _mse[2];		//Mean Square Error due to distance deviations and angle deviations
	string _structure_size;

	void find_largest_side_and_update_vertex_order();
	tetraplMetalicStructure make_transformation();
	void compute_mse();
	void gen_report_vertex_vertex_distances(FILE*, atom*, double[][3]);
	void gen_report_vertex_vertex_angles(FILE*, atom*, double[][4]);
	void write_metal(FILE*);
	void write_residues(FILE*);
	void print(atom, atom*);
	void print(Point3D*);
public:
	tetraplMetalicStructure();
	tetraplMetalicStructure(atom* new_atm, atom metal_atom);
	tetraplMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom);
	tetraplMetalicStructure(triangle tri1, atom atom_tri1[], triangle tri2, atom atom_tri2[], atom metal);
	tetraplMetalicStructure(triplMetalicStructure tps_left, triplMetalicStructure tps_right, atom metal);
	void form_standard_tetragonal_planar_MG();
	void form_standard_tetragonal_planar_FE();
	void form_standard_tetragonal_planar_CD();
	void form_standard_tetragonal_planar_NI();
	string assign_standard_site_type();

	tetraplMetalicStructure apply_transformations();
	Point3D* get_tetrapl_vertices();
	atom get_metal();
	atom* get_tetrapl_atoms();
	triangle get_left_triangle();
	triangle get_right_triangle();
	atom* get_cc_shifted_left_triangle_atoms();
	atom* get_cc_shifted_right_triangle_atoms();

	tetragonalplanar get_tetragonalplanar();
	int* get_vertex_order();
	int get_no_of_binding_atoms();
	double* get_side_deviations();
	double get_distance_deviations_from_cc();
	double* get_angle_deviations_with_cc();
	double* get_angle_deviations_vertex_vertex();
	double* get_mse();
	bool get_nearer();
	string get_structure_size();
	void compare_with_standard(tetraplMetalicStructure &);
	void comparison(tetraplMetalicStructure &, string);
	void write_atoms_pdb_format(string fname);
	void write_atoms_combined_tex_format_left(tetraplMetalicStructure & ctdrn, string fname);
	void write_atoms_combined_tex_format_right(tetraplMetalicStructure & ctdrn, string fname);
	void measure_direction_cosines();
	void gen_report_metal_binding_sites(tetraplMetalicStructure &, tetraplMetalicStructure &, tetraplMetalicStructure &, string, string, string, FILE*);
};

#endif /* TETRAGONALPLANAR_H_ */
