#ifndef COORDINATIONSITE_H_
#define COORDINATIONSITE_H_

#include <iostream>
#include <string>
#include <cmath>

#include "submolecule.h"
#include "ligand.h"

class metal_binding_site{
protected:
	string _accn;
	atom _metal_atom;
	vector<atom> _non_metals;
	vector<atom> _nearest_atoms;
	string _site_type;
	vector<direct_ligand> _inner_coord_sphere;	//Excluding Atoms for Alternate Conformation
	double _mean_distance_direct_ligands_from_metal;
	double _mean_inter_atomic_distance;
	double _mean_angle_atom_metal_atom;
	vector<string> _metal_list;	//Unique Metals binded to protein
	vector<bond> _bond_info;
	vector<atom> _water_group;
	int _inner_coord_no;
	string _shape_of_inner_sphere;
	vector<long int> _planar_atoms;
	vector<long int> _non_planar_atoms;
	double _mean_planar_base_atoms_angles_with_metal;
	double _mean_non_planar_atoms_angles_with_metal;
	double _mean_planar_and_non_planar_atoms_angles_with_metal;
	int _no_of_base_atoms;
	double det_threshold;
	double dist_threshold;

	void compute_nearest_atom();
	void determine_shape();
	void make_order_direct_atoms();
	void separate_atoms();
	void assign_site_type();
	void delete_ligand_tree();
	//string form_children_string(ligand *, vector<long int>);
	vector<int> count_planar_atoms();
	void find_non_planar_atom_indices(vector<int>, int*);
	bool is_same_side(int*, vector<int>);
	bool check_side(vector<Point3D>, vector<Point3D>);
	double compute_determinant(Point3D, Point3D, Point3D);
	double compute_determinant(Point3D, Point3D, Point3D, Point3D);
	bool check_for_coplanarity(vector<Point3D>, vector<double> &);
	bool check_for_coplanarity(vector<Point3D>);
	//long int find_combination(int, int);
	//long int factorial(int);
	void store_combination(vector<int>, int, int, vector<vector<int> > &);
	void all_combinations(vector<int>, int[], int, int, int, int, vector<vector<int> > &);
	bool is_square(Point3D, Point3D, Point3D, Point3D);
	bool is_cube(int*, vector<int>);
	bool is_square_prism(int*, vector<int>);
	bool check_angles_and_distances(Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, Point3D);
	bool check_equilateral_faces(Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, Point3D);
	bool check_tri_angles(Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, Point3D);
	bool check_for_trigonal_prismatic(int*, vector<int>);
	bool check_for_trigonal_bipyramidal(int*, vector<int>);
	bool check_for_hemispheric(vector<Point3D>, vector<Point3D>);
	bool check_for_capped_trigonal_prismatic(int, vector<int>);
	bool check_for_octahedral(int*, vector<int>);
	bool check_for_capped_octahedral(int, int*, vector<int>);
	bool check_side(Point3D, Point3D, Point3D, Point3D);
	bool test_for_perfect_seesaw();
	bool test_for_seesaw();
	bool check_for_capped_trigonal_prismatic(int*, vector<int>);
	bool is_parallel(Point3D, Point3D, Point3D, Point3D);
	bool check_position_for_tetragonal_face_bicapped(vector<int>, vector<int>, int*);
	bool check_cap_sides(Point3D, Point3D, Point3D, Point3D, Point3D, Point3D, vector<int>);
	bool check_position_for_triangular_face_bicapped(vector<int>, vector<int>, int*);
	bool check_intersection_with_plane(Point3D, Point3D, Point3D, Point3D, Point3D);
	bool check_position_for_square_antiprism(int*, vector<int>);
	bool check_position_for_tricapped(vector<int>, vector<int>, int*);
	bool check_for_dodecahedron(int*, vector<int>);
	bool check_for_capped_square_anti_prism(int*, vector<int>, int);
	void compute_mean_angle_planar_atoms_with_metal(vector<int>);
	void compute_mean_angle_non_planar_atoms_with_metal(int *, int);
	void compute_mean_angle_planar_and_non_planar_atoms_with_metal(int *, int, vector<int>);
	void compute_inter_atomic_distances();
	void compute_atom_to_metal_distances();
	void compute_atom_metal_atom_angles();
	void compute_atom_atom_atom_angles();

public:
	metal_binding_site(){}
	metal_binding_site(string, atom, vector<string>, vector<atom>, vector<bond>);
	vector<atom> get_nearest_atoms();
	string get_site_type();
	vector<direct_ligand> get_inner_coord_sphere();
	vector<atom> get_water_molecules();
	int get_inner_coord_no();
	string get_shape();
	vector<long int> get_planar_atoms();
	vector<long int> get_non_planar_atoms();
	atom & get_metal_atom();
	void print();
};

#endif /* COORDINATIONSITE_H_ */
