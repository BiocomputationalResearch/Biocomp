#ifndef TETRA_STRUCTURE_CA_H_
#define TETRA_STRUCTURE_CA_H_

#include <assert.h>
#include "coordinationsite.h"
#include "tetrahedral.h"
#include "stdio.h"

class tetraStructCa{
	string _accn;
	metal_binding_site _tetra_site;
	static int _no_of_binding_atoms;
	static int _no_of_sides;
	static int _no_of_angles_with_circum_centre;
	static int _no_of_vertex_vertex_angles;
	static int _no_of_features_to_compare;
	static tetraMetalicStructure* _tetra_struct_standard_ca;
	tetraMetalicStructure* _tetra_struct_comparable_ca;
	string _site_type_comparable;
	static double _avg_dist_devt_from_cc;
	static double _avg_devt_dist_vert_vert[];	//Average deviation of edges of two tetrahedrons
	static double _avg_devt_angle_with_cc[];	//Average deviation of angles of two tetrahedrons
	static double _avg_devt_angle_vertex_vertex[];
	static double _avg_mse[];	//Average Mean Square Error due to distance deviations and angle deviations
	static int _nearer;
	static vector<string> _nearer_mol_type;
	static vector<string> nearer_protein_name;
	static int _far_apart;
	static vector<string> _far_mol_type;
	static vector<string> far_apart_protein_name;
	static int _larger_size;
	static vector<string> _larger_mol_type;
	static vector<string> larger_protein_name;
	static int _smaller_size;
	static vector<string> _smaller_mol_type;
	static vector<string> smaller_protein_name;
	static int _same_size;
	static vector<string> _same_mol_type;
	static vector<string> same_protein_name;

	void compute_average_deviations_percentages_in_all_tetra_molecules(int);
	static void sort_and_count_max_sites(vector<string>, string, string, FILE*);
public:
	tetraStructCa(){}
	tetraStructCa(metal_binding_site, bool, string);
	tetraMetalicStructure* form_tetra_site(int);
	static void gen_summary_report_tetra_ca(int);
};

#endif /* TETRA_STRUCTURE_CA_H_ */
