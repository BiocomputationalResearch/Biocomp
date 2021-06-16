#ifndef OCTA_STRUCTURE_ZN_H_
#define OCTA_STRUCTURE_ZN_H_

#include <assert.h>
#include "coordinationsite.h"
#include "octahedral.h"
#include "stdio.h"

class octaStructZn{
	string _accn;
	metal_binding_site _octa_site;
	static string _metal_name;
	static int _no_of_binding_atoms;
	static int _no_of_sides;
	static int _no_of_angles_with_circum_centre;
	static int _no_of_vertex_vertex_angles;
	static int _no_of_features_to_compare;
	static octaMetalicStructure* _octa_struct_standard_zn;
	octaMetalicStructure* _octa_struct_comparable_zn;
	string _site_type_comparable;
	static double _avg_dist_devt_from_cc;
	static double _avg_devt_dist_vert_vert[12];	//Average deviation of edges of two octahedral sites
	static double _avg_devt_angle_vertex_vertex[28];
	static double _avg_mse[2];	//Average Mean Square Error due to distance deviations and angle deviations
	static int _larger_size;
	static vector<string> _larger_mol_type;
	static vector<string> larger_protein_name;
	static int _smaller_size;
	static vector<string> _smaller_mol_type;
	static vector<string> smaller_protein_name;
	static int _same_size;
	static vector<string> _same_mol_type;
	static vector<string> same_protein_name;
	static int _irregular;
	static vector<string> _irregular_mol_type;
	static vector<string> irregular_protein_name;

	void compute_average_deviations_percentages_in_all_octa_molecules(int);
	static void sort_and_count_max_sites(vector<string>, string, string, FILE*);
public:
	octaStructZn(){}
	octaStructZn(metal_binding_site, bool, string);
	octaMetalicStructure* form_octa_site(int);
	static void gen_summary_report_octa_zn(int);
};

#endif /* OCTA_STRUCTURE_ZN_H_ */
