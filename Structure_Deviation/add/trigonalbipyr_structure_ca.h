#ifndef TRIGONALBIPYR_STRUCTURE_CA_H_
#define TRIGONALBIPYR_STRUCTURE_CA_H_

#include <assert.h>
#include "coordinationsite.h"
#include "trigonal_bipyramidal.h"
#include "stdio.h"

class tribipyrStructCa{
	string _accn;
	metal_binding_site _tribipyr_site;
	static string _metal_name;
	static int _no_of_binding_atoms;
	static int _no_of_sides;
	static int _no_of_angles_with_circum_centre;
	static int _no_of_vertex_vertex_angles;
	static int _no_of_features_to_compare;
	static tribipyrMetalicStructure* _tribipyr_struct_standard_ca;
	tribipyrMetalicStructure* _tribipyr_struct_comparable_ca;
	string _site_type_comparable;
	static double _avg_dist_devt_from_cc;
	static double _avg_devt_dist_vert_vert[9];	//Average deviation of edges of two trigonal bipyramidal sites
	static double _avg_devt_angle_vertex_vertex[21];
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

	void compute_average_deviations_percentages_in_all_tribipyr_molecules(int);
	static void sort_and_count_max_sites(vector<string>, string, string, FILE*);
public:
	tribipyrStructCa(){}
	tribipyrStructCa(metal_binding_site, bool, string);
	tribipyrMetalicStructure* form_tribipyr_site(int);
	static void gen_summary_report_tribipyr_ca(int);
};

#endif /* TRIGONALBIPYR_STRUCTURE_CA_H_ */
