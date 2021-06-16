
#ifndef TETRAGONALPL_STRUCTURE_CD_H_
#define TETRAGONALPL_STRUCTURE_CD_H_

#include <assert.h>
#include "coordinationsite.h"
#include "tetragonalplanar.h"
#include "stdio.h"

class tetraplStructCd{
	string _accn;
	metal_binding_site _tetrapl_site;
	static int _no_of_binding_atoms;
	static int _no_of_sides;
	static int _no_of_angles_with_circum_centre;
	static int _no_of_vertex_vertex_angles;
	static int _no_of_features_to_compare;
	static tetraplMetalicStructure* _tetrapl_struct_standard_cd;
	tetraplMetalicStructure* _tetrapl_struct_comparable_cd;
	string _site_type_comparable;
	static double _avg_dist_devt_from_cc;
	static double _avg_devt_dist_vert_vert[];	//Average deviation of edges of two triangles
	static double _avg_devt_angle_with_cc[];	//Average deviation of angles of two triangles
	static double _avg_devt_angle_vertex_vertex[];
	static double _avg_mse[];	//Average Mean Square Error due to distance deviations and angle deviations
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

	void compute_average_deviations_in_all_tetrapl_molecules(int);
	static void sort_and_count_max_sites(vector<string>, string, string, FILE*);
public:
	tetraplStructCd(){}
	tetraplStructCd(metal_binding_site, bool, string);
	tetraplMetalicStructure* form_tetrapl_site(int);
	static void gen_summary_report_tetrapl_cd(int);
};

#endif /* TRIGONALPL_STRUCTURE_CD_H_ */
