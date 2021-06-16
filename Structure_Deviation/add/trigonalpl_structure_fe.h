
#ifndef TRIGONALPL_STRUCTURE_FE_H_
#define TRIGONALPL_STRUCTURE_FE_H_

#include <assert.h>
#include "coordinationsite.h"
#include "trigonalplanar.h"
#include "stdio.h"

class triplStructFe{
	string _accn;
	metal_binding_site _tripl_site;
	static int _no_of_binding_atoms;
	static int _no_of_sides;
	static int _no_of_angles_with_circum_centre;
	static int _no_of_vertex_vertex_angles;
	static int _no_of_features_to_compare;
	static triplMetalicStructure* _tripl_struct_standard_fe;
	triplMetalicStructure* _tripl_struct_comparable_fe;
	string _site_type_comparable;
	static double _avg_dist_devt_from_cc;
	static double _avg_devt_dist_vert_vert[];	//Average deviation of edges of two triangles
	static double _avg_devt_angle_with_cc[];	//Average deviation of angles of two triangles
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

	void compute_average_deviations_percentages_in_all_tripl_molecules(int);
	static void sort_and_count_max_sites(vector<string>, string, string, FILE*);
public:
	triplStructFe(){}
	triplStructFe(metal_binding_site, bool, string);
	triplMetalicStructure* form_tripl_site(int);
	static void gen_summary_report_tripl_fe(int);
};

#endif /* TRIGONALPL_STRUCTURE_FE_H_ */
