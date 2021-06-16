#include "trigonalbipyr_structure_ni.h"

tribipyrMetalicStructure* tribipyrStructNi::_tribipyr_struct_standard_ni=new tribipyrMetalicStructure[2];	//To store Standard, Standard_transformed
string tribipyrStructNi::_metal_name;
int tribipyrStructNi::_no_of_binding_atoms=5;
int tribipyrStructNi::_no_of_sides=9;
int tribipyrStructNi::_no_of_vertex_vertex_angles=21;
int tribipyrStructNi::_no_of_features_to_compare=2;
double tribipyrStructNi::_avg_devt_dist_vert_vert[9]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double tribipyrStructNi::_avg_devt_angle_vertex_vertex[21]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double tribipyrStructNi::_avg_mse[2]={0.0, 0.0};
int tribipyrStructNi::_larger_size=0;
vector<string> tribipyrStructNi::_larger_mol_type;
vector<string> tribipyrStructNi::larger_protein_name;
int tribipyrStructNi::_smaller_size=0;
vector<string> tribipyrStructNi::_smaller_mol_type;
vector<string> tribipyrStructNi::smaller_protein_name;
int tribipyrStructNi::_same_size=0;
vector<string> tribipyrStructNi::_same_mol_type;
vector<string> tribipyrStructNi::same_protein_name;
int tribipyrStructNi::_irregular=0;
vector<string> tribipyrStructNi::_irregular_mol_type;
vector<string> tribipyrStructNi::irregular_protein_name;

tribipyrStructNi::tribipyrStructNi(metal_binding_site met_site, bool first_time_tribipyr_ni, string accn){
	_accn=accn;
	_tribipyr_site=met_site;
	_metal_name=_tribipyr_site.get_metal_atom().get_atom_symbol();

	assert(_no_of_binding_atoms==5 && "No. of Binding Atoms should be 5!!");
	_tribipyr_struct_comparable_ni=new tribipyrMetalicStructure[2];	//Comparable and Comparable_Transformed Trigonal_Bipyramidal_Metal Structure

	assert(_no_of_sides==9 && "No. of Sides should be 9!!");
	assert(_no_of_vertex_vertex_angles==21 && "For Trigonal Bipyramidal No. of Angles of Two Vertices with Another Vertex should be 21!!");
	assert(_no_of_features_to_compare==2 && "For Trigonal Bipyramidal No. of Features to Compare should be 2!!");

	if (first_time_tribipyr_ni==true){
		cout<<"\n\nForm Standard Trigonal Bipyramidal for NI:";
		cout<<"\n\nOriginal Standard Atoms:";
		_tribipyr_struct_standard_ni[0].form_standard_trigonal_bipyramid_NI();
		string metal_name=_tribipyr_site.get_metal_atom().get_atom_symbol();
		string pdb_out_standard_name="standard_original_trigonal_bipyramidal_"+metal_name+".stdn";
		_tribipyr_struct_standard_ni[0].write_atoms_pdb_format(pdb_out_standard_name);	//Generate PDB File for Standard Atoms

		_tribipyr_struct_standard_ni[1]=_tribipyr_struct_standard_ni[0].apply_transformations();
	}
}

tribipyrMetalicStructure* tribipyrStructNi::form_tribipyr_site(int total_tribipyr_sites){
	_site_type_comparable=_tribipyr_site.get_site_type();

	cout<<"\n\nForm Comparable Trigonal Bipyramidal for "<<_metal_name<<":";
	cout<<"\n\nOriginal Comparable Atoms:";

	_tribipyr_struct_comparable_ni[0]=tribipyrMetalicStructure(_tribipyr_site.get_inner_coord_sphere(), _tribipyr_site.get_metal_atom(), _tribipyr_site.get_planar_atoms(), _tribipyr_site.get_non_planar_atoms());
	string pdb_out_comparable_name=_accn+"_comparable_original_trigonal_bipyramidal_"+_metal_name+"_"+_tribipyr_site.get_metal_atom().get_chain_name()+to_string(_tribipyr_site.get_metal_atom().get_residue_no())+".cmp";
	_tribipyr_struct_comparable_ni[0].write_atoms_pdb_format(pdb_out_comparable_name);	//Generate PDB File for Comparable Atoms
	cout<<"\nCheck Point Start";
	_tribipyr_struct_comparable_ni[1]=_tribipyr_struct_comparable_ni[0].apply_transformations();

	string tex_out_name=_accn+"_comparable_transformed_combined_tribipyr_upper_"+_metal_name+"_"+_tribipyr_site.get_metal_atom().get_chain_name()+to_string(_tribipyr_site.get_metal_atom().get_residue_no());
	_tribipyr_struct_standard_ni[1].write_atoms_combined_tex_format_upper(_tribipyr_struct_comparable_ni[1], tex_out_name);

	tex_out_name=_accn+"_comparable_transformed_combined_tribipyr_lower_"+_metal_name+"_"+_tribipyr_site.get_metal_atom().get_chain_name()+to_string(_tribipyr_site.get_metal_atom().get_residue_no());
	_tribipyr_struct_standard_ni[1].write_atoms_combined_tex_format_lower(_tribipyr_struct_comparable_ni[1], tex_out_name);

	_tribipyr_struct_standard_ni[1].compare(_tribipyr_struct_comparable_ni[1]);

	bool check_larger_accn=false, check_smaller_accn=false, check_same_accn=false, check_irregular_accn=false;;

	if (_tribipyr_struct_standard_ni[1].get_structure_size()=="larger"){
		_larger_size++;
		_larger_mol_type.push_back(_tribipyr_site.get_site_type());

		if (check_larger_accn==false){
			int p;
			for (p=0; p<larger_protein_name.size(); p++){
				if (_accn==larger_protein_name.at(p)){
					break;
				}
			}

			if (p==larger_protein_name.size()){
				larger_protein_name.push_back(_accn);
			}

			check_larger_accn=true;
		}
	}else if (_tribipyr_struct_standard_ni[1].get_structure_size()=="smaller"){
		_smaller_size++;
		_smaller_mol_type.push_back(_tribipyr_site.get_site_type());

		if (check_smaller_accn==false){
			int p;
			for (p=0; p<smaller_protein_name.size(); p++){
				if (_accn==smaller_protein_name.at(p)){
					break;
				}
			}

			if (p==smaller_protein_name.size()){
				smaller_protein_name.push_back(_accn);
			}

			check_smaller_accn=true;
		}
	}else if (_tribipyr_struct_standard_ni[1].get_structure_size()=="same"){
		_same_size++;
		_same_mol_type.push_back(_tribipyr_site.get_site_type());

		if (check_same_accn==false){
			int p;
			for (p=0; p<same_protein_name.size(); p++){
				if (_accn==same_protein_name.at(p)){
					break;
				}
			}

			if (p==same_protein_name.size()){
				same_protein_name.push_back(_accn);
			}

			check_same_accn=true;
		}
	}else if (_tribipyr_struct_standard_ni[1].get_structure_size()=="irregular"){
		_irregular++;
		_irregular_mol_type.push_back(_tribipyr_site.get_site_type());

		if (check_irregular_accn==false){
			int p;
			for (p=0; p<irregular_protein_name.size(); p++){
				if (_accn==irregular_protein_name.at(p)){
					break;
				}
			}

			if (p==irregular_protein_name.size()){
				irregular_protein_name.push_back(_accn);
			}

			check_irregular_accn=true;
		}
	}

	compute_average_deviations_percentages_in_all_tribipyr_molecules(total_tribipyr_sites);

	pdb_out_comparable_name=_accn+"_comparable_different_trigonal_bipyramidal_"+_metal_name+"_"+_tribipyr_site.get_metal_atom().get_chain_name()+to_string(_tribipyr_site.get_metal_atom().get_residue_no())+".comp";
	_tribipyr_struct_standard_ni[0].comparison(_tribipyr_struct_comparable_ni[0], pdb_out_comparable_name);

	string out_file_name="Deviation_Report_Trigonal_Biyramidal_"+_accn+"_"+_metal_name+"_"+_tribipyr_site.get_metal_atom().get_chain_name()+to_string(_tribipyr_site.get_metal_atom().get_residue_no())+".mbsi";
	FILE* fp1=fopen(out_file_name.c_str(), "w");
	cout<<"\n\nStarts reports generation for "<<_metal_name<<" Ions:\n";
	cout <<"\nOutput Report File Name: " <<out_file_name<< endl;

	string site_type_standard=_tribipyr_struct_standard_ni[0].assign_standard_site_type();
	_tribipyr_struct_standard_ni[1].gen_report_metal_binding_sites(_tribipyr_struct_standard_ni[0], _tribipyr_struct_comparable_ni[0], _tribipyr_struct_comparable_ni[1], _metal_name, site_type_standard, _tribipyr_site.get_site_type(), fp1);
	fclose(fp1);

	cout<<"\nEnds report generation for trigonal_bipyramidal "+_accn+"_"+_metal_name+"_"+_tribipyr_site.get_metal_atom().get_chain_name()+to_string(_tribipyr_site.get_metal_atom().get_residue_no())+" Ions\n";

	return _tribipyr_struct_comparable_ni;
}

void tribipyrStructNi::compute_average_deviations_percentages_in_all_tribipyr_molecules(int total_tribipyr_sites){
	double* temp_side_dev=_tribipyr_struct_standard_ni[1].get_side_deviations();
	double* temp_angle_vv=_tribipyr_struct_standard_ni[1].get_angle_deviations_vertex_vertex();
	double* temp_mse=_tribipyr_struct_standard_ni[1].get_mse();

	if (total_tribipyr_sites==1){
		for (int p=0; p<_no_of_sides; p++){
			_avg_devt_dist_vert_vert[p]=temp_side_dev[p];
		}

		for (int p=0; p<_no_of_vertex_vertex_angles; p++){
			_avg_devt_angle_vertex_vertex[p]=temp_angle_vv[p];
		}

		for (int p=0; p<_no_of_features_to_compare; p++){
			_avg_mse[p]=temp_mse[p];
		}
	}else{
		for (int p=0; p<_no_of_sides; p++){
			_avg_devt_dist_vert_vert[p]=(_avg_devt_dist_vert_vert[p]+temp_side_dev[p])/2.0;
		}

		for (int p=0; p<_no_of_vertex_vertex_angles; p++){
			_avg_devt_angle_vertex_vertex[p]=(_avg_devt_angle_vertex_vertex[p]+temp_angle_vv[p])/2.0;
		}

		for (int p=0; p<_no_of_features_to_compare; p++){
			_avg_mse[p]=(_avg_mse[p]+temp_mse[p])/2.0;
		}
	}
}

void tribipyrStructNi::gen_summary_report_tribipyr_ni(int total_tribipyr_sites){
	string summary_file_tribipyr_ni="Summary_Deviation_Report_NI_trigonal_bipyramidal.mbsi";
	FILE *fp_tribipyr_ni=fopen(summary_file_tribipyr_ni.c_str(), "w");

	fprintf(fp_tribipyr_ni, "HEADER     PROGRAM NAME: StructureDeviation   VERSION: 1.1");
	fprintf(fp_tribipyr_ni, "\nTITLE      Summary Report for Computation of Structural Deviations of Two Trigonal Bipyramidal NI Binding Sites with Coordination Number 5");

	fprintf(fp_tribipyr_ni, "\nREMARK     1");
	fprintf(fp_tribipyr_ni, "\nREMARK     1  No. of Trigonal Bipyramidal Sites Compared with Standard NI ION: %d", total_tribipyr_sites);

	if (larger_protein_name.size()!=0){
		fprintf(fp_tribipyr_ni, "\nREMARK     2");
		fprintf(fp_tribipyr_ni, "\nREMARK     2  No. of Trigonal Bipyramidal Sites whose Size is Larger than the Standard Site: %d", _larger_size);
		fprintf(fp_tribipyr_ni, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<larger_protein_name.size(); p++){
			fprintf(fp_tribipyr_ni, "%s ", larger_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_larger_mol_type, "Larger", "REMARK     2", fp_tribipyr_ni);
	}else{
		fprintf(fp_tribipyr_ni, "\nREMARK     2  No Such Site Exists  whose Size is Larger than the Standard Site");
	}

	if (smaller_protein_name.size()!=0){
		fprintf(fp_tribipyr_ni, "\nREMARK     2");
		fprintf(fp_tribipyr_ni, "\nREMARK     2  No. of Trigonal Bipyramidal Sites whose Size is Smaller than the Standard Site: %d", _smaller_size);
		fprintf(fp_tribipyr_ni, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<smaller_protein_name.size(); p++){
			fprintf(fp_tribipyr_ni, "%s ", smaller_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_smaller_mol_type, "Smaller", "REMARK     2", fp_tribipyr_ni);
	}else{
		fprintf(fp_tribipyr_ni, "\nREMARK     2  No Such Site Exists  whose Size is Smaller than the Standard Site");
	}

	if (same_protein_name.size()!=0){
		fprintf(fp_tribipyr_ni, "\nREMARK     2");
		fprintf(fp_tribipyr_ni, "\nREMARK     2  No. of Trigonal Bipyramidal Sites whose Size is Same as that of the Standard Site: %d", _same_size);
		fprintf(fp_tribipyr_ni, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<same_protein_name.size(); p++){
			fprintf(fp_tribipyr_ni, "%s ", same_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_same_mol_type, "Same", "REMARK     2", fp_tribipyr_ni);
	}else{
		fprintf(fp_tribipyr_ni, "\nREMARK     2  No Such Site Exists  whose Size is Same as that of the Standard Site");
	}

	if (irregular_protein_name.size()!=0){
		fprintf(fp_tribipyr_ni, "\nREMARK     2");
		fprintf(fp_tribipyr_ni, "\nREMARK     2  No. of Trigonal Bipyramidal Sites whose Size is irregular in comparison to the Standard Site: %d", _irregular);
		fprintf(fp_tribipyr_ni, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<irregular_protein_name.size(); p++){
			fprintf(fp_tribipyr_ni, "%s ", irregular_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_irregular_mol_type, "Irregular", "REMARK     2", fp_tribipyr_ni);
	}else{
		fprintf(fp_tribipyr_ni, "\nREMARK     2  No Such Site Exists  whose Size is Irregular in comparison to the Standard Site");
	}

	fprintf(fp_tribipyr_ni, "\nREMARK     3");
	fprintf(fp_tribipyr_ni, "\nREMARK     3  Average Deviation of Distances of Atom-to-Atom:                         ");
	for (int p=0; p<_no_of_sides; p++){
		fprintf(fp_tribipyr_ni, "\nSide %d: %8.3lf", (p+1), _avg_devt_dist_vert_vert[p]);
	}

	fprintf(fp_tribipyr_ni, "\nREMARK     3         Average MSE Distance AA: %8.3lf", _avg_mse[0]);

	fprintf(fp_tribipyr_ni, "\nREMARK     3");
	fprintf(fp_tribipyr_ni, "\nREMARK     3  Average  Deviation of Angles of Made by Two Atoms with Another Atom:                         ");
	for (int p=0; p<_no_of_vertex_vertex_angles; p++){
		fprintf(fp_tribipyr_ni, "\nAngle %d: %8.3lf", (p+1), _avg_devt_angle_vertex_vertex[p]);
	}

	fprintf(fp_tribipyr_ni, "\nREMARK     3         Average MSE Angle AA: %8.3lf", _avg_mse[1]);

	fclose(fp_tribipyr_ni);
}

void tribipyrStructNi::sort_and_count_max_sites(vector<string> mol_type, string qualifier, string remark, FILE* file_name){
	vector<int> count_type;
	vector<string> unique_type;

	for (int x=0; x<mol_type.size(); x++){
		int k;
		for (k=0; k<unique_type.size(); k++){
			if (mol_type.at(x)==unique_type.at(k)){
				break;
			}
		}

		if (k==unique_type.size()){
			unique_type.push_back(mol_type.at(x));
			int c=1;
			for (int y=0; y<mol_type.size(); y++){
				if (x<y){
					if (mol_type.at(x)==mol_type.at(y)){
						c++;
					}
				}
			}

			count_type.push_back(c);
		}
	}

	for (int x=0; x<count_type.size()-1; x++){
		for (int y=(x+1); y<count_type.size(); y++){
			if (count_type.at(x)<count_type.at(y)){
				int c=count_type.at(x);
				count_type.at(x)=count_type.at(y);
				count_type.at(y)=c;

				string st=unique_type.at(x);
				unique_type.at(x)=unique_type.at(y);
				unique_type.at(y)=st;
			}
		}
	}

	fprintf(file_name, "\n%s Max. %s Molecule Type: %s; No. of Such Sites: %d", remark.c_str(), qualifier.c_str(), unique_type.at(0).c_str(), count_type.at(0));
}
