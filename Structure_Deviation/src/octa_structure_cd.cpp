#include "octa_structure_cd.h"

octaMetalicStructure* octaStructCd::_octa_struct_standard_cd=new octaMetalicStructure[2];	//To store Standard, Standard_transformed
string octaStructCd::_metal_name;
int octaStructCd::_no_of_binding_atoms=6;
int octaStructCd::_no_of_sides=12;
int octaStructCd::_no_of_vertex_vertex_angles=28;
int octaStructCd::_no_of_features_to_compare=2;
double octaStructCd::_avg_devt_dist_vert_vert[12]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double octaStructCd::_avg_devt_angle_vertex_vertex[28]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double octaStructCd::_avg_mse[2]={0.0, 0.0};
int octaStructCd::_larger_size=0;
vector<string> octaStructCd::_larger_mol_type;
vector<string> octaStructCd::larger_protein_name;
int octaStructCd::_smaller_size=0;
vector<string> octaStructCd::_smaller_mol_type;
vector<string> octaStructCd::smaller_protein_name;
int octaStructCd::_same_size=0;
vector<string> octaStructCd::_same_mol_type;
vector<string> octaStructCd::same_protein_name;
int octaStructCd::_irregular=0;
vector<string> octaStructCd::_irregular_mol_type;
vector<string> octaStructCd::irregular_protein_name;

octaStructCd::octaStructCd(metal_binding_site met_site, bool first_time_octa_cd, string accn){
	_accn=accn;
	_octa_site=met_site;
	_metal_name=_octa_site.get_metal_atom().get_atom_symbol();

	assert(_no_of_binding_atoms==6 && "No. of Binding Atoms should be 6!!");
	_octa_struct_comparable_cd=new octaMetalicStructure[2];	//Comparable and Comparable_Transformed Octahedral Metal Structure

	assert(_no_of_sides==12 && "No. of Sides should be 12!!");
	assert(_no_of_vertex_vertex_angles==28 && "For Octahedral No. of Angles of Two Vertices with Another Vertex should be 28!!");
	assert(_no_of_features_to_compare==2 && "For Octahedral No. of Features to Compare should be 2!!");

	if (first_time_octa_cd==true){
		cout<<"\n\nForm Standard Octahedron for CD:";
		cout<<"\n\nOriginal Standard Atoms:";
		_octa_struct_standard_cd[0].form_standard_octa_CD();
		string metal_name=_octa_site.get_metal_atom().get_atom_symbol();
		string pdb_out_standard_name="standard_original_octahedral_"+metal_name+".stdn";
		_octa_struct_standard_cd[0].write_atoms_pdb_format(pdb_out_standard_name);	//Generate PDB File for Standard Atoms

		_octa_struct_standard_cd[1]=_octa_struct_standard_cd[0].apply_transformations();
	}
}

octaMetalicStructure* octaStructCd::form_octa_site(int total_octa_sites){
	_site_type_comparable=_octa_site.get_site_type();

	cout<<"\n\nForm Comparable Octahedral for "<<_metal_name<<":";
	cout<<"\n\nOriginal Comparable Atoms:";

	_octa_struct_comparable_cd[0]=octaMetalicStructure(_octa_site.get_inner_coord_sphere(), _octa_site.get_metal_atom(), _octa_site.get_planar_atoms(), _octa_site.get_non_planar_atoms());
	string pdb_out_comparable_name=_accn+"_comparable_original_octahedral_"+_metal_name+"_"+_octa_site.get_metal_atom().get_chain_name()+to_string(_octa_site.get_metal_atom().get_residue_no())+".cmp";
	_octa_struct_comparable_cd[0].write_atoms_pdb_format(pdb_out_comparable_name);	//Generate PDB File for Comparable Atoms

	_octa_struct_comparable_cd[1]=_octa_struct_comparable_cd[0].apply_transformations();

	string tex_out_name=_accn+"_comparable_transformed_combined_octa_upper_left_"+_metal_name+"_"+_octa_site.get_metal_atom().get_chain_name()+to_string(_octa_site.get_metal_atom().get_residue_no());
	_octa_struct_standard_cd[1].write_atoms_combined_tex_format_upper_left(_octa_struct_comparable_cd[1], tex_out_name);

	tex_out_name=_accn+"_comparable_transformed_combined_octa_upper_right_"+_metal_name+"_"+_octa_site.get_metal_atom().get_chain_name()+to_string(_octa_site.get_metal_atom().get_residue_no());
	_octa_struct_standard_cd[1].write_atoms_combined_tex_format_upper_right(_octa_struct_comparable_cd[1], tex_out_name);

	tex_out_name=_accn+"_comparable_transformed_combined_octa_lower_left_"+_metal_name+"_"+_octa_site.get_metal_atom().get_chain_name()+to_string(_octa_site.get_metal_atom().get_residue_no());
	_octa_struct_standard_cd[1].write_atoms_combined_tex_format_lower_left(_octa_struct_comparable_cd[1], tex_out_name);

	tex_out_name=_accn+"_comparable_transformed_combined_octa_lower_right_"+_metal_name+"_"+_octa_site.get_metal_atom().get_chain_name()+to_string(_octa_site.get_metal_atom().get_residue_no());
	_octa_struct_standard_cd[1].write_atoms_combined_tex_format_lower_right(_octa_struct_comparable_cd[1], tex_out_name);

	_octa_struct_standard_cd[1].compare(_octa_struct_comparable_cd[1]);

	bool check_larger_accn=false, check_smaller_accn=false, check_same_accn=false, check_irregular_accn=false;;

	if (_octa_struct_standard_cd[1].get_structure_size()=="larger"){
		_larger_size++;
		_larger_mol_type.push_back(_octa_site.get_site_type());

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
	}else if (_octa_struct_standard_cd[1].get_structure_size()=="smaller"){
		_smaller_size++;
		_smaller_mol_type.push_back(_octa_site.get_site_type());

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
	}else if (_octa_struct_standard_cd[1].get_structure_size()=="same"){
		_same_size++;
		_same_mol_type.push_back(_octa_site.get_site_type());

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
	}else if (_octa_struct_standard_cd[1].get_structure_size()=="irregular"){
		_irregular++;
		_irregular_mol_type.push_back(_octa_site.get_site_type());

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

	compute_average_deviations_percentages_in_all_octa_molecules(total_octa_sites);

	pdb_out_comparable_name=_accn+"_comparable_different_octahedral_"+_metal_name+"_"+_octa_site.get_metal_atom().get_chain_name()+to_string(_octa_site.get_metal_atom().get_residue_no())+".comp";
	_octa_struct_standard_cd[0].comparison(_octa_struct_comparable_cd[0], pdb_out_comparable_name);

	string out_file_name="Deviation_Report_Octahedral_"+_accn+"_"+_metal_name+"_"+_octa_site.get_metal_atom().get_chain_name()+to_string(_octa_site.get_metal_atom().get_residue_no())+".mbsi";
	FILE* fp1=fopen(out_file_name.c_str(), "w");
	cout<<"\n\nStarts reports generation for "<<_metal_name<<" Ions:\n";
	cout <<"\nOutput Report File Name: " <<out_file_name<< endl;

	string site_type_standard=_octa_struct_standard_cd[0].assign_standard_site_type();
	_octa_struct_standard_cd[1].gen_report_metal_binding_sites(_octa_struct_standard_cd[0], _octa_struct_comparable_cd[0], _octa_struct_comparable_cd[1], _metal_name, site_type_standard, _octa_site.get_site_type(), fp1);
	fclose(fp1);

	cout<<"\nEnds report generation for Octahedral "+_accn+"_"+_metal_name+"_"+_octa_site.get_metal_atom().get_chain_name()+to_string(_octa_site.get_metal_atom().get_residue_no())+" Ions\n";

	return _octa_struct_comparable_cd;
}

void octaStructCd::compute_average_deviations_percentages_in_all_octa_molecules(int total_octa_sites){
	double* temp_side_dev=_octa_struct_standard_cd[1].get_side_deviations();
	double* temp_angle_vv=_octa_struct_standard_cd[1].get_angle_deviations_vertex_vertex();
	double* temp_mse=_octa_struct_standard_cd[1].get_mse();

	if (total_octa_sites==1){
		for (int p=0; p<_no_of_sides; p++){
			_avg_devt_dist_vert_vert[p]=temp_side_dev[p];
		}

		for (int p=0; p<_no_of_vertex_vertex_angles; p++){
			cout<<"\nIN "<<_metal_name<<" SITE:"<<temp_angle_vv[p];
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
			cout<<"\nIN "<<_metal_name<<" SITE:"<<temp_angle_vv[p];
			_avg_devt_angle_vertex_vertex[p]=(_avg_devt_angle_vertex_vertex[p]+temp_angle_vv[p])/2.0;
		}

		for (int p=0; p<_no_of_features_to_compare; p++){
			_avg_mse[p]=(_avg_mse[p]+temp_mse[p])/2.0;
		}
	}
}

void octaStructCd::gen_summary_report_octa_cd(int total_octa_sites){
	string summary_file_octa_cd="Summary_Deviation_Report_CD_octahedral.mbsi";
	FILE *fp_octa_cd=fopen(summary_file_octa_cd.c_str(), "w");

	fprintf(fp_octa_cd, "HEADER     PROGRAM NAME: StructureDeviation  VERSION: 1.1");
	fprintf(fp_octa_cd, "\nTITLE      Summary Report for Computation of Structural Deviations of Two Octahedral CD Binding Sites with Coordination Number 6");

	fprintf(fp_octa_cd, "\nREMARK     1");
	fprintf(fp_octa_cd, "\nREMARK     1  No. of Octahedral Sites Compared with Standard CD ION: %d", total_octa_sites);

	if (larger_protein_name.size()!=0){
		fprintf(fp_octa_cd, "\nREMARK     2");
		fprintf(fp_octa_cd, "\nREMARK     2  No. of Octahedral Sites whose Size is Larger than the Standard Site: %d", _larger_size);
		fprintf(fp_octa_cd, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<larger_protein_name.size(); p++){
			fprintf(fp_octa_cd, "%s ", larger_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_larger_mol_type, "Larger", "REMARK     2", fp_octa_cd);
	}else{
		fprintf(fp_octa_cd, "\nREMARK     2  No Such Site Exists  whose Size is Larger than the Standard Site");
	}

	if (smaller_protein_name.size()!=0){
		fprintf(fp_octa_cd, "\nREMARK     2");
		fprintf(fp_octa_cd, "\nREMARK     2  No. of Octahedral Sites whose Size is Smaller than the Standard Site: %d", _smaller_size);
		fprintf(fp_octa_cd, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<smaller_protein_name.size(); p++){
			fprintf(fp_octa_cd, "%s ", smaller_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_smaller_mol_type, "Smaller", "REMARK     2", fp_octa_cd);
	}else{
		fprintf(fp_octa_cd, "\nREMARK     2  No Such Site Exists  whose Size is Smaller than the Standard Site");
	}

	if (same_protein_name.size()!=0){
		fprintf(fp_octa_cd, "\nREMARK     2");
		fprintf(fp_octa_cd, "\nREMARK     2  No. of Octahedral Sites whose Size is Same as that of the Standard Site: %d", _same_size);
		fprintf(fp_octa_cd, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<same_protein_name.size(); p++){
			fprintf(fp_octa_cd, "%s ", same_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_same_mol_type, "Same", "REMARK     2", fp_octa_cd);
	}else{
		fprintf(fp_octa_cd, "\nREMARK     2  No Such Site Exists  whose Size is Same as that of the Standard Site");
	}

	if (irregular_protein_name.size()!=0){
		fprintf(fp_octa_cd, "\nREMARK     2");
		fprintf(fp_octa_cd, "\nREMARK     2  No. of Octahedral Sites whose Size is irregular in comparison to the Standard Site: %d", _irregular);
		fprintf(fp_octa_cd, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<irregular_protein_name.size(); p++){
			fprintf(fp_octa_cd, "%s ", irregular_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_irregular_mol_type, "Irregular", "REMARK     2", fp_octa_cd);
	}else{
		fprintf(fp_octa_cd, "\nREMARK     2  No Such Site Exists  whose Size is Irregular in comparison to the Standard Site");
	}

	fprintf(fp_octa_cd, "\nREMARK     3");
	fprintf(fp_octa_cd, "\nREMARK     3  Average Deviation of Distances of Atom-to-Atom:                         ");
	for (int p=0; p<_no_of_sides; p++){
		fprintf(fp_octa_cd, "\nSide %d: %8.3lf", (p+1), _avg_devt_dist_vert_vert[p]);
	}

	fprintf(fp_octa_cd, "\nREMARK     3         Average MSE Distance AA: %8.3lf", _avg_mse[0]);

	fprintf(fp_octa_cd, "\nREMARK     3");
	fprintf(fp_octa_cd, "\nREMARK     3  Average  Deviation of Angles Made by Two Atoms with Another Atom:                         ");
	for (int p=0; p<_no_of_vertex_vertex_angles; p++){
		fprintf(fp_octa_cd, "\nAngle %d: %8.3lf", (p+1), _avg_devt_angle_vertex_vertex[p]);
	}

	fprintf(fp_octa_cd, "\nREMARK     3         Average MSE Angle AA: %8.3lf", _avg_mse[1]);

	fclose(fp_octa_cd);
}

void octaStructCd::sort_and_count_max_sites(vector<string> mol_type, string qualifier, string remark, FILE* file_name){
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

