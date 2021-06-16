#include "tetragonalpyr_structure_na.h"

tetrapyrMetalicStructure* tetrapyrStructNa::_tetrapyr_struct_standard_na=new tetrapyrMetalicStructure[2];	//To store Standard, Standard_transformed
string tetrapyrStructNa::_metal_name;
int tetrapyrStructNa::_no_of_binding_atoms=5;
int tetrapyrStructNa::_no_of_sides=8;
int tetrapyrStructNa::_no_of_vertex_vertex_angles=16;
int tetrapyrStructNa::_no_of_features_to_compare=2;
double tetrapyrStructNa::_avg_devt_dist_vert_vert[8]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double tetrapyrStructNa::_avg_devt_angle_vertex_vertex[16]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double tetrapyrStructNa::_avg_mse[2]={0.0, 0.0};
int tetrapyrStructNa::_larger_size=0;
vector<string> tetrapyrStructNa::_larger_mol_type;
vector<string> tetrapyrStructNa::larger_protein_name;
int tetrapyrStructNa::_smaller_size=0;
vector<string> tetrapyrStructNa::_smaller_mol_type;
vector<string> tetrapyrStructNa::smaller_protein_name;
int tetrapyrStructNa::_same_size=0;
vector<string> tetrapyrStructNa::_same_mol_type;
vector<string> tetrapyrStructNa::same_protein_name;
int tetrapyrStructNa::_irregular=0;
vector<string> tetrapyrStructNa::_irregular_mol_type;
vector<string> tetrapyrStructNa::irregular_protein_name;

tetrapyrStructNa::tetrapyrStructNa(metal_binding_site met_site, bool first_time_tetrapyr_na, string accn){
	_accn=accn;
	_tetrapyr_site=met_site;
	_metal_name=_tetrapyr_site.get_metal_atom().get_atom_symbol();

	assert(_no_of_binding_atoms==5 && "No. of Binding Atoms should be 5!!");
	_tetrapyr_struct_comparable_na=new tetrapyrMetalicStructure[2];	//Comparable and Comparable_Transformed Tetragonal_Pyramidal_Metal Structure

	assert(_no_of_sides==8 && "No. of Sides should be 8!!");
	assert(_no_of_vertex_vertex_angles==16 && "For Tetragonal Pyramidal No. of Angles of Two Vertices with Another Vertex should be 16!!");
	assert(_no_of_features_to_compare==2 && "For Tetragonal Pyramidal No. of Features to Compare should be 2!!");

	if (first_time_tetrapyr_na==true){
		cout<<"\n\nForm Standard Tetragonal Pyramid for NA:";
		cout<<"\n\nOriginal Standard Atoms:";
		_tetrapyr_struct_standard_na[0].form_standard_tetragonal_pyramid_NA();
		string metal_name=_tetrapyr_site.get_metal_atom().get_atom_symbol();
		string pdb_out_standard_name="standard_original_tetragonal_pyramidal_"+metal_name+".stdn";
		_tetrapyr_struct_standard_na[0].write_atoms_pdb_format(pdb_out_standard_name);	//Generate PDB File for Standard Atoms

		_tetrapyr_struct_standard_na[1]=_tetrapyr_struct_standard_na[0].apply_transformations();
	}
}

tetrapyrMetalicStructure* tetrapyrStructNa::form_tetrapyr_site(int total_tetrapyr_sites){
	_site_type_comparable=_tetrapyr_site.get_site_type();

	cout<<"\n\nForm Comparable Tetragonal Pyramidal for "<<_metal_name<<":";
	cout<<"\n\nOriginal Comparable Atoms:";

	_tetrapyr_struct_comparable_na[0]=tetrapyrMetalicStructure(_tetrapyr_site.get_inner_coord_sphere(), _tetrapyr_site.get_metal_atom(), _tetrapyr_site.get_planar_atoms(), _tetrapyr_site.get_non_planar_atoms());
	string pdb_out_comparable_name=_accn+"_comparable_original_tetragonal_pyramidal_"+_metal_name+"_"+_tetrapyr_site.get_metal_atom().get_chain_name()+to_string(_tetrapyr_site.get_metal_atom().get_residue_no())+".cmp";
	_tetrapyr_struct_comparable_na[0].write_atoms_pdb_format(pdb_out_comparable_name);	//Generate PDB File for Comparable Atoms

	_tetrapyr_struct_comparable_na[1]=_tetrapyr_struct_comparable_na[0].apply_transformations();

	string tex_out_name=_accn+"_comparable_transformed_combined_tetrapyr_left_"+_metal_name+"_"+_tetrapyr_site.get_metal_atom().get_chain_name()+to_string(_tetrapyr_site.get_metal_atom().get_residue_no());
	_tetrapyr_struct_standard_na[1].write_atoms_combined_tex_format_left(_tetrapyr_struct_comparable_na[1], tex_out_name);

	tex_out_name=_accn+"_comparable_transformed_combined_tetrapyr_right_"+_metal_name+"_"+_tetrapyr_site.get_metal_atom().get_chain_name()+to_string(_tetrapyr_site.get_metal_atom().get_residue_no());
	_tetrapyr_struct_standard_na[1].write_atoms_combined_tex_format_right(_tetrapyr_struct_comparable_na[1], tex_out_name);

	_tetrapyr_struct_standard_na[1].compare(_tetrapyr_struct_comparable_na[1]);

	bool check_larger_accn=false, check_smaller_accn=false, check_same_accn=false, check_irregular_accn=false;;

	if (_tetrapyr_struct_standard_na[1].get_structure_size()=="larger"){
		_larger_size++;
		_larger_mol_type.push_back(_tetrapyr_site.get_site_type());

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
	}else if (_tetrapyr_struct_standard_na[1].get_structure_size()=="smaller"){
		_smaller_size++;
		_smaller_mol_type.push_back(_tetrapyr_site.get_site_type());

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
	}else if (_tetrapyr_struct_standard_na[1].get_structure_size()=="same"){
		_same_size++;
		_same_mol_type.push_back(_tetrapyr_site.get_site_type());

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
	}else if (_tetrapyr_struct_standard_na[1].get_structure_size()=="irregular"){
		_irregular++;
		_irregular_mol_type.push_back(_tetrapyr_site.get_site_type());

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

	compute_average_deviations_percentages_in_all_tetrapyr_molecules(total_tetrapyr_sites);

	pdb_out_comparable_name=_accn+"_comparable_different_tetragonal_pyramidal_"+_metal_name+"_"+_tetrapyr_site.get_metal_atom().get_chain_name()+to_string(_tetrapyr_site.get_metal_atom().get_residue_no())+".comp";
	_tetrapyr_struct_standard_na[0].comparison(_tetrapyr_struct_comparable_na[0], pdb_out_comparable_name);

	string out_file_name="Deviation_Report_Tetragonal_Pyramidal_"+_accn+"_"+_metal_name+"_"+_tetrapyr_site.get_metal_atom().get_chain_name()+to_string(_tetrapyr_site.get_metal_atom().get_residue_no())+".mbsi";
	FILE* fp1=fopen(out_file_name.c_str(), "w");
	cout<<"\n\nStarts reports generation for "<<_metal_name<<" Ions:\n";
	cout <<"\nOutput Report File Name: " <<out_file_name<< endl;

	string site_type_standard=_tetrapyr_struct_standard_na[0].assign_standard_site_type();
	_tetrapyr_struct_standard_na[1].gen_report_metal_binding_sites(_tetrapyr_struct_standard_na[0], _tetrapyr_struct_comparable_na[0], _tetrapyr_struct_comparable_na[1], _metal_name, site_type_standard, _tetrapyr_site.get_site_type(), fp1);
	fclose(fp1);

	cout<<"\nEnds report generation for Tetragonal Pyramidal "+_accn+"_"+_metal_name+"_"+_tetrapyr_site.get_metal_atom().get_chain_name()+to_string(_tetrapyr_site.get_metal_atom().get_residue_no())+" Ions\n";

	return _tetrapyr_struct_comparable_na;
}

void tetrapyrStructNa::compute_average_deviations_percentages_in_all_tetrapyr_molecules(int total_tetrapyr_sites){
	double* temp_side_dev=_tetrapyr_struct_standard_na[1].get_side_deviations();
	double* temp_angle_vv=_tetrapyr_struct_standard_na[1].get_angle_deviations_vertex_vertex();
	double* temp_mse=_tetrapyr_struct_standard_na[1].get_mse();

	if (total_tetrapyr_sites==1){
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

void tetrapyrStructNa::gen_summary_report_tetrapyr_na(int total_tetrapyr_sites){
	string summary_file_tetrapyr_na="Summary_Deviation_Report_NA_tetragonal_pyramidal.mbsi";
	FILE *fp_tetrapyr_na=fopen(summary_file_tetrapyr_na.c_str(), "w");

	fprintf(fp_tetrapyr_na, "HEADER     PROGRAM NAME: StructureDeviation  VERSION: 1.1");
	fprintf(fp_tetrapyr_na, "\nTITLE      Summary Report for Computation of Structural Deviations of Two Tetragonal Pyramidal NA Binding Sites with Coordination Number 5");

	fprintf(fp_tetrapyr_na, "\nREMARK     1");
	fprintf(fp_tetrapyr_na, "\nREMARK     1  No. of Tetragonal Pyramidal Sites Compared with Standard NA ION: %d", total_tetrapyr_sites);

	if (larger_protein_name.size()!=0){
		fprintf(fp_tetrapyr_na, "\nREMARK     2");
		fprintf(fp_tetrapyr_na, "\nREMARK     2  No. of Tetragonal Pyramidal Sites whose Size is Larger than the Standard Site: %d", _larger_size);
		fprintf(fp_tetrapyr_na, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<larger_protein_name.size(); p++){
			fprintf(fp_tetrapyr_na, "%s ", larger_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_larger_mol_type, "Larger", "REMARK     2", fp_tetrapyr_na);
	}else{
		fprintf(fp_tetrapyr_na, "\nREMARK     2  No Such Site Exists  whose Size is Larger than the Standard Site");
	}

	if (smaller_protein_name.size()!=0){
		fprintf(fp_tetrapyr_na, "\nREMARK     2");
		fprintf(fp_tetrapyr_na, "\nREMARK     2  No. of Tetragonal Pyramidal Sites whose Size is Smaller than the Standard Site: %d", _smaller_size);
		fprintf(fp_tetrapyr_na, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<smaller_protein_name.size(); p++){
			fprintf(fp_tetrapyr_na, "%s ", smaller_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_smaller_mol_type, "Smaller", "REMARK     2", fp_tetrapyr_na);
	}else{
		fprintf(fp_tetrapyr_na, "\nREMARK     2  No Such Site Exists  whose Size is Smaller than the Standard Site");
	}

	if (same_protein_name.size()!=0){
		fprintf(fp_tetrapyr_na, "\nREMARK     2");
		fprintf(fp_tetrapyr_na, "\nREMARK     2  No. of Tetragonal Pyramidal Sites whose Size is Same as that of the Standard Site: %d", _same_size);
		fprintf(fp_tetrapyr_na, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<same_protein_name.size(); p++){
			fprintf(fp_tetrapyr_na, "%s ", same_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_same_mol_type, "Same", "REMARK     2", fp_tetrapyr_na);
	}else{
		fprintf(fp_tetrapyr_na, "\nREMARK     2  No Such Site Exists  whose Size is Same as that of the Standard Site");
	}

	if (irregular_protein_name.size()!=0){
		fprintf(fp_tetrapyr_na, "\nREMARK     2");
		fprintf(fp_tetrapyr_na, "\nREMARK     2  No. of Tetragonal Pyramidal Sites whose Size is irregular in comparison to the Standard Site: %d", _irregular);
		fprintf(fp_tetrapyr_na, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<irregular_protein_name.size(); p++){
			fprintf(fp_tetrapyr_na, "%s ", irregular_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_irregular_mol_type, "Irregular", "REMARK     2", fp_tetrapyr_na);
	}else{
		fprintf(fp_tetrapyr_na, "\nREMARK     2  No Such Site Exists  whose Size is Irregular in comparison to the Standard Site");
	}

	fprintf(fp_tetrapyr_na, "\nREMARK     3");
	fprintf(fp_tetrapyr_na, "\nREMARK     3  Average Deviation of Distances of Atom-to-Atom:                         ");
	for (int p=0; p<_no_of_sides; p++){
		fprintf(fp_tetrapyr_na, "\nSide %d: %8.3lf", (p+1), _avg_devt_dist_vert_vert[p]);
	}

	fprintf(fp_tetrapyr_na, "\nREMARK     3         Average MSE Distance AA: %8.3lf", _avg_mse[0]);

	fprintf(fp_tetrapyr_na, "\nREMARK     3");
	fprintf(fp_tetrapyr_na, "\nREMARK     3  Average  Deviation of Angles of Made by Two Atoms with Another Atom:                         ");
	for (int p=0; p<_no_of_vertex_vertex_angles; p++){
		fprintf(fp_tetrapyr_na, "\nAngle %d: %8.3lf", (p+1), _avg_devt_angle_vertex_vertex[p]);
	}

	fprintf(fp_tetrapyr_na, "\nREMARK     3         Average MSE Angle AA: %8.3lf", _avg_mse[1]);

	fclose(fp_tetrapyr_na);
}

void tetrapyrStructNa::sort_and_count_max_sites(vector<string> mol_type, string qualifier, string remark, FILE* file_name){
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
