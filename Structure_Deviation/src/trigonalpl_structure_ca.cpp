#include "trigonalpl_structure_ca.h"

triplMetalicStructure* triplStructCa::_tripl_struct_standard_ca=new triplMetalicStructure[2];	//To store Standard, Standard_transformed
int triplStructCa::_no_of_binding_atoms=3;
int triplStructCa::_no_of_sides=3;
int triplStructCa::_no_of_angles_with_circum_centre=3;
int triplStructCa::_no_of_vertex_vertex_angles=3;
int triplStructCa::_no_of_features_to_compare=4;
double triplStructCa::_avg_dist_devt_from_cc=0.0;
double triplStructCa::_avg_devt_dist_vert_vert[]={0.0, 0.0, 0.0};
double triplStructCa::_avg_devt_angle_with_cc[]={0.0, 0.0, 0.0};
double triplStructCa::_avg_devt_angle_vertex_vertex[]={0.0, 0.0, 0.0};
double triplStructCa::_avg_mse[]={0.0, 0.0, 0.0, 0.0};
int triplStructCa::_nearer=0;
vector<string> triplStructCa::_nearer_mol_type;
vector<string> triplStructCa::nearer_protein_name;
int triplStructCa::_far_apart=0;
vector<string> triplStructCa::_far_mol_type;
vector<string> triplStructCa::far_apart_protein_name;
int triplStructCa::_larger_size=0;
vector<string> triplStructCa::_larger_mol_type;
vector<string> triplStructCa::larger_protein_name;
int triplStructCa::_smaller_size=0;
vector<string> triplStructCa::_smaller_mol_type;
vector<string> triplStructCa::smaller_protein_name;
int triplStructCa::_same_size=0;
vector<string> triplStructCa::_same_mol_type;
vector<string> triplStructCa::same_protein_name;

triplStructCa::triplStructCa(metal_binding_site met_site, bool first_time_tripl_ca, string accn){
	_accn=accn;
	_tripl_site=met_site;

	assert(_no_of_binding_atoms==3 && "No. of Binding Atoms should be 3!!");
	_tripl_struct_comparable_ca=new triplMetalicStructure[2];	//Comparable and Comparable_Transformed Trigonal_Planar_Metal Structure

	assert(_no_of_sides==3 && "No. of Sides should be 3!!");
	assert(_no_of_angles_with_circum_centre==3 && "For Trigonal Planar No. of Angles of Two Vertices with Circumcentre should be 3!!");
	assert(_no_of_vertex_vertex_angles==3 && "For Trigonal Planar No. of Angles of Two Vertices with Another Vertex should be 3!!");
	assert(_no_of_features_to_compare==4 && "For Trigonal Planar No. of Features to Compare should be 4!!");

	if (first_time_tripl_ca==true){
		cout<<"\n\nForm Standard Triangle for CA:";
		cout<<"\n\nStandard Triangle:\n";
		cout<<"\n\nOriginal Standard Atoms:";
		_tripl_struct_standard_ca[0].form_standard_triangle_CA();
		string metal_name=_tripl_site.get_metal_atom().get_atom_symbol();
		string pdb_out_standard_name="standard_original_trigonal_planar_"+metal_name+".stdn";
		_tripl_struct_standard_ca[0].write_atoms_pdb_format(pdb_out_standard_name);	//Generate PDB File for Standard Atoms

		_tripl_struct_standard_ca[1]=_tripl_struct_standard_ca[0].apply_transformations();

		cout<<"\n\nStandard Transformed Atoms (Circumcentre shifted to origin):\n";
		pdb_out_standard_name="standard_transformed_trigonal_planar_"+metal_name+".stdn";
		_tripl_struct_standard_ca[1].write_atoms_pdb_format(pdb_out_standard_name);

		string tex_out_name="standard_transformed_trigonal_planar_"+metal_name;
		_tripl_struct_standard_ca[1].write_atoms_single_tex_format(tex_out_name);

		_tripl_struct_standard_ca[1].measure_direction_cosines();
	}
}

triplMetalicStructure* triplStructCa::form_tripl_site(int total_tripl_sites){
	string metal_name=_tripl_site.get_metal_atom().get_atom_symbol();
	_site_type_comparable=_tripl_site.get_site_type();

	cout<<"Form Comparable Trigonal Planar for "<<metal_name<<":";
	cout<<"\nComparable Triangle:\n";
	cout<<"\n\nOriginal Comparable Atoms:";

	_tripl_struct_comparable_ca[0]=triplMetalicStructure(_tripl_site.get_inner_coord_sphere(), _tripl_site.get_metal_atom());
	string pdb_out_comparable_name=_accn+"_comparable_original_trigonal_planar_"+metal_name+"_"+_tripl_site.get_metal_atom().get_chain_name()+to_string(_tripl_site.get_metal_atom().get_residue_no())+".cmp";
	_tripl_struct_comparable_ca[0].write_atoms_pdb_format(pdb_out_comparable_name);	//Generate PDB File for Comparable Atoms

	_tripl_struct_comparable_ca[1]=_tripl_struct_comparable_ca[0].apply_transformations();

	cout<<"\n\nComparable Transformed Atoms (Circumcentre shifted to origin):\n";
	pdb_out_comparable_name=_accn+"_comparable_transformed_trigonal_planar_"+metal_name+"_"+_tripl_site.get_metal_atom().get_chain_name()+to_string(_tripl_site.get_metal_atom().get_residue_no())+".cmp";
	_tripl_struct_comparable_ca[1].write_atoms_pdb_format(pdb_out_comparable_name);

	string tex_out_name=_accn+"_comparable_transformed_trigonal_planar_"+metal_name+"_"+_tripl_site.get_metal_atom().get_chain_name()+to_string(_tripl_site.get_metal_atom().get_residue_no());
	_tripl_struct_comparable_ca[1].write_atoms_single_tex_format(tex_out_name);

	tex_out_name=_accn+"_comparable_transformed_combined_tripl_"+metal_name+"_"+_tripl_site.get_metal_atom().get_chain_name()+to_string(_tripl_site.get_metal_atom().get_residue_no());
	_tripl_struct_standard_ca[1].write_atoms_combined_tex_format(_tripl_struct_comparable_ca[1], tex_out_name);

	_tripl_struct_standard_ca[1].compare(_tripl_struct_comparable_ca[1]);

	bool check_nearer_accn=false, check_far_accn=false, check_larger_accn=false, check_smaller_accn=false, check_same_accn=false;
	if (_tripl_struct_standard_ca[1].get_nearer()==true){
		_nearer++;
		_nearer_mol_type.push_back(_tripl_site.get_site_type());

		if (check_nearer_accn==false){
			int p;
			for (p=0; p<nearer_protein_name.size(); p++){
				if (_accn==nearer_protein_name.at(p)){
					break;
				}
			}

			if (p==nearer_protein_name.size()){
				nearer_protein_name.push_back(_accn);
			}

			check_nearer_accn=true;
		}
	}else{
		_far_apart++;
		_far_mol_type.push_back(_tripl_site.get_site_type());

		if (check_far_accn==false){
			int p;
			for (p=0; p<far_apart_protein_name.size(); p++){
				if (_accn==far_apart_protein_name.at(p)){
					break;
				}
			}

			if (p==far_apart_protein_name.size()){
				far_apart_protein_name.push_back(_accn);
			}

			check_far_accn=true;
		}
	}

	if (_tripl_struct_standard_ca[1].get_structure_size()=="larger"){
		_larger_size++;
		_larger_mol_type.push_back(_tripl_site.get_site_type());

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
	}else if (_tripl_struct_standard_ca[1].get_structure_size()=="smaller"){
		_smaller_size++;
		_smaller_mol_type.push_back(_tripl_site.get_site_type());

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
	}else if (_tripl_struct_standard_ca[1].get_structure_size()=="same"){
		_same_size++;
		_same_mol_type.push_back(_tripl_site.get_site_type());

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
	}

	compute_average_deviations_percentages_in_all_tripl_molecules(total_tripl_sites);

	string out_file_name="Deviation_Report_Trigonal_Planar_"+_accn+"_"+metal_name+"_"+_tripl_site.get_metal_atom().get_chain_name()+to_string(_tripl_site.get_metal_atom().get_residue_no())+".mbsi";
	FILE* fp1=fopen(out_file_name.c_str(), "w");
	cout<<"\n\nStarts reports generation for "<<metal_name<<" Ions:\n";
	cout <<"\nOutput Report File Name: " <<out_file_name<< endl;

	string site_type_standard=_tripl_struct_standard_ca[0].assign_standard_site_type();
	_tripl_struct_standard_ca[1].gen_report_metal_binding_sites(_tripl_struct_standard_ca[0], _tripl_struct_comparable_ca[0], _tripl_struct_comparable_ca[1], metal_name, site_type_standard, _tripl_site.get_site_type(), fp1);
	fclose(fp1);

	pdb_out_comparable_name=_accn+"_comparable_different_trigonal_planar_"+metal_name+"_"+_tripl_site.get_metal_atom().get_chain_name()+to_string(_tripl_site.get_metal_atom().get_residue_no())+".comp";
	_tripl_struct_standard_ca[0].comparison(_tripl_struct_comparable_ca[0], pdb_out_comparable_name);

	return _tripl_struct_comparable_ca;
}

void triplStructCa::compute_average_deviations_percentages_in_all_tripl_molecules(int total_tripl_sites){
	double* temp_side_dev=_tripl_struct_standard_ca[1].get_side_deviations();
	double* temp_angle_cc=_tripl_struct_standard_ca[1].get_angle_deviations_with_cc();
	double* temp_angle_vv=_tripl_struct_standard_ca[1].get_angle_deviations_vertex_vertex();
	double* temp_mse=_tripl_struct_standard_ca[1].get_mse();
	string metal_name=_tripl_site.get_metal_atom().get_atom_symbol();

	if (total_tripl_sites==1){
		_avg_dist_devt_from_cc=_tripl_struct_standard_ca[1].get_distance_deviations_from_cc();

		for (int p=0; p<_no_of_sides; p++){
			_avg_devt_dist_vert_vert[p]=temp_side_dev[p];
		}

		for (int p=0; p<_no_of_angles_with_circum_centre; p++){
			_avg_devt_angle_with_cc[p]=temp_angle_cc[p];
		}

		for (int p=0; p<_no_of_vertex_vertex_angles; p++){
			cout<<"\nIN "<<metal_name<<" SITE:"<<temp_angle_vv[p];
			_avg_devt_angle_vertex_vertex[p]=temp_angle_vv[p];
		}

		for (int p=0; p<_no_of_features_to_compare; p++){
			_avg_mse[p]=temp_mse[p];
		}
	}else{
		_avg_dist_devt_from_cc=(_avg_dist_devt_from_cc+_tripl_struct_standard_ca[1].get_distance_deviations_from_cc())/2.0;

		for (int p=0; p<_no_of_sides; p++){
			_avg_devt_dist_vert_vert[p]=(_avg_devt_dist_vert_vert[p]+temp_side_dev[p])/2.0;
		}

		for (int p=0; p<_no_of_angles_with_circum_centre; p++){
			_avg_devt_angle_with_cc[p]=(_avg_devt_angle_with_cc[p]+temp_angle_cc[p])/2.0;
		}

		for (int p=0; p<_no_of_vertex_vertex_angles; p++){
			cout<<"\nIN "<<metal_name<<" SITE:"<<temp_angle_vv[p];
			_avg_devt_angle_vertex_vertex[p]=(_avg_devt_angle_vertex_vertex[p]+temp_angle_vv[p])/2.0;
		}

		for (int p=0; p<_no_of_features_to_compare; p++){
			_avg_mse[p]=(_avg_mse[p]+temp_mse[p])/2.0;
		}
	}
}

void triplStructCa::gen_summary_report_tripl_ca(int total_tripl_sites){
	string summary_file_tripl_na="Summary_Deviation_Report_CA_trigonal_planar.mbsi";
	FILE *fp_tripl_ca=fopen(summary_file_tripl_na.c_str(), "w");

	fprintf(fp_tripl_ca, "HEADER     PROGRAM NAME: StructureDeviation VERSION: 1.1");
	fprintf(fp_tripl_ca, "\nTITLE      Summary Report for Computation of Structural Deviations of Two Trigonal Planar CA Binding Sites with Coordination Number 3");

	fprintf(fp_tripl_ca, "\nREMARK     1");
	fprintf(fp_tripl_ca, "\nREMARK     1  No. of Trigonal Planar Sites Compared with Standard CA ION: %d", total_tripl_sites);

	if (nearer_protein_name.size()!=0){
		fprintf(fp_tripl_ca, "\nREMARK     2");
		fprintf(fp_tripl_ca, "\nREMARK     2  No. of Trigonal Planar Sites whose Binding Atoms are Nearer to Binding Atoms of Standard Site: %d", _nearer);
		fprintf(fp_tripl_ca, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<nearer_protein_name.size(); p++){
			fprintf(fp_tripl_ca, "%s ", nearer_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_nearer_mol_type, "Nearer", "REMARK     2", fp_tripl_ca);
	}else{
		fprintf(fp_tripl_ca, "\nREMARK     2  No Such Site Exists whose Binding Atoms are Nearer to Binding Atoms of Standard Site");
	}

	if (far_apart_protein_name.size()!=0){
		fprintf(fp_tripl_ca, "\nREMARK     2");
		fprintf(fp_tripl_ca, "\nREMARK     2  No. of Trigonal Planar Sites whose Binding Atoms are Far Apart as Compared to Standard Site: %d", _far_apart);
		fprintf(fp_tripl_ca, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<far_apart_protein_name.size(); p++){
			fprintf(fp_tripl_ca, "%s ", far_apart_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_far_mol_type, "Far_Apart", "REMARK     2", fp_tripl_ca);
	}else{
		fprintf(fp_tripl_ca, "\nREMARK     2  No Such Site Exists whose Binding Atoms are Far Apart as Compared to Standard Site");
	}

	if (larger_protein_name.size()!=0){
		fprintf(fp_tripl_ca, "\nREMARK     2");
		fprintf(fp_tripl_ca, "\nREMARK     2  No. of Trigonal Planar Sites whose Size is Larger than the Standard Site: %d", _larger_size);
		fprintf(fp_tripl_ca, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<larger_protein_name.size(); p++){
			fprintf(fp_tripl_ca, "%s ", larger_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_larger_mol_type, "Larger", "REMARK     2", fp_tripl_ca);
	}else{
		fprintf(fp_tripl_ca, "\nREMARK     2  No Such Site Exists  whose Size is Larger than the Standard Site");
	}

	if (smaller_protein_name.size()!=0){
		fprintf(fp_tripl_ca, "\nREMARK     2");
		fprintf(fp_tripl_ca, "\nREMARK     2  No. of Trigonal Planar Sites whose Size is Smaller than the Standard Site: %d", _smaller_size);
		fprintf(fp_tripl_ca, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<smaller_protein_name.size(); p++){
			fprintf(fp_tripl_ca, "%s ", smaller_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_smaller_mol_type, "Smaller", "REMARK     2", fp_tripl_ca);
	}else{
		fprintf(fp_tripl_ca, "\nREMARK     2  No Such Site Exists  whose Size is Smaller than the Standard Site");
	}

	if (same_protein_name.size()!=0){
		fprintf(fp_tripl_ca, "\nREMARK     2");
		fprintf(fp_tripl_ca, "\nREMARK     2  No. of Trigonal Planar Sites whose Size is Same as that of the Standard Site: %d", _same_size);
		fprintf(fp_tripl_ca, "\nREMARK     2  Name of Proteins: ");
		for (int p=0; p<same_protein_name.size(); p++){
			fprintf(fp_tripl_ca, "%s ", same_protein_name.at(p).c_str());
		}

		sort_and_count_max_sites(_same_mol_type, "Same", "REMARK     2", fp_tripl_ca);
	}else{
		fprintf(fp_tripl_ca, "\nREMARK     2  No Such Site Exists  whose Size is Same as that of the Standard Site");
	}

	fprintf(fp_tripl_ca, "\nREMARK     3");
	fprintf(fp_tripl_ca, "\nREMARK     3  Average Deviation of Distance from Circumcentre to Any Atom: %8.3lf", _avg_dist_devt_from_cc);
	fprintf(fp_tripl_ca, "\nREMARK     3         Average MSE Distance CC: %8.3lf", _avg_mse[0]);

	fprintf(fp_tripl_ca, "\nREMARK     3");
	fprintf(fp_tripl_ca, "\nREMARK     3  Average Deviation of Distances of Atom-to-Atom:                         ");
	for (int p=0; p<_no_of_sides; p++){
		fprintf(fp_tripl_ca, "\nSide %d: %8.3lf", (p+1), _avg_devt_dist_vert_vert[p]);
	}

	fprintf(fp_tripl_ca, "\nREMARK     3         Average MSE Distance AA: %8.3lf", _avg_mse[1]);

	fprintf(fp_tripl_ca, "\nREMARK     3");
	fprintf(fp_tripl_ca, "\nREMARK     3  Average  Deviation of Angles Made by Two Atoms with Circumcentre:                         ");
	for (int p=0; p<_no_of_angles_with_circum_centre; p++){
		fprintf(fp_tripl_ca, "\nAngle %d: %8.3lf", (p+1), _avg_devt_angle_with_cc[p]);
	}

	fprintf(fp_tripl_ca, "\nREMARK     3         Average MSE Angle CC: %8.3lf", _avg_mse[2]);

	fprintf(fp_tripl_ca, "\nREMARK     3");
	fprintf(fp_tripl_ca, "\nREMARK     3  Average  Deviation of Angles of Made by Two Atoms with Another Atom:                         ");
	for (int p=0; p<_no_of_vertex_vertex_angles; p++){
		fprintf(fp_tripl_ca, "\nAngle %d: %8.3lf", (p+1), _avg_devt_angle_vertex_vertex[p]);
	}

	fprintf(fp_tripl_ca, "\nREMARK     3         Average MSE Angle AA: %8.3lf", _avg_mse[3]);

	fclose(fp_tripl_ca);
}

void triplStructCa::sort_and_count_max_sites(vector<string> mol_type, string qualifier, string remark, FILE* file_name){
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
