#include "biomolecule.h"
#include <set>
#include <fstream>

protein::protein(string pname, string mol_file, int metal_threshold){
	int p1 = pname.find_last_of(".");
	string _accn = pname.substr(0, p1);
	string _ext = pname.substr(p1);

	bool read_ok=false;
    //_ext is extension of input file
    if (_ext==".cif"){
    	read_ok=read_cif_atoms(pname, mol_file);
    }else if (_ext==".pdb"){
    	read_ok=read_pdb_atoms(pname, mol_file);
    }

    if (read_ok==true){
    	bool comp_ok=compute_binding_sites(_accn, _ext, metal_threshold);

    	if (comp_ok==true){
    		generate_general_report(_accn);
    	}else{
    		cout<<"Can't Compute Binding Site!!";
    		exit(1);
    	}
    }else{
    	cout<<"Can't Read Mol2 File!!";
    	exit(1);
    }
}

vector<atom> protein::get_na_atoms(){
	return _na_atoms;
}

vector<atom> protein::get_ca_atoms(){
	return _ca_atoms;
}

vector<atom> protein::get_k_atoms(){
	return _k_atoms;
}

vector<atom> protein::get_mg_atoms(){
	return _mg_atoms;
}

vector<atom> protein::get_mn_atoms(){
	return _mn_atoms;
}

vector<atom> protein::get_fe_atoms(){
	return _fe_atoms;
}

vector<atom> protein::get_cu_atoms(){
	return _cu_atoms;
}

vector<atom> protein::get_zn_atoms(){
	return _zn_atoms;
}

vector<atom> protein::get_co_atoms(){
	return _co_atoms;
}

vector<atom> protein::get_mo_atoms(){
	return _mo_atoms;
}

vector<atom> protein::get_ag_atoms(){
	return _ag_atoms;
}

vector<atom> protein::get_pb_atoms(){
	return _pb_atoms;
}

vector<atom> protein::get_hg_atoms(){
	return _hg_atoms;
}

vector<atom> protein::get_as_atoms(){
	return _as_atoms;
}

vector<atom> protein::get_cd_atoms(){
	return _cd_atoms;
}

vector<atom> protein::get_ni_atoms(){
	return _ni_atoms;
}

vector<atom> protein::get_cr_atoms(){
	return _cr_atoms;
}

vector<atom> protein::get_cs_atoms(){
	return _cs_atoms;
}

vector<string> protein::get_metal_list(){
	return _metal_list;
}

vector<string> protein::get_unique_metal_list(){
	return _unique_metal_list;
}

vector<int> protein::get_count_unique_metal(){
	return _count_unique_metal;
}

vector<atom> protein::get_non_metal_atoms(){
	return _non_metals;
}

vector<metal_binding_site> protein::get_na_binding_sites(){
	return _na_sites;
}

vector<metal_binding_site> protein::get_ca_binding_sites(){
	return _ca_sites;
}

vector<metal_binding_site> protein::get_k_binding_sites(){
	return _k_sites;
}

vector<metal_binding_site> protein::get_mg_binding_sites(){
	return _mg_sites;
}

vector<metal_binding_site> protein::get_mn_binding_sites(){
	return _mn_sites;
}

vector<metal_binding_site> protein::get_fe_binding_sites(){
	return _fe_sites;
}

vector<metal_binding_site> protein::get_cu_binding_sites(){
	return _cu_sites;
}

vector<metal_binding_site> protein::get_zn_binding_sites(){
	return _zn_sites;
}

vector<metal_binding_site> protein::get_co_binding_sites(){
	return _co_sites;
}

vector<metal_binding_site> protein::get_mo_binding_sites(){
	return _mo_sites;
}

vector<metal_binding_site> protein::get_ag_binding_sites(){
	return _ag_sites;
}

vector<metal_binding_site> protein::get_pb_binding_sites(){
	return _pb_sites;
}

vector<metal_binding_site> protein::get_hg_binding_sites(){
	return _hg_sites;
}

vector<metal_binding_site> protein::get_as_binding_sites(){
	return _as_sites;
}

vector<metal_binding_site> protein::get_cd_binding_sites(){
	return _cd_sites;
}

vector<metal_binding_site> protein::get_ni_binding_sites(){
	return _ni_sites;
}

vector<metal_binding_site> protein::get_cr_binding_sites(){
	return _cr_sites;
}

vector<metal_binding_site> protein::get_cs_binding_sites(){
	return _cs_sites;
}

vector<bond> protein::get_bond_info(){
	return _bond_info;
}

vector<vector<atom> > protein::get_all_metal_atoms(){
	vector<vector<atom> > _all_metal_atoms;

	_all_metal_atoms.push_back(_na_atoms);
	_all_metal_atoms.push_back(_ca_atoms);
	_all_metal_atoms.push_back(_k_atoms);
	_all_metal_atoms.push_back(_mg_atoms);
	_all_metal_atoms.push_back(_mn_atoms);
	_all_metal_atoms.push_back(_fe_atoms);
	_all_metal_atoms.push_back(_cu_atoms);
	_all_metal_atoms.push_back(_zn_atoms);
	_all_metal_atoms.push_back(_co_atoms);
	_all_metal_atoms.push_back(_mo_atoms);
	_all_metal_atoms.push_back(_ag_atoms);
	_all_metal_atoms.push_back(_pb_atoms);
	_all_metal_atoms.push_back(_hg_atoms);
	_all_metal_atoms.push_back(_as_atoms);
	_all_metal_atoms.push_back(_cd_atoms);
	_all_metal_atoms.push_back(_ni_atoms);
	_all_metal_atoms.push_back(_cr_atoms);
	_all_metal_atoms.push_back(_cs_atoms);

	return _all_metal_atoms;
}

vector<vector<metal_binding_site> > protein::get_all_metal_sites(){
	vector<vector<metal_binding_site> > _all_metal_sites;

	if (_na_sites.size()!=0){
		_all_metal_sites.push_back(_na_sites);
	}

	if (_ca_sites.size()!=0){
		_all_metal_sites.push_back(_ca_sites);
	}

	if (_k_sites.size()!=0){
		_all_metal_sites.push_back(_k_sites);
	}

	if (_mg_sites.size()!=0){
		_all_metal_sites.push_back(_mg_sites);
	}

	if (_mn_sites.size()!=0){
		_all_metal_sites.push_back(_mn_sites);
	}

	if (_fe_sites.size()!=0){
		_all_metal_sites.push_back(_fe_sites);
	}

	if (_cu_sites.size()!=0){
		_all_metal_sites.push_back(_cu_sites);
	}

	if (_zn_sites.size()!=0){
		_all_metal_sites.push_back(_zn_sites);
	}

	if (_co_sites.size()!=0){
		_all_metal_sites.push_back(_co_sites);
	}

	if (_mo_sites.size()!=0){
		_all_metal_sites.push_back(_mo_sites);
	}

	if (_ag_sites.size()!=0){
		_all_metal_sites.push_back(_ag_sites);
	}

	if (_pb_sites.size()!=0){
		_all_metal_sites.push_back(_pb_sites);
	}

	if (_hg_sites.size()!=0){
		_all_metal_sites.push_back(_hg_sites);
	}

	if (_as_sites.size()!=0){
		_all_metal_sites.push_back(_as_sites);
	}

	if (_cd_sites.size()!=0){
		_all_metal_sites.push_back(_cd_sites);
	}

	if (_ni_sites.size()!=0){
		_all_metal_sites.push_back(_ni_sites);
	}

	if (_cr_sites.size()!=0){
		_all_metal_sites.push_back(_cr_sites);
	}

	if (_cs_sites.size()!=0){
		_all_metal_sites.push_back(_cs_sites);
	}

	return _all_metal_sites;
}

vector<vector<string> > protein::read_mol_atom_info(string mol_file){
	vector<vector<string> > mol_atom_info;
	ifstream fm(mol_file.c_str());

	if (fm.is_open()){
		string l;

		while (getline(fm, l)){
			if (l!="@<TRIPOS>ATOM"){
				continue;
			}else{
				break;
			}
		}

		while (getline(fm, l)){
			if (l!="@<TRIPOS>BOND"){
				vector<string> result;
				istringstream iss(l);

				for(string s; iss >> s; ){
					result.push_back(s);
				}

				//if (result.size()==9){
					mol_atom_info.push_back(result);
				//}else{
					//cout<<"\nMol Atom Size: "<<result.size();
					//exit(1);
				//}
			}else{
				break;
			}
		}
	}

	fm.close();
cout<<"\n\nNo. of MOL Atoms: "<<mol_atom_info.size();
	return mol_atom_info;
}

vector<vector<string> > protein::read_mol_bond_info(string mol_file){
	vector<vector<string> > mol_bond_info;
	ifstream fm(mol_file.c_str());

	if (fm.is_open()){
		string l;

		while (getline(fm, l)){
			if (l!="@<TRIPOS>BOND"){
				continue;
			}else{
				break;
			}
		}

		while (getline(fm, l)){
			vector<string> result;
			istringstream iss(l);

			for(string s; iss >> s; ){
				result.push_back(s);
			}

			result.push_back("false");	//Update Status of Bonded Atom1
			result.push_back("false");	//Update Status of Bonded Atom2

			//if (result.size()==6){
				mol_bond_info.push_back(result);
			//}else{
				//cout<<"\nMol Bond Size: "<<result.size();
				//exit(1);
			//}
		}
	}

	fm.close();

	return mol_bond_info;
}

bool protein::read_cif_atoms(string prot_file, string mol_file){
	cif::Document _cif_document = cif::read_file(prot_file);
	const string _tag_atom="_atom_site.";
	string feature_atom[]={"group_PDB", "id", "auth_atom_id", "auth_comp_id", "auth_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "type_symbol", "label_alt_id"};
	vector<string> _feature_atom_vect;
	vector<atom> atom_group;

	vector<vector<string> > mol_atom_info=read_mol_atom_info(mol_file);
	vector<vector<string> > mol_bond_info=read_mol_bond_info(mol_file);

	if (mol_atom_info.size()!=0 && mol_bond_info.size()!=0){
		for (auto ft: feature_atom){
			_feature_atom_vect.push_back(ft);
		}

		int i=0;
		for (cif::Block& block : _cif_document.blocks) {
			for (auto cc: block.find(_tag_atom, _feature_atom_vect)) {
				char* end;

				if((cc.at(0) == "ATOM") || (cc.at(0) == "HETATM")){
					atom_group.push_back(atom(cc.at(0), strtol(cc.at(1).c_str(),&end,10), cc.at(2), cc.at(12), cc.at(3), cc.at(4), strtol(cc.at(5).c_str(),&end,10), strtod(cc.at(6).c_str(),&end), strtod(cc.at(7).c_str(),&end), strtod(cc.at(8).c_str(),&end), strtod(cc.at(9).c_str(),&end), strtod(cc.at(10).c_str(),&end), 0.0, cc.at(11), strtod(mol_atom_info.at(i).at(8).c_str(),&end)));

					for (int j=0; j<mol_bond_info.size(); j++){
						if ((mol_atom_info.at(i).at(0)==mol_bond_info.at(j).at(1)) && (mol_bond_info.at(j).at(4)=="false")){
							long int diff=strtol(cc.at(1).c_str(),&end,10)-strtol(mol_atom_info.at(i).at(0).c_str(), &end, 10);
							long int att=strtol(mol_bond_info.at(j).at(1).c_str(), &end, 10)+diff;
							mol_bond_info.at(j).at(1)=to_string(att);
							mol_bond_info.at(j).at(4)="true";
						}else if ((mol_atom_info.at(i).at(0)==mol_bond_info.at(j).at(2)) && (mol_bond_info.at(j).at(5)=="false")){
							long int diff=strtol(cc.at(1).c_str(),&end,10)-strtol(mol_atom_info.at(i).at(0).c_str(), &end, 10);
							long int att=strtol(mol_bond_info.at(j).at(2).c_str(), &end, 10)+diff;
							mol_bond_info.at(j).at(2)=to_string(att);
							mol_bond_info.at(j).at(5)="true";
						}
					}

					i++;
				}
			}
		}

		cout<<"\nBOND SIZE: "<<mol_bond_info.size();
		cout<<"\nINPUT: "<<prot_file;
		assign_bond_info(atom_group, mol_bond_info);
		make_residue_metal_non_metal_group(atom_group);

		set<string> s(_metal_list.begin(), _metal_list.end());
		for (auto e: s){
			_unique_metal_list.push_back(e);
			int c=count(_metal_list.begin(), _metal_list.end(), e);
			_count_unique_metal.push_back(c);
		}

		cout<<"\nSuccessfully read "<<prot_file<<" file!";
		cout<<"\n\nNo. of Atoms: "<<atom_group.size();
		//for (int i=0; i<atom_group.size(); i++){
    		//atom_group.at(i).fprint(stdout);
		//}

		return true;
	}else{
		return false;
	}
}

void protein::assign_bond_info(vector<atom> atom_group, vector<vector<string> > mol_bond_info){
	atom atm1;
	atom atm2;
	char *end;

	for (auto bnd: mol_bond_info){
		for (auto atm:atom_group){
			if (strtol(bnd.at(1).c_str(), &end, 10)==atm.get_serial_no()){
				atm1=atm;
				break;
			}
		}

		for (auto atm:atom_group){
			if (strtol(bnd.at(2).c_str(), &end, 10)==atm.get_serial_no()){
				atm2=atm;
				break;
			}
		}

		_bond_info.push_back(bond(atm1, atm2, bnd.at(3), false, false));
	}

	//cout<<"\n\nBond Information: "<<_bond_info.size()<<endl;
	//for (int i=0; i<_bond_info.size(); i++){
		//cout<<_bond_info.at(i).get_bonded_atom1().get_serial_no()<<", "<<_bond_info.at(i).get_bonded_atom2().get_serial_no()<<", "<<_bond_info.at(i).get_bond_type()<<", "<<_bond_info.at(i).get_bond_strength()<<endl;
	//}
}

void protein::make_residue_metal_non_metal_group(vector<atom> atom_group){
	vector<atom> sub_group;

	for (int i=1; i<atom_group.size(); i++){
		if (sub_group.size()==0){
			sub_group.push_back(atom_group.at(i-1));
		}

	    long int prev_res_no=atom_group.at(i-1).get_residue_no();
	    string prev_chain_name=atom_group.at(i-1).get_chain_name();

	    if ((prev_chain_name==atom_group.at(i).get_chain_name()) && (prev_res_no==atom_group.at(i).get_residue_no())){
	  		sub_group.push_back(atom_group.at(i));
	    }else{
	    	_residue_group.push_back(residue(sub_group));
	    	separate_metal_and_non_metals(sub_group);
	    	sub_group.clear();

	    	if (i==(atom_group.size()-1)){
	    		sub_group.push_back(atom_group.at(i));
	    		_residue_group.push_back(residue(sub_group));
	    		separate_metal_and_non_metals(sub_group);
	    	}
	    }
	}
}

void protein::separate_metal_and_non_metals(vector<atom> sub_group){
	char *end;
	for (auto atm: sub_group){
		if (atm.get_atom_type()=="HETATM" && ((atm.get_alt_conf_indicator()==" ") || (atm.get_alt_conf_indicator()=="A" || strtol(atm.get_alt_conf_indicator().c_str(), &end, 10)==1))){
			if (atm.get_atom_symbol()=="NA"){
				_na_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="CA"){
				_ca_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="K"){
				_k_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="MG"){
				_mg_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="MN"){
				_mn_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="FE"){
				_fe_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="CU"){
				_cu_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="ZN"){
				_zn_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="CO"){
				_co_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="MO"){
				_mo_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="AG"){
				_ag_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="PB"){
				_pb_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="HG"){
				_hg_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="AS"){
				_as_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="CD"){
				_cd_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="NI"){
				_ni_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="CR"){
				_cr_atoms.push_back(atm);
				_metal_list.push_back(atm.get_atom_symbol());
			}else if (atm.get_atom_symbol()=="CS"){
				_cs_atoms.push_back(atm);
		    	_metal_list.push_back(atm.get_atom_symbol());
			}else{
				_non_metals.push_back(atm);
			}
		}
	}
}

int protein::set_pdb_atom_field_index(vector<vector<int> > & pdb_atom_field_index){
	int no_of_features=max_atom-min_atom+1;

	for (int i=atom_type; i<=atom_symbol; i++){
		vector<int> sub_field;

		//1st col - Field No., 2nd col - start index, 3rd - end index
		switch(i){
			case atom_type: sub_field.push_back(i);
		                    sub_field.push_back(1);
		                    sub_field.push_back(6);
		                    assert((sub_field.at(1)==1 && sub_field.at(2)==6) && "Invalid column nos. for atom type");
		                    pdb_atom_field_index.push_back(sub_field);
		                    break;

			case serial_no: sub_field.push_back(i);
			                sub_field.push_back(7);
			                sub_field.push_back(11);
			                assert((sub_field.at(1)==7 && sub_field.at(2)==11) && "Invalid column nos. for serial no.");
			                pdb_atom_field_index.push_back(sub_field);
			                break;

			case atom_name: sub_field.push_back(i);
						    sub_field.push_back(13);
						    sub_field.push_back(16);
						    assert((sub_field.at(1)==13 && sub_field.at(2)==16) && "Invalid column nos. for atom name");
						    pdb_atom_field_index.push_back(sub_field);
						    break;

			case alt_conf_indicator: sub_field.push_back(i);
									 sub_field.push_back(17);
									 sub_field.push_back(17);
									 assert((sub_field.at(1)==17 && sub_field.at(2)==17) && "Invalid column nos. for Alternate Conformation Indicator");
									 pdb_atom_field_index.push_back(sub_field);
									 break;

			case residue_name: sub_field.push_back(i);
						       sub_field.push_back(18);
						       sub_field.push_back(20);
						       assert((sub_field.at(1)==18 && sub_field.at(2)==20) && "Invalid column nos. for residue name");
						       pdb_atom_field_index.push_back(sub_field);
						       break;

			case chain_id: sub_field.push_back(i);
						   sub_field.push_back(22);
						   sub_field.push_back(22);
						   assert((sub_field.at(1)==22 && sub_field.at(2)==22) && "Invalid column nos. for chain id");
						   pdb_atom_field_index.push_back(sub_field);
						   break;

			case residue_no: sub_field.push_back(i);
		                     sub_field.push_back(23);
		                     sub_field.push_back(26);
		                     assert((sub_field.at(1)==23 && sub_field.at(2)==26) && "Invalid column nos. for residue no.");
		                     pdb_atom_field_index.push_back(sub_field);
		                     break;

			case x_coord: sub_field.push_back(i);
	                      sub_field.push_back(31);
	                      sub_field.push_back(38);
	                      assert((sub_field.at(1)==31 && sub_field.at(2)==38) && "Invalid column nos. for x-coord");
	                      pdb_atom_field_index.push_back(sub_field);
	                      break;

			case y_coord: sub_field.push_back(i);
	                      sub_field.push_back(39);
	                      sub_field.push_back(46);
	                      assert((sub_field.at(1)==39 && sub_field.at(2)==46) && "Invalid column nos. for y-coord");
	                      pdb_atom_field_index.push_back(sub_field);
	                      break;

			case z_coord: sub_field.push_back(i);
	                      sub_field.push_back(47);
	                      sub_field.push_back(54);
	                      assert((sub_field.at(1)==47 && sub_field.at(2)==54) && "Invalid column nos. for z-coord");
	                      pdb_atom_field_index.push_back(sub_field);
	                      break;

			case occupancy: sub_field.push_back(i);
	                        sub_field.push_back(55);
	                        sub_field.push_back(60);
	                        assert((sub_field.at(1)==55 && sub_field.at(2)==60) && "Invalid column nos. for occupancy");
	                        pdb_atom_field_index.push_back(sub_field);
	                        break;

			case t_factor: sub_field.push_back(i);
	                       sub_field.push_back(61);
	                       sub_field.push_back(66);
	                       assert((sub_field.at(1)==61 && sub_field.at(2)==66) && "Invalid column nos. for temp. factor");
	                       pdb_atom_field_index.push_back(sub_field);
	                       break;

			case atom_symbol: sub_field.push_back(i);
	                          sub_field.push_back(77);
	                          sub_field.push_back(78);
	                          assert((sub_field.at(1)==77 && sub_field.at(2)==78) && "Invalid column nos. for atom symbol");
	                          pdb_atom_field_index.push_back(sub_field);
	                          break;

			default: cout<<"\nInvalid Features!!";
			         exit(1);
		}
	}

	return no_of_features;
}

bool protein::read_pdb_atoms(string prot_file, string mol_file){
	vector<string> pdb_atom_lines;
	vector<string> pdb_link_lines;
	vector<vector<int> > pdb_atom_field_index;
	vector<vector<int> > pdb_link_field_index;
	vector<atom> atom_group;
	vector<vector<string> > mol_atom_info=read_mol_atom_info(mol_file);
	vector<vector<string> > mol_bond_info=read_mol_bond_info(mol_file);

	if (mol_atom_info.size()!=0 && mol_bond_info.size()!=0){
		ifstream fp(prot_file.c_str());
		int _no_of_pdb_atom_features=set_pdb_atom_field_index(pdb_atom_field_index);

		if (fp.is_open()){
			string l;
			char *end;
			int i=0;

			while (getline(fp, l)){
				vector<string> str;
				string tag=l.substr(0, 6);
				assert((tag.length()==6) && "Invalid Record Tag!");
				tag=tag.substr(tag.find_first_not_of(" "));
				tag=tag.substr(0, tag.find_last_not_of(" ")+1);

				if ((tag=="ATOM") || (tag=="HETATM")){
					pdb_atom_lines.push_back(l);
					str.push_back(tag);
					split_pdb_atom_str(l, str, pdb_atom_field_index);

					if (str.size()==13){
						string st=mol_atom_info.at(i).at(mol_atom_info.at(i).size()-1);
						atom_group.push_back(atom(str.at(0), strtol(str.at(1).c_str(), &end,10), str.at(2), str.at(3), str.at(4), str.at(5), strtol(str.at(6).c_str(), &end,10), strtod(str.at(7).c_str(), &end), strtod(str.at(8).c_str(), &end), strtod(str.at(9).c_str(), &end), strtod(str.at(10).c_str(), &end), strtod(str.at(11).c_str(), &end), 0.0, str.at(12), strtod(st.c_str(), &end)));
					}else{
						cout<<"\nSTR: "<<str.size();
						exit(1);
					}

					for (int j=0; j<mol_bond_info.size(); j++){
						if ((mol_atom_info.at(i).at(0)==mol_bond_info.at(j).at(1)) && (mol_bond_info.at(j).at(4)=="false")){
							long int diff=strtol(str.at(1).c_str(),&end,10)-strtol(mol_atom_info.at(i).at(0).c_str(), &end, 10);
							long int att=strtol(mol_bond_info.at(j).at(1).c_str(), &end, 10)+diff;
							mol_bond_info.at(j).at(1)=to_string(att);
							mol_bond_info.at(j).at(4)="true";
						}else if ((mol_atom_info.at(i).at(0)==mol_bond_info.at(j).at(2)) && (mol_bond_info.at(j).at(5)=="false")){
							long int diff=strtol(str.at(1).c_str(),&end,10)-strtol(mol_atom_info.at(i).at(0).c_str(), &end, 10);
							long int att=strtol(mol_bond_info.at(j).at(2).c_str(), &end, 10)+diff;
							mol_bond_info.at(j).at(2)=to_string(att);
							mol_bond_info.at(j).at(5)="true";
						}
					}

					i++;
				}
			}

			fp.close();

			cout<<"\nBOND SIZE: "<<mol_bond_info.size();
			cout<<"\nINPUT: "<<prot_file;
	    	assign_bond_info(atom_group, mol_bond_info);
			make_residue_metal_non_metal_group(atom_group);

			set<string> s(_metal_list.begin(), _metal_list.end());
			for (auto e: s){
				_unique_metal_list.push_back(e);
				int c=count(_metal_list.begin(), _metal_list.end(), e);
				_count_unique_metal.push_back(c);
			}

			cout<<"\nSuccessfully read "<<prot_file<<" file!";
			cout<<"\n\nNo. of Atoms: "<<atom_group.size();
			//for (int i=0; i<atom_group.size(); i++){
	    		//atom_group.at(i).fprint(stdout);
	    	//}

			//cout<<"\n\nAfter Change - No. of MOL Bonds: "<<mol_bond_info.size();
			//for (auto bnd: mol_bond_info){
				//cout<<"\n"<<bnd.at(0)<<"  "<<bnd.at(1)<<"  "<<bnd.at(2)<<"  "<<bnd.at(3);
			//}
		}else{
			cout<<"\nUnable to open input file!";
			return false;
		}

		return true;
	}else{
		return false;
	}
}

void protein::split_pdb_atom_str(string l, vector<string> & str, vector<vector<int> > & pdb_field_index){
	cout<<endl;
	for(int i=serial_no; i<=atom_symbol; i++){
		int start=pdb_field_index.at(i).at(1);
		int end=pdb_field_index.at(i).at(2);
		string s=l.substr((start-1), (end-start+1));

		int pos=s.find_first_not_of(" ");
		if (pos!=string::npos){
			string st=s.substr(pos);

			pos=st.find_last_not_of(" ");

			if (pos!=string::npos){
				st=st.substr(0, pos+1);
				str.push_back(st);
			}
		}else{
			str.push_back(s);
		}

		cout<<endl;
	}
}

bool protein::compute_binding_sites(string accn, string ext, int metal_threshold){
	vector<vector<atom> > _all_metal_atoms=get_all_metal_atoms();

	if (_all_metal_atoms.size()==0){
		cout<<"\nThis Protein has No Binding Metal Ion";
		return false;
	}else{
		for (auto _metal_atoms: _all_metal_atoms){
		//Compute binding Sites for Metal ions
			if (_metal_atoms.size()!=0){
				if (_metal_atoms.size()>metal_threshold){
					cout<<"\n Too many metals in "<<accn<<"."<<ext<<"!! System will exhaust!!";
					return false;
				}

		    	cout <<"\nStarting computations for "<<_metal_atoms.at(0).get_atom_symbol()<<" Ions";
		    	for (int i=0; i<_metal_atoms.size(); i++){
		    		metal_binding_site mbs=metal_binding_site(accn, _metal_atoms.at(i), _unique_metal_list, _non_metals, _bond_info);

		    		if (_metal_atoms.at(0).get_atom_symbol()=="NA"){
		    			_na_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="CA"){
		    			_ca_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="K"){
		    			_k_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="MG"){
		    			_mg_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="MN"){
		    			_mn_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="FE"){
		    			_fe_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="CU"){
		    			_cu_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="ZN"){
		    			_zn_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="CO"){
		    			_co_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="MO"){
		    			_mo_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="AG"){
		    			_ag_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="PB"){
		    			_pb_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="HG"){
		    			_hg_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="AS"){
		    			_as_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="CD"){
		    			_cd_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="NI"){
		    			_ni_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="CR"){
		    			_cr_sites.push_back(mbs);
		    		}else if (_metal_atoms.at(0).get_atom_symbol()=="CS"){
		    			_cs_sites.push_back(mbs);
		    		}
		    	}

		    	cout <<"\nComputation Ends for "<<_metal_atoms.at(0).get_atom_symbol()<<" Ions";
			}
		}

		return true;
	}
}
