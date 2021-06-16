#include "coordinationsite.h"
#include <set>

metal_binding_site::metal_binding_site(string accn, atom atm, vector<string> metals, vector<atom> non_metals, vector<bond> bond_info){
	_accn=accn;
	_metal_atom=atm;
	_metal_list=metals;
	_non_metals=non_metals;
	_bond_info=bond_info;
	_inner_coord_no=0;
	_shape_of_inner_sphere="zero_coordinated";
	_mean_distance_direct_ligands_from_metal=0.0;
	_mean_inter_atomic_distance=0.0;
	_mean_angle_atom_metal_atom=0.0;
	_site_type="zero_coordinated";
	det_threshold=3.0;
	dist_threshold=0.5;

	compute_nearest_atom();
	//print();

	//if (_inner_coord_no!=0){
		//delete_ligand_tree();
	//}
}

vector<atom> metal_binding_site::get_nearest_atoms(){
	return _nearest_atoms;
}

string metal_binding_site::get_site_type(){
	return _site_type;
}

vector<direct_ligand> metal_binding_site::get_inner_coord_sphere(){
	return _inner_coord_sphere;
}

vector<long int> metal_binding_site::get_planar_atoms(){
	return _planar_atoms;
}

vector<long int> metal_binding_site::get_non_planar_atoms(){
	return _non_planar_atoms;
}

vector<atom> metal_binding_site::get_water_molecules(){
	return _water_group;
}

int metal_binding_site::get_inner_coord_no(){
	return _inner_coord_no;
}

string metal_binding_site::get_shape(){
	return _shape_of_inner_sphere;
}

atom & metal_binding_site::get_metal_atom(){
	return _metal_atom;
}

void metal_binding_site::print(){
	if (_nearest_atoms.size()!=0){
		cout<<"\n\nBinding Site for "<<_metal_atom.get_atom_name()<<" Ion - Sr. No. "<<_metal_atom.get_serial_no()<<endl;
		for (int i=0; i<_nearest_atoms.size(); i++){
			_nearest_atoms.at(i).fprint(stdout);
		}
	}

	if (_inner_coord_sphere.size()!=0){
		cout<<"\n\nCoordination No.: "<<_inner_coord_no<<endl;
    	cout<<"\n\nAtom in First Coordination Sphere:";
		cout<<"\nDirect Ligand Atoms:";
		for (int i=0; i<_inner_coord_sphere.size(); i++){
	    	_inner_coord_sphere.at(i).print_direct_ligand();
	    }
	}
}

void metal_binding_site::compute_nearest_atom(){
	for (auto bnd: _bond_info){
		if (_metal_atom.get_serial_no()==bnd.get_bonded_atom1().get_serial_no()){
	    	_nearest_atoms.push_back(bnd.get_bonded_atom2());
		}else if (_metal_atom.get_serial_no()==bnd.get_bonded_atom2().get_serial_no()){
	    	_nearest_atoms.push_back(bnd.get_bonded_atom1());
		}
	}

	if (_nearest_atoms.size()!=0){
		cout<<"\nPre Size: "<<_nearest_atoms.size();
		make_order_direct_atoms();
		cout<<"\nPost Size: "<<_nearest_atoms.size();

		separate_atoms();
		assign_site_type();
		_inner_coord_no=_inner_coord_sphere.size();
		cout<<"\n\nSIZE: "<<_inner_coord_no;

		if (_inner_coord_no!=0){
			determine_shape();
			compute_atom_to_metal_distances();


			if (_inner_coord_no>1){
				compute_inter_atomic_distances();

			}
		}
	}else{
		cout<<"\n\n      No binding atom found for this "<<_metal_atom.get_atom_name()<<",  Sr. No. "<<_metal_atom.get_serial_no()<<endl;
	}
}

void metal_binding_site::make_order_direct_atoms(){
	//cout<<"\n\nStart Order";
	/*cout<<"\n\nBefore Sorted X: ";
	for (auto atm: _nearest_atoms){
		cout<<atm.get_serial_no()<<", ";
	}*/

	for (int x=0; x<(_nearest_atoms.size()-1); x++){
		for (int y=(x+1); y<_nearest_atoms.size(); y++){
			if (_nearest_atoms.at(x).X()>_nearest_atoms.at(y).X()){
				atom atmx=_nearest_atoms.at(x);
				atom atmy=_nearest_atoms.at(y);
				_nearest_atoms.at(x)=atmy;
				_nearest_atoms.at(y)=atmx;
			}
		}
	}
/*cout<<"\n\nOK Sorted X: ";
for (auto atm: _nearest_atoms){
	cout<<atm.get_serial_no()<<", ";
}*/

	int x=0, y=1;
	vector<atom> sorted_atoms_y, sort_reverse_y;
	while(x<(_nearest_atoms.size()-1) && y<_nearest_atoms.size()){
		//cout<<"\n\nX: "<<x<<", Y: "<<y;
		if (_nearest_atoms.at(x).Y()>=_nearest_atoms.at(y).Y()){
			atom atmx=_nearest_atoms.at(x);
			sorted_atoms_y.push_back(atmx);
			x=y;
			y++;
		}else{
			atom atmy=_nearest_atoms.at(y);
			sort_reverse_y.push_back(atmy);
			y++;
		}
	}

	/*cout<<"\n\nOK Sorted Y:";
	for (auto atm: sorted_atoms_y){
		cout<<atm.get_serial_no()<<", ";
	}

	cout<<"\n\nOK Sorted Reverse Y:";
	for (auto atm: sort_reverse_y){
		cout<<atm.get_serial_no()<<", ";
	}*/

	atom atmx=_nearest_atoms.at(x);
	sorted_atoms_y.push_back(atmx);

	/*cout<<"\n\nOK Sorted1 Y:";
	for (auto atm: sorted_atoms_y){
		cout<<atm.get_serial_no()<<", ";
	}*/

	if (sort_reverse_y.size()>0){
		for (int z=(sort_reverse_y.size()-1); z>=0; z--){
			sorted_atoms_y.push_back(sort_reverse_y.at(z));
		}
	}
//cout<<"\nSize Sorted Y: "<<sorted_atoms_y.size();

	x=0;
	y=1;
	vector<atom> sorted_atoms_z, sort_rest_z;
	while(x<(sorted_atoms_y.size()-1) && y<sorted_atoms_y.size()){
		//cout<<"\n\nX: "<<x<<", Y: "<<y;
		if (sorted_atoms_y.at(x).Z()>=sorted_atoms_y.at(y).Z()){
			atom atm=sorted_atoms_y.at(x);
			sorted_atoms_z.push_back(atm);
			x=y;
		}else{
			atom atmy=sorted_atoms_y.at(y);
			sort_rest_z.push_back(atmy);
		}

		y++;
	}

	atom atmxx=sorted_atoms_y.at(x);
	sorted_atoms_z.push_back(atmxx);

	/*cout<<"\n\nOK Sorted2 Z:";
	for (auto atm: sorted_atoms_z){
		cout<<atm.get_serial_no()<<", ";
	}*/

	if (sort_rest_z.size()>0){
		for (int z=0; z<sort_rest_z.size(); z++){
			sorted_atoms_z.push_back(sort_rest_z.at(z));
		}
	}

	_nearest_atoms.clear();
	cout<<"\n\nOrdered Atoms: ";
	if (sorted_atoms_z.size()>0){
		//cout<<"\nSize Sorted Z: "<<sorted_atoms_z.size();
		for (auto atm: sorted_atoms_z){
			_nearest_atoms.push_back(atm);
			cout<<atm.get_serial_no()<<", ";
		}
	}else{
		//cout<<"\nSize Sorted Z: "<<sorted_atoms_y.size();
		for (auto atm: sorted_atoms_y){
			_nearest_atoms.push_back(atm);
			cout<<atm.get_serial_no()<<", ";
		}
	}

	//cout<<"\n\nEnd Order";
}

void metal_binding_site::separate_atoms(){
	char *end;

	for (auto atm: _nearest_atoms){
		if ((atm.get_alt_conf_indicator()==" ") || (atm.get_alt_conf_indicator()=="A" || strtol(atm.get_alt_conf_indicator().c_str(), &end, 10)==1)){
			direct_ligand dl=direct_ligand(atm, _metal_atom, _metal_list, _bond_info);
			_inner_coord_sphere.push_back(dl);
		}
	}
}

void metal_binding_site::assign_site_type(){
	int aa=0, dna=0, rna=0, other=0;
	for (int x=0; x<_inner_coord_sphere.size(); x++){
		if (_inner_coord_sphere.at(x).get_atom_direct_contact().get_molecule_type()=="amino_acid"){
			aa++;
		}else if (_inner_coord_sphere.at(x).get_atom_direct_contact().get_molecule_type()=="dna_molecule"){
			dna++;
		}else if (_inner_coord_sphere.at(x).get_atom_direct_contact().get_molecule_type()=="rna_molecule"){
			rna++;
		}else if (_inner_coord_sphere.at(x).get_atom_direct_contact().get_molecule_type()=="other_biom"){
			other++;
		}
	}

	if (_water_group.size()==_inner_coord_sphere.size()){
		_site_type="all_water";
	}else if ((aa!=0) && (dna==0) && (rna==0)){
		_site_type="protein";
	}else if ((aa==0) && (dna!=0) && (rna==0)){
		_site_type="dna";
	}else if ((aa==0) && (dna==0) && (rna!=0)){
		_site_type="rna";
	}else if ((aa!=0) && (dna!=0) && (rna==0)){
		_site_type="protein_dna_complex";
	}else if ((aa!=0) && (dna==0) && (rna!=0)){
		_site_type="protein_rna_complex";
	}else if ((aa==0) && (dna!=0) && (rna!=0)){
		_site_type="dna_rna_complex";
	}else if ((aa!=0) && (dna!=0) && (rna!=0)){
		_site_type="protein_dna_rna_complex";
	}else{
		_site_type="other";
	}
}

/*void metal_binding_site::delete_ligand_tree(){
	ligand *root=_inner_coord_sphere.at(0).get_ligand_atoms()->_parent.at(0);

	for (auto dlig: _inner_coord_sphere){
		ligand *start_pt=dlig.get_ligand_atoms();
		queue<ligand *> visit_pt;
		visit_pt.push(start_pt);
		vector<ligand *> delete_pt;

		while (!visit_pt.empty()){
			ligand *curr_lig=visit_pt.front();
			bool found=false;
			for (auto pt: delete_pt){
				if (pt==curr_lig){
					found=true;
				}

			}

			if (found==false){
				delete_pt.push_back(curr_lig);
			}

			visit_pt.pop();

			if (curr_lig->_children.size()!=0){
				for (auto ch: curr_lig->_children){
					visit_pt.push(ch);
				}
			}
		}

		for (auto del_lig: delete_pt){
			ligand *temp=del_lig;

			if (temp->_nbh_list.size()!=0){
				for (auto atm: temp->_nbh_list){
					atom *at=atm;
					delete at;
				}
			}

			delete temp;
		}
	}

	delete root;
}*/
