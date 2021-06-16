#include "ligand.h"

direct_ligand::direct_ligand(atom atm, atom mt, vector<string> metal_list, vector<bond> bond_info){
	_parent_metal=mt;
	_metal_list=metal_list;

	_atom_direct_contact=atm;
	_distance_from_parent_metal=_parent_metal.dist(&_atom_direct_contact);
	_bond_info=bond_info;
	//_ligand_atoms=new ligand;
	//_ligand_atoms->_atm=_atom_direct_contact;
	//_no_of_ligand_atoms=1;
	//find_lig_nbh_atom(_ligand_atoms);
	//find_ligand_tree(_ligand_atoms);
}

atom & direct_ligand::get_atom_direct_contact(){
	return _atom_direct_contact;
}

double direct_ligand::get_distance_direct_contact(){
	return _distance_from_parent_metal;
}

//ligand *direct_ligand::get_ligand_atoms(){
	//return _ligand_atoms;
//}

atom direct_ligand::get_parent_metal(){
	return _parent_metal;
}

vector<string> direct_ligand::get_metal_list(){
	return _metal_list;
}

//int direct_ligand::get_no_of_ligand_atoms(){
	//return _no_of_ligand_atoms;
//}

vector<bond> direct_ligand::get_bond_info(){
	return _bond_info;
}

void direct_ligand::print_direct_ligand(){
	cout<<"\n\nAtom in Direct Contact with Metal "<<_parent_metal.get_atom_name()<<"_"<<_parent_metal.get_serial_no()<<": "<<_atom_direct_contact.get_serial_no()<<"_"<<_atom_direct_contact.get_atom_name()<<"_"<<_atom_direct_contact.get_chain_name()<<"_"<<_atom_direct_contact.get_residue_name()<<"_"<<_atom_direct_contact.get_residue_no();
    cout<<"\nDistance from "<<_parent_metal.get_atom_symbol()<<"_"<<_parent_metal.get_serial_no()<<": "<<_distance_from_parent_metal;
}

/*void direct_ligand::print_direct_ligand(){
	cout<<"\n\nAtom in Direct Contact with Metal "<<_parent_metal.get_atom_name()<<"_"<<_parent_metal.get_serial_no()<<": "<<_atom_direct_contact.get_serial_no()<<"_"<<_atom_direct_contact.get_atom_name()<<"_"<<_atom_direct_contact.get_chain_name()<<"_"<<_atom_direct_contact.get_residue_name()<<"_"<<_atom_direct_contact.get_residue_no();
    cout<<"\nDistance from "<<_parent_metal.get_atom_symbol()<<"_"<<_parent_metal.get_serial_no()<<": "<<_distance_from_parent_metal;
	cout<<"\nTotal No. of Ligand Atoms: "<<_no_of_ligand_atoms;
	cout<<"\nLigand Atoms:";

	queue<ligand *> node_list;
	node_list.push(_ligand_atoms);
	ligand *tnode=node_list.front();

	while (!node_list.empty()){
		atom atm=tnode->_atm;
		cout<<"\nNode Atom: "<<atm.get_serial_no()<<"_"<<atm.get_atom_name()<<"_"<<atm.get_chain_name()<<"_"<<atm.get_residue_name()<<"_"<<atm.get_residue_no();

		cout<<"\n\tNo. of Parent Atom: "<<tnode->_parent.size();
		if (tnode->_parent.size()>0){
			cout<<"\n\tParent Atom(s): ";
			for (auto pn: tnode->_parent){
				cout<<pn->_atm.get_serial_no()<<"_"<<pn->_atm.get_atom_name()<<"_"<<pn->_atm.get_chain_name()<<"_"<<pn->_atm.get_residue_name()<<"_"<<pn->_atm.get_residue_no()<<"--";
			}
		}

		cout<<"\n\tNo. of Children Atom: "<<tnode->_children.size();
		if (tnode->_children.size()>0){
			cout<<"\n\tChildren Atom(s): ";
			for (auto cn: tnode->_children){
				cout<<cn->_atm.get_serial_no()<<"_"<<cn->_atm.get_atom_name()<<"_"<<cn->_atm.get_chain_name()<<"_"<<cn->_atm.get_residue_name()<<"_"<<cn->_atm.get_residue_no()<<"--";
				node_list.push(cn);
			}
		}

		node_list.pop();
		if (!node_list.empty()){
			tnode=node_list.front();
		}
	}
}

void direct_ligand:: find_lig_nbh_atom(ligand* ligand_node){
	atom atm=ligand_node->_atm;

	for (auto bnd: _bond_info){
		if (atm.get_serial_no()==bnd.get_bonded_atom1().get_serial_no()){
			if ((atm.get_chain_name()==bnd.get_bonded_atom2().get_chain_name()) && (atm.get_residue_name()==bnd.get_bonded_atom2().get_residue_name()) && (atm.get_residue_no()==bnd.get_bonded_atom2().get_residue_no())){
				atom* at=new atom;
				copy_atom(at, bnd.get_bonded_atom2());
				ligand_node->_nbh_list.push_back(at);
			}else{
				int k;
				for (k=0; k<_metal_list.size(); k++){
					if (bnd.get_bonded_atom2().get_atom_name()==_metal_list.at(k)){
						break;
					}
				}

				if (k<_metal_list.size()){
					atom* at=new atom;
					copy_atom(at, bnd.get_bonded_atom2());
					ligand_node->_nbh_list.push_back(at);
				}
			}
		}else if (atm.get_serial_no()==bnd.get_bonded_atom2().get_serial_no()){
			if ((atm.get_chain_name()==bnd.get_bonded_atom1().get_chain_name()) && (atm.get_residue_name()==bnd.get_bonded_atom1().get_residue_name()) && (atm.get_residue_no()==bnd.get_bonded_atom1().get_residue_no())){
				atom* at=new atom;
				copy_atom(at, bnd.get_bonded_atom1());
				ligand_node->_nbh_list.push_back(at);
			}else{
				int k;
				for (k=0; k<_metal_list.size(); k++){
					if (bnd.get_bonded_atom1().get_atom_name()==_metal_list.at(k)){
						break;
					}
				}

				if (k<_metal_list.size()){
					atom* at=new atom;
					copy_atom(at, bnd.get_bonded_atom1());
					ligand_node->_nbh_list.push_back(at);
				}
			}
		}
	}
}

void direct_ligand::find_ligand_tree(ligand* node){
	queue<ligand *> node_list;
	node_list.push(node);
	vector<ligand *> visited_node;
	ligand* temp_node=node_list.front();
	ligand* root=NULL;

	while (!node_list.empty()){
		int k;
		do{
			for (k=0; k<visited_node.size(); k++){
				if (temp_node->_atm.get_serial_no()==visited_node.at(k)->_atm.get_serial_no()){
					break;
				}
			}

			node_list.pop();
			if (k==visited_node.size()){
				break;
			}else{
				if (!node_list.empty()){
					temp_node=node_list.front();
				}else{
					temp_node=NULL;
				}
			}
		}while (!node_list.empty());

		if (temp_node!=NULL){
			atom atm=temp_node->_atm;

			for (auto bnd: _bond_info){
				if (atm.get_serial_no()==bnd.get_bonded_atom1().get_serial_no()){
					if (bnd.get_bonded_atom2().get_serial_no()==_parent_metal.get_serial_no()){
						root=new ligand;
						root->_atm=_parent_metal;
						root->_parent.push_back(NULL);
						root->_children.push_back(temp_node);
						temp_node->_parent.push_back(root);
					}else{
						int k;
						ligand *found_visited;
						for (k=0; k<visited_node.size(); k++){
							if (bnd.get_bonded_atom2().get_serial_no()==visited_node.at(k)->_atm.get_serial_no()){
								found_visited=visited_node.at(k);
								break;
							}
						}

						if (k<visited_node.size()){
							if ((root!=NULL) && (found_visited!=root->_children.at(0))){
								temp_node->_parent.push_back(found_visited);
							}
						}else if (k==visited_node.size()){
							vector<ligand *> temp_lig;
							while (node_list.size()>0){
								temp_lig.push_back(node_list.front());
								node_list.pop();
							}

							int p;
							ligand *found_node_list;
							for (p=0; p<temp_lig.size(); p++){
								if (bnd.get_bonded_atom2().get_serial_no()==temp_lig.at(p)->_atm.get_serial_no()){
									found_node_list=temp_lig.at(p);
									break;
								}
							}

							if (p==temp_lig.size()){
								if (((atm.get_chain_name()==bnd.get_bonded_atom2().get_chain_name()) && (atm.get_residue_name()==bnd.get_bonded_atom2().get_residue_name()) && (atm.get_residue_no()==bnd.get_bonded_atom2().get_residue_no())) || (bnd.get_bonded_atom2().get_atom_type()=="HETATM")){
									ligand *temp=new ligand;
									temp->_atm= bnd.get_bonded_atom2();
									find_lig_nbh_atom(temp);
									temp_node->_children.push_back(temp);

									if ((temp_node->_parent.size()>0) && (temp_node->_parent.at(0)==root)){
										temp->_parent.push_back(temp_node);
									}

									temp_lig.push_back(temp);
									_no_of_ligand_atoms++;
								}
							}else{
								temp_node->_children.push_back(found_node_list);
								temp_lig.push_back(found_node_list);
							}

							for (p=0; p<temp_lig.size(); p++){
								node_list.push(temp_lig.at(p));
							}
						}
					}
				}else if (atm.get_serial_no()==bnd.get_bonded_atom2().get_serial_no()){
					if (bnd.get_bonded_atom1().get_serial_no()==_parent_metal.get_serial_no()){
						root=new ligand;
						root->_atm=_parent_metal;
						root->_children.push_back(temp_node);
						temp_node->_parent.push_back(root);
					}else{
						int k;
						ligand *found_visited;
						for (k=0; k<visited_node.size(); k++){
							if (bnd.get_bonded_atom1().get_serial_no()==visited_node.at(k)->_atm.get_serial_no()){
								found_visited=visited_node.at(k);
								break;
							}
						}

						if (k<visited_node.size()){
							if ((root!=NULL) && (found_visited!=root->_children.at(0))){
								temp_node->_parent.push_back(found_visited);
							}
						}else if (k==visited_node.size()){
							vector<ligand *> temp_lig;
							while (node_list.size()>0){
								temp_lig.push_back(node_list.front());
								node_list.pop();
							}

							int p;
							ligand *found_node_list;
							for (p=0; p<temp_lig.size(); p++){
								if (bnd.get_bonded_atom1().get_serial_no()==temp_lig.at(p)->_atm.get_serial_no()){
									found_node_list=temp_lig.at(p);
									break;
								}
							}

							if (p==temp_lig.size()){
								if (((atm.get_chain_name()==bnd.get_bonded_atom1().get_chain_name()) && (atm.get_residue_name()==bnd.get_bonded_atom1().get_residue_name()) && (atm.get_residue_no()==bnd.get_bonded_atom1().get_residue_no())) || (bnd.get_bonded_atom1().get_atom_type()=="HETATM")){
									ligand *temp=new ligand;
									temp->_atm= bnd.get_bonded_atom1();
									find_lig_nbh_atom(temp);
									temp_node->_children.push_back(temp);

									if ((temp_node->_parent.size()>0) && (temp_node->_parent.at(0)==root)){
										temp->_parent.push_back(temp_node);
									}

									temp_lig.push_back(temp);
									_no_of_ligand_atoms++;
								}
							}else{
								temp_node->_children.push_back(found_node_list);
								temp_lig.push_back(found_node_list);
							}

							for (p=0; p<temp_lig.size(); p++){
								node_list.push(temp_lig.at(p));
							}
						}
					}
				}
			}

			visited_node.push_back(temp_node);
			if (!node_list.empty()){
				temp_node=node_list.front();
			}else{
				temp_node=NULL;
			}
		}
	}
}*/

void direct_ligand::copy_atom(atom *atm1, atom atm2){
	atm1->set_atom_type(atm2.get_atom_type());
	atm1->set_serial_no(atm2.get_serial_no());
	atm1->set_atom_name(atm2.get_atom_name());
	atm1->set_alt_conf_indicator(atm2.get_alt_conf_indicator());
	atm1->set_residue_name(atm2.get_residue_name());
	atm1->set_chain_name(atm2.get_chain_name());
	atm1->set_residue_no(atm2.get_residue_no());
	atm1->set_occupancy(atm2.get_occupancy());
	atm1->set_temp_factor(atm2.get_temp_factor());
	atm1->set_atom_symbol(atm2.get_atom_symbol());
	atm1->set_partial_charge(atm2.get_partial_charge());
	atm1->set_molecule_type(atm2.get_molecule_type());
	atm1->set_residue_type(atm2.get_residue_type());
	atm1->set_X(atm2.X());
	atm1->set_Y(atm2.Y());
	atm1->set_Z(atm2.Z());
}
