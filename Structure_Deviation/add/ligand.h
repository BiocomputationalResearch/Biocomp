#ifndef LIGAND_H_
#define LIGAND_H_

#include <iostream>
#include <string>
#include<stack>
#include<queue>

#include "submolecule.h"

/*struct ligand_denticity{
	int _denticity;
	long int _sr_no;
	string _chain_id;
	string _res_name;
	long int _res_no;
	vector<long int> _dented_atoms_sr_no;
	vector<string> _dented_atoms_name;
	int _ligated_atoms_ring_size;
};

struct ligand{
	atom _atm;
	vector<atom *> _nbh_list;
	vector<ligand *> _parent;
	vector<ligand *> _children;
};*/

class direct_ligand{
	atom _parent_metal;
	vector<string> _metal_list;
	atom _atom_direct_contact;
	double _distance_from_parent_metal;
	//ligand* _ligand_atoms;
	//int _no_of_ligand_atoms;
	vector<bond> _bond_info;
	//vector<indirect_ligand> _indirect_res_gr;

	/*void print(ligand *lig_atm){
		atom atm=lig_atm->_atm;
		cout<<atm.get_serial_no()<<"_"<<atm.get_atom_name()<<"_"<<atm.get_chain_name()<<"_"<<atm.get_residue_name()<<"_"<<atm.get_residue_no()<<"--";

		if (lig_atm->_children.size()!=0){
			for (auto nd: lig_atm->_children){
				print(nd);
			}
		}else{
			return;
		}
	}*/

	void copy_atom(atom *, atom);
public:
	direct_ligand(atom, atom, vector<string>, vector<bond>);
	//direct_ligand(atom, atom, vector<string>, vector<bond>, double);
	//void find_ligand_tree(ligand *);
	//void find_lig_nbh_atom(ligand *);
	atom & get_atom_direct_contact();
	double get_distance_direct_contact();
	//ligand *get_ligand_atoms();
	atom get_parent_metal();
	vector<string> get_metal_list();
	//int get_no_of_ligand_atoms();
	vector<bond> get_bond_info();
	void print_direct_ligand();
};

#endif /* LIGAND_H_ */
