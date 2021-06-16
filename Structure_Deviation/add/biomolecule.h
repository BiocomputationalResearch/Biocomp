#ifndef BIOMOLECULE_H_
#define BIOMOLECULE_H_

#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <gemmi/cif.hpp>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>

#include "submolecule.h"
#include "coordinationsite.h"

namespace cif = gemmi::cif;

class protein{
	vector<residue> _residue_group;
	vector<atom> _na_atoms;
	vector<atom> _ca_atoms;
	vector<atom> _k_atoms;
	vector<atom> _mg_atoms;
	vector<atom> _mn_atoms;
	vector<atom> _fe_atoms;
	vector<atom> _cu_atoms;
	vector<atom> _zn_atoms;
	vector<atom> _co_atoms;
	vector<atom> _mo_atoms;
	vector<atom> _ag_atoms;
	vector<atom> _pb_atoms;
	vector<atom> _hg_atoms;
	vector<atom> _as_atoms;
	vector<atom> _cd_atoms;
	vector<atom> _ni_atoms;
	vector<atom> _cr_atoms;
	vector<atom> _cs_atoms;
	vector<string> _metal_list;
	vector<string> _unique_metal_list;
	vector<int> _count_unique_metal;
	vector<atom> _non_metals;
	vector<metal_binding_site> _na_sites;
	vector<metal_binding_site> _ca_sites;
	vector<metal_binding_site> _k_sites;
	vector<metal_binding_site> _mg_sites;
	vector<metal_binding_site> _mn_sites;
	vector<metal_binding_site> _fe_sites;
	vector<metal_binding_site> _cu_sites;
	vector<metal_binding_site> _zn_sites;
	vector<metal_binding_site> _co_sites;
	vector<metal_binding_site> _mo_sites;
	vector<metal_binding_site> _ag_sites;
	vector<metal_binding_site> _pb_sites;
	vector<metal_binding_site> _hg_sites;
	vector<metal_binding_site> _as_sites;
	vector<metal_binding_site> _cd_sites;
	vector<metal_binding_site> _ni_sites;
	vector<metal_binding_site> _cr_sites;
	vector<metal_binding_site> _cs_sites;
	vector<bond> _bond_info;

	enum feature_index_atom{min_atom=0, atom_type=min_atom, serial_no, atom_name, alt_conf_indicator, residue_name, chain_id, residue_no, x_coord, y_coord, z_coord, occupancy, t_factor, atom_symbol, max_atom=atom_symbol};

	vector<vector<string> > read_mol_atom_info(string);
	vector<vector<string> > read_mol_bond_info(string);
	bool read_cif_atoms(string, string);
	void assign_bond_info(vector<atom>, vector<vector<string> >);
	void make_residue_metal_non_metal_group(vector<atom>);
	int set_pdb_atom_field_index(vector<vector<int> > &);
	bool read_pdb_atoms(string, string);
	void split_pdb_atom_str(string, vector<string> &, vector<vector<int> > &);
	void separate_metal_and_non_metals(vector<atom>);
	bool compute_binding_sites(string, string, int);
	void compute_inter_atomic_distances(FILE *, int, string, long int, string, vector<direct_ligand>);
	void compute_atom_to_metal_angles(FILE *, int, string, long int, string, vector<direct_ligand>);
	void generate_general_report(string);

public:
	protein(string, string, int);
	vector<atom> get_na_atoms();
	vector<atom> get_ca_atoms();
	vector<atom> get_k_atoms();
	vector<atom> get_mg_atoms();
	vector<atom> get_mn_atoms();
	vector<atom> get_fe_atoms();
	vector<atom> get_cu_atoms();
	vector<atom> get_zn_atoms();
	vector<atom> get_co_atoms();
	vector<atom> get_mo_atoms();
	vector<atom> get_ag_atoms();
	vector<atom> get_pb_atoms();
	vector<atom> get_hg_atoms();
	vector<atom> get_as_atoms();
	vector<atom> get_cd_atoms();
	vector<atom> get_ni_atoms();
	vector<atom> get_cr_atoms();
	vector<atom> get_cs_atoms();
	vector<vector<atom> > get_all_metal_atoms();
	vector<string> get_metal_list();
	vector<string> get_unique_metal_list();
	vector<int> get_count_unique_metal();
	vector<atom> get_non_metal_atoms();
	vector<metal_binding_site> get_na_binding_sites();
	vector<metal_binding_site> get_ca_binding_sites();
	vector<metal_binding_site> get_k_binding_sites();
	vector<metal_binding_site> get_mg_binding_sites();
	vector<metal_binding_site> get_mn_binding_sites();
	vector<metal_binding_site> get_fe_binding_sites();
	vector<metal_binding_site> get_cu_binding_sites();
	vector<metal_binding_site> get_zn_binding_sites();
	vector<metal_binding_site> get_co_binding_sites();
	vector<metal_binding_site> get_mo_binding_sites();
	vector<metal_binding_site> get_ag_binding_sites();
	vector<metal_binding_site> get_pb_binding_sites();
	vector<metal_binding_site> get_hg_binding_sites();
	vector<metal_binding_site> get_as_binding_sites();
	vector<metal_binding_site> get_cd_binding_sites();
	vector<metal_binding_site> get_ni_binding_sites();
	vector<metal_binding_site> get_cr_binding_sites();
	vector<metal_binding_site> get_cs_binding_sites();
	vector<vector<metal_binding_site> > get_all_metal_sites();
	vector<bond> get_bond_info();
};

#endif /* BIOMOLECULE_H_ */
