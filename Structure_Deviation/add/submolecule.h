#ifndef SUBMOLECULE_H_
#define SUBMOLECULE_H_

#include <iostream>
#include <iostream>
#include <string>
#include <cmath>

#include "geometry.h"

using namespace geom;

class atom: public Sphere{
    string _type; // ATOM, HETATM etc.
    long int _serial_no;
    string _atom_name;
    string _alt_conf_indicator;
    string _residue_name;
    string _chain_name;
    long int _residue_no;
    double _x_coord;
    double _y_coord;
    double _z_coord;
    double _occupancy;
    double _temp_factor;
    string _atom_symbol;
    double _partial_charge;
	string _molecule_type;
	string _residue_type;

	void assign_molecule_and_residue_type();

public:
	atom(){}
	atom(string type, long int serial, string atom_name, string alt_ind, string res_name, string chain_name, long int res_no, double x, double y, double z, double occu, double t_fact, double radius, string atom_symbol, double partial_charge);
	void fprint(FILE*);
	void set_partial_charge(double);
	double get_partial_charge();
	void set_atom_type(string);
	string get_atom_type();
	void set_serial_no(long int);
	long int get_serial_no();
	void set_atom_name(string);
	string get_atom_name();
	void set_alt_conf_indicator(string);
	string get_alt_conf_indicator();
	void set_residue_name(string);
	string get_residue_name();
	void set_chain_name(string);
	string get_chain_name();
	void set_residue_no(long int);
	long int get_residue_no();
	void set_occupancy(double);
	double get_occupancy();
	void set_temp_factor(double);
	double get_temp_factor();
	void set_atom_symbol(string);
	string get_atom_symbol();
	void set_molecule_type(string);
	string get_molecule_type();
	void set_residue_type(string);
	string get_residue_type();
	void set_X(double);
	void set_Y(double);
	void set_Z(double);
};

class residue{
	vector<atom> _atom_group;

public:
	residue(vector<atom>);
	int get_no_of_atoms();
	vector<atom> get_atoms();
};

class bond{
	atom _atom1;
	atom _atom2;
	string _bond_type;
	bool _visited_atom1;	//Visiting Status of Atom1 while detecting Nbh List
	bool _visited_atom2;	//Visiting Status of Atom2 while detecting Nbh List
	double _bond_length;
	double _bond_strength;

	void compute_bond_strength();

public:
	bond(){}
	bond(atom, atom, string, bool, bool);
	atom get_bonded_atom1();
	atom get_bonded_atom2();
	string get_bond_type();
	void set_visiting_status_atom1(bool);
	bool get_visiting_status_atom1();
	void set_visiting_status_atom2(bool);
	bool get_visiting_status_atom2();
	double get_bond_length();
	double get_bond_strength();
};

#endif /* SUBMOLECULE_H_ */
