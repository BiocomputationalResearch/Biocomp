//====================================================================================================================================================================================================================================================================
// Name        : Structure_Deviation.cpp
// Author      : Parthajit Roy and Swati Adhikari
// Version     : 1.1
// Copyright   : Any software that uses our code must abide by terms of the GNU Public License version 2 and must refer our papers 'A Geometry Based Algorithm for Comparison of Tetrahedral Metal Binding Sites' and 'Algorithms to compare some metal binding sites'
// Description : Provide threshold value for number of each metals in an input file that can be processed with -M option and provide input file in .pdb/.cif format while executing the program
//====================================================================================================================================================================================================================================================================

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <bits/stdc++.h>
#include <sstream>
#include "./biomolecule.h"
#include "coordinationsite.h"
#include "globalMain.h"
#include "mainIncludeFunc.h"
#define GetCurrentDir getcwd
#define default_metal_threshold 5

using namespace std;

int main(int argc,char **argv)
{
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	string current_working_dir(buff);
	string _file_path = current_working_dir;
	cout << "Current path is " << _file_path << '\n';

	string protein_name="";
	vector<vector<metal_binding_site> > all_metal_sites;
	vector<string> _metal_name;
	int metal_threshold, ip;

	if (argc==1){
		cout<<"Provide Input Protein's .pdb/.cif File and -M option for max. no. of each metal to be processed";
		exit(1);
	}else if (argc==2){
		int len_option=strlen(argv[1]);

		if (len_option==2){
			char option[2]={argv[1][0], argv[1][1]};

			if (option[0]=='-'){
				if (option[1]=='M'){
					cout<<"\nProvide value for -M option and the Input .pdb/.cif File";
					exit(1);
				}else{
					cout<<"\nWrong Option!!";
					exit(1);
				}
			}else{
				cout<<"\nWrong Option!!";
				exit(1);
			}
		}else if (len_option>2){
			metal_threshold=default_metal_threshold;
			ip=0;
		}else{
			cout<<"\nWrong Option!!";
			exit(1);
		}
	}else if (argc==3){
		int len_option=strlen(argv[1]);

		if (len_option==2){
			char option[2]={argv[1][0], argv[1][1]};

			if (option[0]=='-'){
				if (option[1]=='M'){
					bool isNumber=testNumber(argv[2]);

					if (isNumber==true){
						cout<<"\nProvide Input Protein's .pdb/.cif File";
						exit(1);
					}else{
						cout<<"\nProvide integer value for -M option";
						exit(1);
					}
				}else{
					cout<<"\nWrong Option!!";
					exit(1);
				}
			}else{
					cout<<"\nWrong Option!!";
				exit(1);
			}
		}else if (len_option>2){
			metal_threshold=default_metal_threshold;
			ip=0;
		}else{
			cout<<"\nWrong Option!!";
			exit(1);
		}
	}else if (argc>=4){
		int len_option=strlen(argv[1]);

		if (len_option==2){
			char option[2]={argv[1][0], argv[1][1]};

			if (option[0]=='-'){
				if (option[1]=='M'){
					bool isNumber=testNumber(argv[2]);

					if (isNumber==true){
						metal_threshold=atoi(argv[2]);
						ip=2;
					}else{
						cout<<"\nProvide integer value for -M option";
						exit(1);
					}
				}else{
					cout<<"\nWrong Option!!";
					exit(1);
				}
			}else{
					cout<<"\nWrong Option!!";
				exit(1);
			}
		}else if (len_option>2){
			metal_threshold=default_metal_threshold;
			ip=0;
		}else{
			cout<<"\nWrong Option!!";
			exit(1);
		}
	}

	while(++ip<argc){
		protein_name=argv[ip];
		fflush(stdin);
		int result=access(protein_name.c_str(), F_OK);

		if (result==0){
			int p1 = protein_name.find_last_of(".");

			if(p1 == -1){
				cout<<"\nInvalid File name: "<<protein_name<<endl;
			}else{
				string _accn = protein_name.substr(0, p1);
				string _ext = protein_name.substr(p1);
				string _save_ext=_ext;

				//Convert to MOL File
			    transform(_ext.begin(), _ext.end(), _ext.begin(), ::toupper);
			    string mol_file=convert_to_mol2(protein_name, _accn, _ext);

			    if (mol_file!="NOT_CONVERTED"){
			    	if ((_save_ext!=".cif") && (_save_ext!=".pdb")){
			    		cout<<"\nInvalid file extension... please provide cif or pdb files: "<<protein_name<<endl;
			    	}else{
			    		cout<<"\n\nProtein To Be Considered:";
			    		protein prot=protein(protein_name, mol_file, metal_threshold);

			    		if (prot.get_all_metal_sites().size()!=0){
			    			all_metal_sites=prot.get_all_metal_sites();

			    			for (auto met_sites: all_metal_sites){
			    				_metal_name.push_back(met_sites.at(0).get_metal_atom().get_atom_symbol());
			    				vector<metal_binding_site> tripl;
			    				vector<metal_binding_site> tetrapl;
			    				vector<metal_binding_site> tetra;
			    				vector<metal_binding_site> tetrapyr;
			    				vector<metal_binding_site> tribipyr;
			    				vector<metal_binding_site> octa;

			    				for (int i=0; i<met_sites.size(); i++){
			    					string shape=met_sites.at(i).get_shape();

			    					if ((shape=="trigonal_planar") || (shape=="trigonal_pyramidal") || (shape=="t_shaped")){
			    						tripl.push_back(met_sites.at(i));
			    					}else if((shape=="square_planar_metal_on_plane") || (shape=="square_planar_metal_outside") || (shape=="tetragonal_planar_metal_on_plane") || (shape=="tetragonal_planar_metal_outside")){
			    						tetrapl.push_back(met_sites.at(i));
			    					}else if((shape=="tetrahedral") || (shape=="regular_tetrahedral")){
			    						tetra.push_back(met_sites.at(i));
			    					}else if(shape=="tetragonal_pyramidal"){
			    						tetrapyr.push_back(met_sites.at(i));
			    					}else if(shape=="trigonal_bipyramidal"){
			    						tribipyr.push_back(met_sites.at(i));
			    					}else if((shape=="octahedral") || (shape=="tetragonal_bipyramidal")){
			    						octa.push_back(met_sites.at(i));
			    					}
			    				}

			    				if (tripl.size()!=0){
			    					form_tripl_sites(tripl, _metal_name.at(_metal_name.size()-1), _accn);
			    				}

			    				if (tetrapl.size()!=0){
			    					form_tetrapl_sites(tetrapl, _metal_name.at(_metal_name.size()-1), _accn);
			    				}

			    				if (tetra.size()!=0){
			    					form_tetra_sites(tetra, _metal_name.at(_metal_name.size()-1), _accn);
			    				}

			    				if (tetrapyr.size()!=0){
			    					form_tetrapyr_sites(tetrapyr, _metal_name.at(_metal_name.size()-1), _accn);
			    				}

			    				if (tribipyr.size()!=0){
			    					form_tribipyr_sites(tribipyr, _metal_name.at(_metal_name.size()-1), _accn);
			    				}

			    				if (octa.size()!=0){
			    					form_octa_sites(octa, _metal_name.at(_metal_name.size()-1), _accn);
			    				}
			    			}
			    		}
			    	}
			    }else{
			    	cout<<"\nThis Protein "<<protein_name<<" can't be converted to MOL File!!\n";
			    	exit(1);
			    }
			}
		}else{
			cout<<"\nFile Does Not Exist!";
			exit(1);
		}

		//cout<<"\n\nPress Enter to Continue:";
		//int ch=cin.get();
	}

	if ((no_of_compared_tripl_na!=0) || (no_of_compared_tripl_ca!=0) || (no_of_compared_tripl_k!=0) || (no_of_compared_tripl_mg!=0) || (no_of_compared_tripl_mn!=0) || (no_of_compared_tripl_fe!=0) || (no_of_compared_tripl_cu!=0) || (no_of_compared_tripl_zn!=0) || (no_of_compared_tripl_co!=0) || (no_of_compared_tripl_pb!=0) || (no_of_compared_tripl_hg!=0) || (no_of_compared_tripl_as!=0) || (no_of_compared_tripl_cd!=0) || (no_of_compared_tripl_ni!=0)){
		generate_summary_report_tripl();
	}else {
		cout<<"No Trigonal Planar Site for Test Metals Exists";
	}

	if ((no_of_compared_tetrapl_mg!=0) || (no_of_compared_tetrapl_fe!=0) || (no_of_compared_tetrapl_cd!=0) || (no_of_compared_tetrapl_ni!=0)){
		generate_summary_report_tetrapl();
	}else {
		cout<<"No Tetragonal Planar Site for Test Metals Exists";
	}

	if ((no_of_compared_tetra_na!=0) || (no_of_compared_tetra_ca!=0) || (no_of_compared_tetra_k!=0) || (no_of_compared_tetra_mg!=0) || (no_of_compared_tetra_mn!=0) || (no_of_compared_tetra_fe!=0) || (no_of_compared_tetra_cu!=0) || (no_of_compared_tetra_zn!=0) || (no_of_compared_tetra_co!=0) || (no_of_compared_tetra_mo!=0) || (no_of_compared_tetra_ag!=0) || (no_of_compared_tetra_pb!=0) || (no_of_compared_tetra_hg!=0) || (no_of_compared_tetra_as!=0) || (no_of_compared_tetra_cd!=0) || (no_of_compared_tetra_ni!=0) || (no_of_compared_tetra_cr!=0) || (no_of_compared_tetra_cs!=0)){
		generate_summary_report_tetra();
	}else {
		cout<<"No Tetrahedral Site for Test Metals Exists";
	}

	if ((no_of_compared_tetrapyr_na!=0) || (no_of_compared_tetrapyr_ca!=0) || (no_of_compared_tetrapyr_mg!=0) || (no_of_compared_tetrapyr_cd!=0) || (no_of_compared_tetrapyr_ni!=0) || (no_of_compared_tetrapyr_zn!=0)){
		generate_summary_report_tetrapyr();
	}else {
		cout<<"No Tetragonal Pyramidal Site for Test Metals Exists";
	}

	if ((no_of_compared_tribipyr_na!=0) || (no_of_compared_tribipyr_ca!=0) || (no_of_compared_tribipyr_mg!=0) || (no_of_compared_tribipyr_mn!=0) || (no_of_compared_tribipyr_fe!=0) || (no_of_compared_tribipyr_zn!=0) || (no_of_compared_tribipyr_co!=0) || (no_of_compared_tribipyr_hg!=0) || (no_of_compared_tribipyr_cd!=0) || (no_of_compared_tribipyr_ni!=0)){
		generate_summary_report_tribipyr();
	}else {
		cout<<"No Trigonal Bipyramidal Site for Test Metals Exists";
	}

	if ((no_of_compared_octa_na!=0)  || (no_of_compared_octa_ca!=0) || (no_of_compared_octa_mg!=0) || (no_of_compared_octa_mn!=0) || (no_of_compared_octa_fe!=0) || (no_of_compared_octa_cd!=0) || (no_of_compared_octa_ni!=0) || (no_of_compared_octa_zn!=0) || (no_of_compared_octa_co!=0)){
		generate_summary_report_octa();
	}else {
		cout<<"No Octahedral Site for Test Metals Exists";
	}

	return 0;
}
