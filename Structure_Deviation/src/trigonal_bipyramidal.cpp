#include "trigonal_bipyramidal.h"

void tribipyrMetalicStructure::print(atom mt, atom* non_mt, int no_of_atoms){
	mt.fprint(stdout);

	for (int i=0; i<no_of_atoms; i++){
		non_mt[i].fprint(stdout);
	}
}

void tribipyrMetalicStructure::print(Point3D* vert, int no_of_atoms){
	for (int i=0; i<no_of_atoms; i++){
		cout<<"\nVertex "<<i<<": "<<vert[i].X()<<"  "<<vert[i].Y()<<"  "<<vert[i].Z();
	}
}

tribipyrMetalicStructure::tribipyrMetalicStructure(){
	_no_of_binding_atoms=5;
	assert(_no_of_binding_atoms==5 && "For Trigonal Bipyramidal No. of Vertices should be 5!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=9;
	assert(_no_of_sides==9 && "For Trigonal Bipyramidal No. of Sides should be 9!!");

	_no_of_planar_atoms=3;
	assert(_no_of_planar_atoms==3 && "For Trigonal Bipyramidal No. of Planar Atoms should be 3!!");

	_no_of_non_planar_atoms=2;
	assert(_no_of_non_planar_atoms==2 && "For Trigonal Bipyramidal No. of Non Planar Atoms should be 2!!");

	_no_of_vertex_vertex_angles=21;
	assert(_no_of_vertex_vertex_angles==21 && "For Trigonal Bipyramidal No. of Angles of Two Vertices with Another Vertex should be 21!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Trigonal Bipyramidal No. of Features to Compare should be 2!!");
}

tribipyrMetalicStructure::tribipyrMetalicStructure(tetrahedron upper, atom upper_atoms[], tetrahedron lower, atom lower_atoms[], atom metal_atom){
	_no_of_binding_atoms=5;
	assert(_no_of_binding_atoms==5 && "For Trigonal Bipyramidal No. of Vertices should be 5!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=9;
	assert(_no_of_sides==9 && "For Trigonal Bipyramidal No. of Sides should be 9!!");

	_no_of_planar_atoms=3;
	assert(_no_of_planar_atoms==3 && "For Trigonal Bipyramidal No. of Planar Atoms should be 3!!");

	_no_of_non_planar_atoms=2;
	assert(_no_of_non_planar_atoms==2 && "For Trigonal Bipyramidal No. of Non Planar Atoms should be 2!!");

	_no_of_vertex_vertex_angles=21;
	assert(_no_of_vertex_vertex_angles==21 && "For Trigonal Bipyramidal No. of Angles of Two Vertices with Another Vertex should be 21!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Trigonal Bipyramidal No. of Features to Compare should be 2!!");

	_metal_atom=metal_atom;

	_upper_tetra=upper;
	for (int i=0; i<4; i++){
		_upper_tetra_atoms[i]=upper_atoms[i];
	}

	_lower_tetra=lower;
	for (int i=0; i<4; i++){
		_lower_tetra_atoms[i]=lower_atoms[i];
	}
}

	tribipyrMetalicStructure::tribipyrMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom, vector<long int> platom, vector<long int> nonplatom){
		_no_of_binding_atoms=5;
		assert(_no_of_binding_atoms==5 && "For Trigonal Bipyramidal No. of Vertices should be 5!!");

		_no_of_coords=3;
		assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

		_no_of_sides=9;
		assert(_no_of_sides==9 && "For Trigonal Bipyramidal No. of Sides should be 9!!");

		_no_of_planar_atoms=3;
		assert(_no_of_planar_atoms==3 && "For Trigonal Bipyramidal No. of Planar Atoms should be 3!!");

		_no_of_non_planar_atoms=2;
		assert(_no_of_non_planar_atoms==2 && "For Trigonal Bipyramidal No. of Non Planar Atoms should be 2!!");


		_no_of_vertex_vertex_angles=21;
		assert(_no_of_vertex_vertex_angles==21 && "For Trigonal Bipyramidal No. of Angles of Two Vertices with Another Vertex should be 21!!");

		_no_of_features_to_compare=2;
		assert(_no_of_features_to_compare==2 && "For Trigonal Bipyramidal No. of Features to Compare should be 2!!");

		_metal_atom=metal_atom;

		atom all_tribipyr_atoms[_no_of_binding_atoms];
		for (int i=0; i<new_atm.size(); i++){
			all_tribipyr_atoms[i]=new_atm.at(i).get_atom_direct_contact();
		}

		_planar_nos=platom;
		_non_planar_nos=nonplatom;

		int p=0;
		for (auto pl: _planar_nos){
			cout<<"\nPL: "<<pl;
			for (int i=0; i<_no_of_binding_atoms; i++){
				if (pl==all_tribipyr_atoms[i].get_serial_no()){
					_tribipyr_atoms_pl[p]=all_tribipyr_atoms[i];
					_tribipyr_vertices_pl[p]=Point3D(_tribipyr_atoms_pl[p].X(), _tribipyr_atoms_pl[p].Y(), _tribipyr_atoms_pl[p].Z());
					cout<<"\nSr. No.: "<<_tribipyr_atoms_pl[p].get_serial_no()<<"; Tri Bi Pyr Pl: "<<_tribipyr_atoms_pl[p].X()<<", "<<_tribipyr_atoms_pl[p].Y()<<", "<<_tribipyr_atoms_pl[p].Z();
					p++;
				}
			}
		}

		p=0;
		for (auto npl: _non_planar_nos){
			cout<<"\nNPL: "<<npl;
			for (int i=0; i<_no_of_binding_atoms; i++){
				if (npl==all_tribipyr_atoms[i].get_serial_no()){
					_tribipyr_atoms_npl[p]=all_tribipyr_atoms[i];
					_tribipyr_vertices_npl[p]=Point3D(_tribipyr_atoms_npl[p].X(), _tribipyr_atoms_npl[p].Y(), _tribipyr_atoms_npl[p].Z());
					cout<<"\nSr. No.: "<<_tribipyr_atoms_npl[p].get_serial_no()<<"; Tri Bi Pyr Npl: "<<_tribipyr_atoms_npl[p].X()<<", "<<_tribipyr_atoms_npl[p].Y()<<", "<<_tribipyr_atoms_npl[p].Z();
					p++;
				}
			}
		}

		_tribipyr=trigonal_bipyramidal(_tribipyr_vertices_pl, _tribipyr_vertices_npl);

		print(_metal_atom, _tribipyr_atoms_pl, _no_of_planar_atoms);
		print(_tribipyr_atoms_npl, _no_of_non_planar_atoms);
	}

	tribipyrMetalicStructure tribipyrMetalicStructure::apply_transformations(){
		cout<<"\n Original Vertices Planar";
		print(_tribipyr_vertices_pl, _no_of_planar_atoms);

		cout<<"\n Original Vertices Non-Planar";
		print(_tribipyr_vertices_npl, _no_of_non_planar_atoms);

		for (int i=0; i<3; i++){
			_pl_vertex_order[i]=_tribipyr.get_pl_vertex_order()[i];
		}

		for (int i=0; i<2; i++){
			_npl_vertex_order[i]=_tribipyr.get_npl_vertex_order()[i];
		}

		find_largest_side_and_update_vertex_order();
		tribipyrMetalicStructure tms=make_transformation();

		cout<<"\nVERTEX ORDER Planar After Orientation: "<<_pl_vertex_order[0]<<"  "<<_pl_vertex_order[1]<<"  "<<_pl_vertex_order[2];
		cout<<"\nVERTEX ORDER Non-planar After Orientation: "<<_npl_vertex_order[0]<<"  "<<_npl_vertex_order[1];

		return tms;
	}

	atom tribipyrMetalicStructure::get_metal(){
		return _metal_atom;
	}

	trigonal_bipyramidal tribipyrMetalicStructure::get_tribipyr(){
		return _tribipyr;
	}

	Point3D* tribipyrMetalicStructure::get_tribipyr_vertices_pl(){
		return _tribipyr_vertices_pl;
	}

	Point3D* tribipyrMetalicStructure::get_tribipyr_vertices_npl(){
		return _tribipyr_vertices_npl;
	}

	atom* tribipyrMetalicStructure::get_tribipyr_atoms_pl(){
		return _tribipyr_atoms_pl;
	}

	atom* tribipyrMetalicStructure::get_tribipyr_atoms_npl(){
		return _tribipyr_atoms_npl;
	}

	tetrahedron tribipyrMetalicStructure::get_upper_tetra(){
		return _upper_tetra;
	}

	atom* tribipyrMetalicStructure::get_upper_tetra_atoms(){
		return _upper_tetra_atoms;
	}

	atom* tribipyrMetalicStructure::get_lower_tetra_atoms(){
		return _lower_tetra_atoms;
	}

	tetrahedron tribipyrMetalicStructure::get_lower_tetra(){
		return _lower_tetra;
	}

	int tribipyrMetalicStructure::get_no_of_binding_atoms(){
		return _no_of_binding_atoms;
	}

	double* tribipyrMetalicStructure::get_side_deviations(){
		return _devt_dist_vert_vert;
	}

	double* tribipyrMetalicStructure::get_angle_deviations_vertex_vertex(){
		return _devt_angle_vertex_vertex;
	}

	double* tribipyrMetalicStructure::get_mse(){
		return _mse;
	}

	string tribipyrMetalicStructure::get_structure_size(){
		return _structure_size;
	}

	void tribipyrMetalicStructure::comparison(tribipyrMetalicStructure & comp, string out_comparable_name){
		double dist_stand_planar[3][3];
		for (int i=0; i<3; i++){
			for (int j=0; j<3; j++){
				dist_stand_planar[i][j]=_tribipyr.get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		double dist_stand_non_planar[6][3];
		for (int i=0; i<6; i++){
			for (int j=0; j<3; j++){
				dist_stand_non_planar[i][j]=_tribipyr.get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		Point3D* comp_vertices_pl=comp.get_tribipyr_vertices_pl();
		Point3D* comp_vertices_npl=comp.get_tribipyr_vertices_npl();
		measure_direction_cosines();
		double proj_pl_pl[_no_of_planar_atoms];
		double proj_npl_pl[6];

		int k=0;

		for (int i=0; i<(_no_of_planar_atoms-1); i++){
			for (int j=i+1; j<_no_of_planar_atoms; j++){
				double prx=	_direction_cosines_pl[k][0]*(comp_vertices_pl[j].X()-comp_vertices_pl[i].X());
				double pry= _direction_cosines_pl[k][1]*(comp_vertices_pl[j].Y()-comp_vertices_pl[i].Y());
				double prz=	_direction_cosines_pl[k][2]*(comp_vertices_pl[j].Z()-comp_vertices_pl[i].Z());

				proj_pl_pl[k]=(prx+pry+prz);

				k++;
			}
		}

		k=0;
		for (int i=0; i<6; i++){
			for (int j=0; j<_no_of_planar_atoms; j++){
				double prx=	_direction_cosines_npl[k][0]*(comp_vertices_pl[j].X()-comp_vertices_npl[i].X());
				double pry= _direction_cosines_npl[k][1]*(comp_vertices_pl[j].Y()-comp_vertices_npl[i].Y());
				double prz=	_direction_cosines_npl[k][2]*(comp_vertices_pl[j].Z()-comp_vertices_npl[i].Z());

				proj_npl_pl[k]=(prx+pry+prz);

				k++;
			}
		}

		double dist_dev_vert_vert_pl_pl[_no_of_planar_atoms];
		double dist_dev_vert_vert_npl_pl[6];

		FILE* outf=fopen(out_comparable_name.c_str(), "w");

		cout<<"\n\nComparison Dist:\n";
		for (int i=0; i<_no_of_planar_atoms; i++){
			double diff=dist_stand_planar[i][2]-proj_pl_pl[i];
			dist_dev_vert_vert_pl_pl[i]=diff;
			cout<<"Side "<<(i+1)<<": "<<dist_dev_vert_vert_pl_pl[i]<<endl;
		}

		for (int i=0; i<6; i++){
			double diff=dist_stand_non_planar[i][2]-proj_npl_pl[i];
			dist_dev_vert_vert_npl_pl[i]=diff;
			cout<<"Side "<<(i+1)<<": "<<dist_dev_vert_vert_npl_pl[i]<<endl;
		}

		double sum1=0.0;
		for (int p=0; p<_no_of_planar_atoms; p++){
			sum1+=(dist_dev_vert_vert_pl_pl[p] * dist_dev_vert_vert_pl_pl[p]);
		}

		double sum2=0.0;
		for (int p=0; p<6; p++){
			sum2+=(dist_dev_vert_vert_npl_pl[p] * dist_dev_vert_vert_npl_pl[p]);
		}

		double mse_dist[2];
		mse_dist[0]=(sum1/_no_of_planar_atoms);
		mse_dist[1]=(sum2/6);
		double mse_distance=mse_dist[0]+mse_dist[1];
		cout<<"\n MSE Distance: "<<mse_distance;
		fprintf(outf, "\nOriginal MSE Dist: %8.3lf", mse_distance);

		cout<<"\n\nComparison Angle:\n";
		double angle_stand_planar[3][4];
		for (int i=0; i<3; i++){
			for (int j=0; j<4; j++){
				angle_stand_planar[i][j]=_tribipyr.get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		double angle_stand_planar_non_planar[18][4];
		for (int i=0; i<18; i++){
			for (int j=0; j<4; j++){
				angle_stand_planar_non_planar[i][j]=_tribipyr.get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}


		double angle_comp_planar[3][4];
		for (int i=0; i<3; i++){
			for (int j=0; j<4; j++){
				angle_comp_planar[i][j]=comp.get_tribipyr().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		double angle_comp_planar_non_planar[18][4];
		for (int i=0; i<18; i++){
			for (int j=0; j<4; j++){
				angle_comp_planar_non_planar[i][j]=comp.get_tribipyr().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		double angle_dev_vert_vert_pl_pl[3];
		double angle_dev_vert_vert_npl_pl[18];

		for (int p=0; p<_no_of_planar_atoms; p++){
			angle_dev_vert_vert_pl_pl[p]=angle_stand_planar[p][3]-angle_comp_planar[p][3];
			cout<<"\nAngle "<<(p+1)<<": "<<angle_dev_vert_vert_pl_pl[p];
		}

		for (int p=0; p<18; p++){
			angle_dev_vert_vert_npl_pl[p]=angle_stand_planar_non_planar[p][3]-angle_comp_planar_non_planar[p][3];
			cout<<"\nAngle "<<(p+1)<<": "<<angle_dev_vert_vert_npl_pl[p];
		}

		sum1=0.0;
		for (int p=0; p<_no_of_planar_atoms; p++){
			sum1+=(angle_dev_vert_vert_pl_pl[p] * angle_dev_vert_vert_pl_pl[p]);
		}

		sum2=0.0;
		for (int p=0; p<18; p++){
			sum2+=(angle_dev_vert_vert_npl_pl[p] * angle_dev_vert_vert_npl_pl[p]);
		}

		double mse_angle[2];
		mse_angle[0]=(sum1/_no_of_planar_atoms);
		mse_angle[1]=(sum2/18);
		double mse_ang=mse_angle[0]+mse_angle[1];
		cout<<"\n MSE Angle: "<<mse_ang;
		fprintf(outf, "\nMSE Angle: %8.3lf", mse_ang);

		fclose(outf);
	}

	void tribipyrMetalicStructure:: find_largest_side_and_update_vertex_order(){
		double dist[3][3];
		for (int i=0; i<3; i++){
			for (int j=0; j<3; j++){
				dist[i][j]=_tribipyr.get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		int* _initial_vertex_order=new int[3];

		for (int i=0; i<3; i++){
			_initial_vertex_order[i]=_tribipyr.get_pl_vertex_order()[i];
		}

		cout<<"\nInitial V. Order: "<<_pl_vertex_order[0]<<", "<<_pl_vertex_order[1]<<", "<<_pl_vertex_order[2];
		if ((dist[0][2]>=dist[1][2]) && (dist[0][2]>=dist[2][2])){
			_pl_vertex_order[0]=_initial_vertex_order[0];
			_pl_vertex_order[1]=_initial_vertex_order[1];
			_pl_vertex_order[2]=_initial_vertex_order[2];
		}else if ((dist[1][2]>=dist[0][2]) && (dist[1][2]>=dist[2][2])){
			_pl_vertex_order[0]=_initial_vertex_order[1];
			_pl_vertex_order[1]=_initial_vertex_order[2];
			_pl_vertex_order[2]=_initial_vertex_order[0];
		}else if ((dist[2][2]>=dist[0][2]) && (dist[2][2]>=dist[1][2])){
			_pl_vertex_order[0]=_initial_vertex_order[2];
			_pl_vertex_order[1]=_initial_vertex_order[0];
			_pl_vertex_order[2]=_initial_vertex_order[1];
		}

		delete [] _initial_vertex_order;

		cout<<"\nLar Size Up V Order: "<<_pl_vertex_order[0]<<", "<<_pl_vertex_order[1]<<", "<<_pl_vertex_order[2];
	}

	tribipyrMetalicStructure tribipyrMetalicStructure::make_transformation(){
		Point3D* upper_tetra_vertices=new Point3D[4];
		atom upper_tetra_atoms[4];
		Point3D* lower_tetra_vertices=new Point3D[4];
		atom lower_tetra_atoms[4];

		upper_tetra_vertices[0]=_tribipyr_vertices_pl[_pl_vertex_order[0]];
		upper_tetra_vertices[1]=_tribipyr_vertices_pl[_pl_vertex_order[1]];
		upper_tetra_vertices[2]=_tribipyr_vertices_pl[_pl_vertex_order[2]];
		upper_tetra_vertices[3]=_tribipyr_vertices_npl[0];
		upper_tetra_atoms[0]=_tribipyr_atoms_pl[_pl_vertex_order[0]];
		upper_tetra_atoms[1]=_tribipyr_atoms_pl[_pl_vertex_order[1]];
		upper_tetra_atoms[2]=_tribipyr_atoms_pl[_pl_vertex_order[2]];
		upper_tetra_atoms[3]=_tribipyr_atoms_npl[0];

		tetraMetalicStructure upper_struct=tetraMetalicStructure(upper_tetra_atoms, upper_tetra_vertices, _metal_atom, 0, 1, 2, 3);
		delete [] upper_tetra_vertices;

		tetraMetalicStructure  upper_transformed=upper_struct.transform_tetrahedron();

		lower_tetra_vertices[0]=_tribipyr_vertices_pl[_pl_vertex_order[0]];
		lower_tetra_vertices[1]=_tribipyr_vertices_pl[_pl_vertex_order[1]];
		lower_tetra_vertices[2]=_tribipyr_vertices_pl[_pl_vertex_order[2]];
		lower_tetra_vertices[3]=_tribipyr_vertices_npl[1];
		lower_tetra_atoms[0]=_tribipyr_atoms_pl[_pl_vertex_order[0]];
		lower_tetra_atoms[1]=_tribipyr_atoms_pl[_pl_vertex_order[1]];
		lower_tetra_atoms[2]=_tribipyr_atoms_pl[_pl_vertex_order[2]];
		lower_tetra_atoms[3]=_tribipyr_atoms_npl[1];

		tetraMetalicStructure lower_struct=tetraMetalicStructure(lower_tetra_atoms, lower_tetra_vertices, _metal_atom, 0, 1, 2, 3);
		delete [] lower_tetra_vertices;

		tetraMetalicStructure  lower_transformed=lower_struct.transform_tetrahedron();

		tribipyrMetalicStructure tps=tribipyrMetalicStructure(upper_transformed.get_tetrahedron(), upper_tetra_atoms, lower_transformed.get_tetrahedron(), lower_tetra_atoms, _metal_atom);

		return tps;
	}

	void tribipyrMetalicStructure::measure_direction_cosines(){
		int k=0;
		double x, y, z, l, m, n, r;
		for (int i=0; i<_no_of_planar_atoms-1; i++){
			for (int j=i+1; j<_no_of_planar_atoms; j++){
				x=_tribipyr_vertices_pl[j].X() - _tribipyr_vertices_pl[i].X();
				y=_tribipyr_vertices_pl[j].Y() - _tribipyr_vertices_pl[i].Y();
				z=_tribipyr_vertices_pl[j].Z() - _tribipyr_vertices_pl[i].Z();
				r=sqrt(x*x+y*y+z*z);

				l=x/r;
				m=y/r;
				n=z/r;

				_direction_cosines_pl[k][0]=l;
				_direction_cosines_pl[k][1]=m;
				_direction_cosines_pl[k][2]=n;
				k++;
			}
		}

		k=0;
		for (int i=0; i<_no_of_non_planar_atoms; i++){
			for (int j=0; j<_no_of_planar_atoms; j++){

				x=_tribipyr_vertices_pl[j].X() - _tribipyr_vertices_npl[i].X();
				y=_tribipyr_vertices_pl[j].Y() - _tribipyr_vertices_npl[i].Y();
				z=_tribipyr_vertices_pl[j].Z() - _tribipyr_vertices_npl[i].Z();
				r=sqrt(x*x+y*y+z*z);

				l=x/r;
				m=y/r;
				n=z/r;

				_direction_cosines_npl[k][0]=l;
				_direction_cosines_npl[k][1]=m;
				_direction_cosines_npl[k][2]=n;
				k++;
			}
		}
	}

	void tribipyrMetalicStructure::write_atoms_pdb_format(string fname){
		FILE* outf=fopen(fname.c_str(), "w");

		for (int i=0; i<_no_of_planar_atoms; i++){
			fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _tribipyr_atoms_pl[i].get_atom_type().c_str(), _tribipyr_atoms_pl[i].get_serial_no(), _tribipyr_atoms_pl[i].get_atom_name().c_str(), _tribipyr_atoms_pl[i].get_residue_name().c_str(), _tribipyr_atoms_pl[i].get_chain_name().c_str(), _tribipyr_atoms_pl[i].get_residue_no(), _tribipyr_atoms_pl[i].X(), _tribipyr_atoms_pl[i].Y(), _tribipyr_atoms_pl[i].Z(), _tribipyr_atoms_pl[i].get_occupancy(), _tribipyr_atoms_pl[i].get_temp_factor(), _tribipyr_atoms_pl[i].get_atom_symbol().c_str());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _tribipyr_atoms_npl[i].get_atom_type().c_str(), _tribipyr_atoms_npl[i].get_serial_no(), _tribipyr_atoms_npl[i].get_atom_name().c_str(), _tribipyr_atoms_npl[i].get_residue_name().c_str(), _tribipyr_atoms_npl[i].get_chain_name().c_str(), _tribipyr_atoms_npl[i].get_residue_no(), _tribipyr_atoms_npl[i].X(), _tribipyr_atoms_npl[i].Y(), _tribipyr_atoms_npl[i].Z(), _tribipyr_atoms_npl[i].get_occupancy(), _tribipyr_atoms_npl[i].get_temp_factor(), _tribipyr_atoms_npl[i].get_atom_symbol().c_str());
		}

		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[0].get_serial_no(), _tribipyr_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[1].get_serial_no(), _tribipyr_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[1].get_serial_no(), _tribipyr_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[2].get_serial_no(), _tribipyr_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[2].get_serial_no(), _tribipyr_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[0].get_serial_no(), _tribipyr_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_npl[0].get_serial_no(), _tribipyr_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[0].get_serial_no(), _tribipyr_atoms_npl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_npl[0].get_serial_no(), _tribipyr_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[1].get_serial_no(), _tribipyr_atoms_npl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_npl[0].get_serial_no(), _tribipyr_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[2].get_serial_no(), _tribipyr_atoms_npl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_npl[1].get_serial_no(), _tribipyr_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[0].get_serial_no(), _tribipyr_atoms_npl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_npl[1].get_serial_no(), _tribipyr_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[1].get_serial_no(), _tribipyr_atoms_npl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_npl[1].get_serial_no(), _tribipyr_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tribipyr_atoms_pl[2].get_serial_no(), _tribipyr_atoms_npl[1].get_serial_no());

		fclose(outf);
	}

	void tribipyrMetalicStructure::write_atoms_combined_tex_format_upper(tribipyrMetalicStructure & comp, string fname){
		string out_file_name=fname+".tex";
		FILE* outf=fopen(out_file_name.c_str(), "w");
		fprintf(outf, "\\documentclass[a4paper]{article}\n");
		fprintf(outf, "\\usepackage[utf8]{inputenc}\n");
		fprintf(outf, "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}");
		fprintf(outf, "\\usepackage{chemfig}\n");
		fprintf(outf, "\\begin{document}\n");
		fprintf(outf, "\\setbondoffset{3pt}\n");
		string st_upper[4], st1, str;
		double bond_angle_standard_upper[4];
		double bond_angle_comparable_upper[4];
		//atom** _comp_tribipyr_atoms=comp.get_tribipyr_atoms();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		tetrahedron upper_comp=comp.get_upper_tetra();

		Point3D* upper_vertices_standard=_upper_tetra.get_vertices();
		Point3D* upper_vertices_comparable=upper_comp.get_vertices();

		int* stand_vertex_order_upper=_upper_tetra.get_vertex_order();
		int* comp_vertex_order_upper=upper_comp.get_vertex_order();

		double bond_len_standard_upper=_upper_tetra.get_circum_radius();
		double bond_len_comparable_upper=upper_comp.get_circum_radius();

		string adjust_stand_upper="";
		string adjust_comp_upper="";

		if (bond_len_standard_upper>=8.0){
			bond_len_standard_upper=3.5;
			adjust_stand_upper="adjusted";
		}

		if (bond_len_comparable_upper>=8.0){
			bond_len_comparable_upper=3.5;
			adjust_comp_upper="adjusted";
		}

		str="\\center\\chemname[30pt]{\\chemfig{\\chembelow[4pt]{Circumcentre}{Origin}";
		fprintf(outf, "%s", str.c_str());

		st_upper[0]="-";
		st_upper[1]="-";
		st_upper[2]="-";
		st_upper[3]="<";

		//Compute bond angles for Standard Upper Tetrahedron
		double ang1=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_standard[stand_vertex_order_upper[0]], &upper_vertices_standard[stand_vertex_order_upper[1]]);
		bond_angle_standard_upper[0]=180.0+ang1;
		double ang2=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_standard[stand_vertex_order_upper[1]], &upper_vertices_standard[stand_vertex_order_upper[0]]);
		bond_angle_standard_upper[1]=-ang2;

		if (upper_vertices_standard[stand_vertex_order_upper[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&upper_vertices_standard[stand_vertex_order_upper[1]], &_centre, &upper_vertices_standard[stand_vertex_order_upper[2]]);
			bond_angle_standard_upper[2]=ang-ang2;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&upper_vertices_standard[stand_vertex_order_upper[0]], &_centre, &upper_vertices_standard[stand_vertex_order_upper[2]]);
			bond_angle_standard_upper[2]=ang-ang1;
		}

		if ((upper_vertices_standard[stand_vertex_order_upper[3]].X()>=0.0) && (upper_vertices_standard[stand_vertex_order_upper[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_standard[stand_vertex_order_upper[3]], &upper_vertices_standard[stand_vertex_order_upper[2]]);
			bond_angle_standard_upper[3]=90-ang/2.0;
		}else if ((upper_vertices_standard[stand_vertex_order_upper[3]].X()<0.0) && (upper_vertices_standard[stand_vertex_order_upper[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_standard[stand_vertex_order_upper[3]], &upper_vertices_standard[stand_vertex_order_upper[2]]);
			bond_angle_standard_upper[3]=90+ang/2.0;
		}else if ((upper_vertices_standard[stand_vertex_order_upper[3]].X()<0.0) && (upper_vertices_standard[stand_vertex_order_upper[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_standard[stand_vertex_order_upper[3]], &upper_vertices_standard[stand_vertex_order_upper[0]]);
			bond_angle_standard_upper[3]=180+ang/2.0;
		}else  if ((upper_vertices_standard[stand_vertex_order_upper[3]].X()>=0.0) && (upper_vertices_standard[stand_vertex_order_upper[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_standard[stand_vertex_order_upper[1]], &upper_vertices_standard[stand_vertex_order_upper[3]]);
			bond_angle_standard_upper[3]=360-ang/2.0;
		}

		//Compute bond angles for Comparable Upper Tetrahedron
		double ang3=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_comparable[comp_vertex_order_upper[0]], &upper_vertices_comparable[comp_vertex_order_upper[1]]);
		bond_angle_comparable_upper[0]=180.0+ang3;
		double ang4=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_comparable[comp_vertex_order_upper[1]], &upper_vertices_comparable[comp_vertex_order_upper[0]]);
		bond_angle_comparable_upper[1]=-ang4;

		if (upper_vertices_comparable[comp_vertex_order_upper[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&upper_vertices_comparable[comp_vertex_order_upper[1]], &_centre, &upper_vertices_comparable[comp_vertex_order_upper[2]]);
			bond_angle_comparable_upper[2]=ang-ang4;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&upper_vertices_comparable[comp_vertex_order_upper[0]], &_centre, &upper_vertices_comparable[comp_vertex_order_upper[2]]);
			bond_angle_comparable_upper[2]=ang-ang3;
		}

		if ((upper_vertices_comparable[comp_vertex_order_upper[3]].X()>=0.0) && (upper_vertices_comparable[comp_vertex_order_upper[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_comparable[comp_vertex_order_upper[3]], &upper_vertices_comparable[comp_vertex_order_upper[2]]);
			bond_angle_comparable_upper[3]=90-ang/2.0;
		}else if ((upper_vertices_comparable[comp_vertex_order_upper[3]].X()<0.0) && (upper_vertices_comparable[comp_vertex_order_upper[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_comparable[comp_vertex_order_upper[3]], &upper_vertices_comparable[comp_vertex_order_upper[2]]);
			bond_angle_comparable_upper[3]=90+ang/2.0;
		}else if ((upper_vertices_comparable[comp_vertex_order_upper[3]].X()<0.0) && (upper_vertices_comparable[comp_vertex_order_upper[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_comparable[comp_vertex_order_upper[3]], &upper_vertices_comparable[comp_vertex_order_upper[0]]);
			bond_angle_comparable_upper[3]=180+ang/2.0;
		}else  if ((upper_vertices_comparable[comp_vertex_order_upper[3]].X()>=0.0) && (upper_vertices_comparable[comp_vertex_order_upper[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &upper_vertices_comparable[comp_vertex_order_upper[1]], &upper_vertices_comparable[comp_vertex_order_upper[3]]);
			bond_angle_comparable_upper[3]=360-ang/2.0;
		}

		//string main_string="";
		for (int p=0; p<4; p++){
			//Write Standard Tetrahedron Upper Vertices
			str="";

			if (_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_name()=="HOH"){
				st1=_upper_tetra_atoms[stand_vertex_order_upper[p]].get_chain_name()+":"+_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_name()+"_{"+to_string(_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_upper_tetra_atoms[stand_vertex_order_upper[p]].get_atom_name()=="O"){
					st1=_upper_tetra_atoms[stand_vertex_order_upper[p]].get_chain_name()+":"+_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_name()+"_{"+to_string(_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_tetra_atoms[stand_vertex_order_upper[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_upper_tetra_atoms[stand_vertex_order_upper[p]].get_atom_name()=="OXT"){
					st1=_upper_tetra_atoms[stand_vertex_order_upper[p]].get_chain_name()+":"+_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_name()+"_{"+to_string(_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_tetra_atoms[stand_vertex_order_upper[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_upper_tetra_atoms[stand_vertex_order_upper[p]].get_chain_name()+":"+_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_name()+"_{"+to_string(_upper_tetra_atoms[stand_vertex_order_upper[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_tetra_atoms[stand_vertex_order_upper[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_upper[p].c_str(), bond_angle_standard_upper[p], bond_len_standard_upper, str.c_str());

			//Write Comparable Tetrahedron Upper Vertices
			str="";
			atom* comp_upper_tetra_atoms=comp.get_upper_tetra_atoms();
			if (comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_name()=="HOH"){
				st1=comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_chain_name()+":"+comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_name()+"_{"+to_string(comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{blue}H_2O}{\\color{blue}"+st1+"}";
			}else{
				if (comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_atom_name()=="O"){
					st1=comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_chain_name()+":"+comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_name()+"_{"+to_string(comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_atom_symbol()+"1}{\\color{blue}"+st1+"}";
				}else if (comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_atom_name()=="OXT"){
					st1=comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_chain_name()+":"+comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_name()+"_{"+to_string(comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_atom_symbol()+"2}{\\color{blue}"+st1+"}";
				}else{
					st1=comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_chain_name()+":"+comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_name()+"_{"+to_string(comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_upper_tetra_atoms[comp_vertex_order_upper[p]].get_atom_name()+"}{\\color{blue}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,blue]%s)", st_upper[p].c_str(), bond_angle_comparable_upper[p], bond_len_comparable_upper, str.c_str());
		}

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms Upper}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Trigonal Bipyramidal Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Trigonal Bipyramidal Site to Any Atom\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); If Bond Length Exceeds Page Width, It Is Set Equal to the Page Width\\end{flushleft}";

		fprintf(outf, "%s", str.c_str());

		if (adjust_stand_upper!=""){
			fprintf(outf, "\n\\vspace{2ex}Here Bond Length for Standard Site is Adjusted to the Page Width");
		}

		if (adjust_comp_upper!=""){
			fprintf(outf, "\n\\vspace{2ex}Here Bond Length for Comparable Site is Adjusted to the Page Width");
		}

		fprintf(outf, "\\end{document}");
		fclose(outf);

		string mycommand = "latex " + out_file_name;
		system(mycommand.c_str());

		mycommand = "dvipdf " + fname + ".dvi";
		system(mycommand.c_str());

		mycommand = "rm " + fname + ".dvi";
		system(mycommand.c_str());

		mycommand = "rm " + fname + ".aux";
		system(mycommand.c_str());

		mycommand = "rm " + fname + ".log";
		system(mycommand.c_str());

		mycommand = "rm " + fname + ".tex";
		system(mycommand.c_str());
	}

	void tribipyrMetalicStructure::write_atoms_combined_tex_format_lower(tribipyrMetalicStructure & comp, string fname){
		string out_file_name=fname+".tex";
		FILE* outf=fopen(out_file_name.c_str(), "w");
		fprintf(outf, "\\documentclass[a4paper]{article}\n");
		fprintf(outf, "\\usepackage[utf8]{inputenc}\n");
		fprintf(outf, "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}");
		fprintf(outf, "\\usepackage{chemfig}\n");
		fprintf(outf, "\\begin{document}\n");
		fprintf(outf, "\\setbondoffset{3pt}\n");
		string st_lower[4], st1, str;
		double bond_angle_standard_lower[4];
		double bond_angle_comparable_lower[4];
		//atom** _comp_tribipyr_atoms=comp.get_tribipyr_atoms();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		tetrahedron lower_comp=comp.get_lower_tetra();

		Point3D* lower_vertices_standard=_lower_tetra.get_vertices();
		Point3D* lower_vertices_comparable=lower_comp.get_vertices();

		int* stand_vertex_order_lower=_lower_tetra.get_vertex_order();
		int* comp_vertex_order_lower=lower_comp.get_vertex_order();

		double bond_len_standard_lower=_lower_tetra.get_circum_radius();
		double bond_len_comparable_lower=lower_comp.get_circum_radius();

		string adjust_stand_lower="";
		string adjust_comp_lower="";

		if (bond_len_standard_lower>=8.0){
			bond_len_standard_lower=3.5;
			adjust_stand_lower="adjusted";
		}

		if (bond_len_comparable_lower>=8.0){
			bond_len_comparable_lower=3.5;
			adjust_comp_lower="adjusted";
		}

		str="\\center\\chemname[30pt]{\\chemfig{\\chembelow[4pt]{Circumcentre}{Origin}";
		fprintf(outf, "%s", str.c_str());

		st_lower[0]="-";
		st_lower[1]="-";
		st_lower[2]="-";
		st_lower[3]="<:";

		//Compute bond angles for Standard Lower Tetrahedron
		double ang11=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_standard[stand_vertex_order_lower[0]], &lower_vertices_standard[stand_vertex_order_lower[1]]);
		bond_angle_standard_lower[0]=180.0+ang11;
		double ang12=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_standard[stand_vertex_order_lower[1]], &lower_vertices_standard[stand_vertex_order_lower[0]]);
		bond_angle_standard_lower[1]=-ang12;

		if (lower_vertices_standard[stand_vertex_order_lower[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&lower_vertices_standard[stand_vertex_order_lower[1]], &_centre, &lower_vertices_standard[stand_vertex_order_lower[2]]);
			bond_angle_standard_lower[2]=ang-ang12;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&lower_vertices_standard[stand_vertex_order_lower[0]], &_centre, &lower_vertices_standard[stand_vertex_order_lower[2]]);
			bond_angle_standard_lower[2]=ang-ang11;
		}

		if ((lower_vertices_standard[stand_vertex_order_lower[3]].X()>=0.0) && (lower_vertices_standard[stand_vertex_order_lower[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_standard[stand_vertex_order_lower[3]], &lower_vertices_standard[stand_vertex_order_lower[2]]);
			bond_angle_standard_lower[3]=90-ang/2.0;
		}else if ((lower_vertices_standard[stand_vertex_order_lower[3]].X()<0.0) && (lower_vertices_standard[stand_vertex_order_lower[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_standard[stand_vertex_order_lower[3]], &lower_vertices_standard[stand_vertex_order_lower[2]]);
			bond_angle_standard_lower[3]=90+ang/2.0;
		}else if ((lower_vertices_standard[stand_vertex_order_lower[3]].X()<0.0) && (lower_vertices_standard[stand_vertex_order_lower[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_standard[stand_vertex_order_lower[3]], &lower_vertices_standard[stand_vertex_order_lower[0]]);
			bond_angle_standard_lower[3]=180+ang/2.0;
		}else  if ((lower_vertices_standard[stand_vertex_order_lower[3]].X()>=0.0) && (lower_vertices_standard[stand_vertex_order_lower[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_standard[stand_vertex_order_lower[1]], &lower_vertices_standard[stand_vertex_order_lower[3]]);
			bond_angle_standard_lower[3]=360-ang/2.0;
		}

		//Compute bond angles for Comparable Lower Tetrahedron
		double ang13=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_comparable[comp_vertex_order_lower[0]], &lower_vertices_comparable[comp_vertex_order_lower[1]]);
		bond_angle_comparable_lower[0]=180.0+ang13;
		double ang14=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_comparable[comp_vertex_order_lower[1]], &lower_vertices_comparable[comp_vertex_order_lower[0]]);
		bond_angle_comparable_lower[1]=-ang14;

		if (lower_vertices_comparable[comp_vertex_order_lower[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&lower_vertices_comparable[comp_vertex_order_lower[1]], &_centre, &lower_vertices_comparable[comp_vertex_order_lower[2]]);
			bond_angle_comparable_lower[2]=ang-ang14;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&lower_vertices_comparable[comp_vertex_order_lower[0]], &_centre, &lower_vertices_comparable[comp_vertex_order_lower[2]]);
			bond_angle_comparable_lower[2]=ang-ang13;
		}

		if ((lower_vertices_comparable[comp_vertex_order_lower[3]].X()>=0.0) && (lower_vertices_comparable[comp_vertex_order_lower[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_comparable[comp_vertex_order_lower[3]], &lower_vertices_comparable[comp_vertex_order_lower[2]]);
			bond_angle_comparable_lower[3]=90-ang/2.0;
		}else if ((lower_vertices_comparable[comp_vertex_order_lower[3]].X()<0.0) && (lower_vertices_comparable[comp_vertex_order_lower[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_comparable[comp_vertex_order_lower[3]], &lower_vertices_comparable[comp_vertex_order_lower[2]]);
			bond_angle_comparable_lower[3]=90+ang/2.0;
		}else if ((lower_vertices_comparable[comp_vertex_order_lower[3]].X()<0.0) && (lower_vertices_comparable[comp_vertex_order_lower[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_comparable[comp_vertex_order_lower[3]], &lower_vertices_comparable[comp_vertex_order_lower[0]]);
			bond_angle_comparable_lower[3]=180+ang/2.0;
		}else  if ((lower_vertices_comparable[comp_vertex_order_lower[3]].X()>=0.0) && (lower_vertices_comparable[comp_vertex_order_lower[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &lower_vertices_comparable[comp_vertex_order_lower[1]], &lower_vertices_comparable[comp_vertex_order_lower[3]]);
			bond_angle_comparable_lower[3]=360-ang/2.0;
		}

		for (int p=0; p<4; p++){
			//Write Standard Tetrahedron Lower Vertices
			str="";

			if (_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_name()=="HOH"){
				st1=_lower_tetra_atoms[stand_vertex_order_lower[p]].get_chain_name()+":"+_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_name()+"_{"+to_string(_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_lower_tetra_atoms[stand_vertex_order_lower[p]].get_atom_name()=="O"){
					st1=_lower_tetra_atoms[stand_vertex_order_lower[p]].get_chain_name()+":"+_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_name()+"_{"+to_string(_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_tetra_atoms[stand_vertex_order_lower[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_lower_tetra_atoms[stand_vertex_order_lower[p]].get_atom_name()=="OXT"){
					st1=_lower_tetra_atoms[stand_vertex_order_lower[p]].get_chain_name()+":"+_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_name()+"_{"+to_string(_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_tetra_atoms[stand_vertex_order_lower[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_lower_tetra_atoms[stand_vertex_order_lower[p]].get_chain_name()+":"+_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_name()+"_{"+to_string(_lower_tetra_atoms[stand_vertex_order_lower[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_tetra_atoms[stand_vertex_order_lower[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_lower[p].c_str(), bond_angle_standard_lower[p], bond_len_standard_lower, str.c_str());

			//Write Comparable Tetrahedron Lower Vertices
			str="";
			atom* comp_lower_tetra_atoms=comp.get_lower_tetra_atoms();
			if (comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_name()=="HOH"){
				st1=comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_chain_name()+":"+comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_name()+"_{"+to_string(comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{blue}H_2O}{\\color{blue}"+st1+"}";
			}else{
				if (comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_atom_name()=="O"){
					st1=comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_chain_name()+":"+comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_name()+"_{"+to_string(comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_atom_symbol()+"1}{\\color{blue}"+st1+"}";
				}else if (comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_atom_name()=="OXT"){
					st1=comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_chain_name()+":"+comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_name()+"_{"+to_string(comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_atom_symbol()+"2}{\\color{blue}"+st1+"}";
				}else{
					st1=comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_chain_name()+":"+comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_name()+"_{"+to_string(comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_lower_tetra_atoms[comp_vertex_order_lower[p]].get_atom_name()+"}{\\color{blue}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,blue]%s)", st_lower[p].c_str(), bond_angle_comparable_lower[p], bond_len_comparable_lower, str.c_str());
		}

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms Lower}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Trigonal Bipyramidal Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Trigonal Bipyramidal Site to Any Atom\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); If Bond Length Exceeds Page Width, It Is Set Equal to the Page Width\\end{flushleft}";

		fprintf(outf, "%s", str.c_str());

		if (adjust_stand_lower!=""){
			fprintf(outf, "\n\\vspace{2ex}Here Bond Length for Standard Site is Adjusted to the Page Width");
		}

		if (adjust_comp_lower!=""){
			fprintf(outf, "\n\\vspace{2ex}Here Bond Length for Comparable Site is Adjusted to the Page Width");
		}

		fprintf(outf, "\\end{document}");
		fclose(outf);

		string mycommand = "latex " + out_file_name;
		system(mycommand.c_str());

		mycommand = "dvipdf " + fname + ".dvi";
		system(mycommand.c_str());

		mycommand = "rm " + fname + ".dvi";
		system(mycommand.c_str());

		mycommand = "rm " + fname + ".aux";
		system(mycommand.c_str());

		mycommand = "rm " + fname + ".log";
		system(mycommand.c_str());

		mycommand = "rm " + fname + ".tex";
		system(mycommand.c_str());
	}

	void tribipyrMetalicStructure::compare(tribipyrMetalicStructure & comp){
		double deviate_dist_upper[6];
		double deviate_dist_lower[6];

		_upper_tetra.deviation_distances(comp.get_upper_tetra(), _upper_tetra.measure_direction_cosines_of_sides(), deviate_dist_upper);
		_lower_tetra.deviation_distances(comp.get_lower_tetra(), _lower_tetra.measure_direction_cosines_of_sides(), deviate_dist_lower);

		for (int p=0; p<6; p++){
			cout<<"\nUpper Dist:"<<deviate_dist_upper[p]<<", Lower Dist: "<<deviate_dist_lower[p];
			_devt_dist_vert_vert_upper[p]=deviate_dist_upper[p];
			_devt_dist_vert_vert_lower[p]=deviate_dist_lower[p];
			_devt_dist_vert_vert[p]=_devt_dist_vert_vert_upper[p];
		}

		int k=6;
		for (int p=0; p<6; p++){
			if (p==0 || p==1 || p==3){
				continue;
			}else{
				_devt_dist_vert_vert[k]=_devt_dist_vert_vert_lower[p];
				k++;
			}
		}

		double deviate_angle_upper[12];
		double deviate_angle_lower[12];

		_upper_tetra.deviation_angles_vertex_vertex(comp.get_upper_tetra(), deviate_angle_upper);
		_lower_tetra.deviation_angles_vertex_vertex(comp.get_lower_tetra(), deviate_angle_lower);

		for (int p=0; p<12; p++){
			_devt_angle_vertex_vertex_upper[p]=deviate_angle_upper[p];
			_devt_angle_vertex_vertex_lower[p]=deviate_angle_lower[p];
			_devt_angle_vertex_vertex[p]=_devt_angle_vertex_vertex_upper[p];
			//cout<<"\nVert-Vert Angle IN COMPARE: "<<(p+1)<<": "<<_devt_angle_vertex_vertex_upper[p]<<", "<<_devt_angle_vertex_vertex_lower[p];
		}

		k=12;
		for (int p=0; p<12; p++){
			if (p==0 || p==3 || p==6){
				continue;
			}else{
				_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_lower[p];
				k++;
			}
		}

		compute_mse();

		_dist_devt_from_cc_upper=(_upper_tetra.get_circum_radius()-comp.get_upper_tetra().get_circum_radius());
		_dist_devt_from_cc_lower=(_lower_tetra.get_circum_radius()-comp.get_lower_tetra().get_circum_radius());

		if ((_dist_devt_from_cc_upper==0.0) && (_dist_devt_from_cc_lower==0.0)){
			_structure_size="same";
		}else if ((_dist_devt_from_cc_upper>0.0) && (_dist_devt_from_cc_lower>0.0)){
			_structure_size="smaller";
		}else if ((_dist_devt_from_cc_upper<0.0) && (_dist_devt_from_cc_lower<0.0)){
			_structure_size="larger";
		}else{
			_structure_size="irregular";
		}

		cout<<"\nCompare End: "<<_structure_size<<", Upper CCRad Stand: "<<_upper_tetra.get_circum_radius()<<", Upper CCRad Comp: "<<comp.get_upper_tetra().get_circum_radius()<<", Diff: "<<_dist_devt_from_cc_upper;
		cout<<"\nCompare End: "<<_structure_size<<", Lower CCRad Stand: "<<_lower_tetra.get_circum_radius()<<", Lower CCRad Comp: "<<comp.get_lower_tetra().get_circum_radius()<<", Diff: "<<_dist_devt_from_cc_lower;
	}

	void tribipyrMetalicStructure::compute_mse(){
		double sum=0.0;
		for (int i=0; i<_no_of_sides; i++){
			sum+=2*(_devt_dist_vert_vert[i]/2.0) * (_devt_dist_vert_vert[i]/2.0);
		}

		_mse[0]=sqrt((sum/(2*_no_of_sides)));
		cout<<"\n MSE Distance: "<<_mse[0];

		sum=0.0;
		for (int i=0; i<_no_of_vertex_vertex_angles; i++){
			sum+=2*(_devt_angle_vertex_vertex[i]/2.0) * (_devt_angle_vertex_vertex[i]/2.0);
		}

		_mse[1]=sqrt((sum/(2*_no_of_vertex_vertex_angles)));
		cout<<"\n MSE Angle (between vertex to vertex): "<<_mse[1];
	}

	void tribipyrMetalicStructure::gen_report_metal_binding_sites(tribipyrMetalicStructure & ori_stand, tribipyrMetalicStructure & ori_comp, tribipyrMetalicStructure & trans_comp, string mt, string standard_site_type, string comparable_site_type, FILE* fp){
		fprintf(fp, "HEADER     PROGRAM NAME: StructureDeviation VERSION: 1.1");
		string str="TITLE      Deviation Report for Computation of Structural Deviations of Two Trigonal Bipyramidal "+mt+" ION Binding Sites";
		fprintf(fp, "\n%s", str.c_str());

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Standard Trigonal Bipyramid:                   ");
		fprintf(fp, "\nREMARK     1  Molecule Type of Standard Trigonal Bipyramid: %s", standard_site_type.c_str());

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Distances of Atoms-Atoms:                         ");

		for (int i=0; i<3; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_stand_pl_pl[i][j]=ori_stand.get_tribipyr().get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_stand.get_tribipyr_atoms_pl(), _vertex_vertex_distances_stand_pl_pl, 3);

		for (int i=0; i<6; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_stand_npl_pl[i][j]=ori_stand.get_tribipyr().get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_stand.get_tribipyr_atoms_pl(), ori_stand.get_tribipyr_atoms_npl(), _vertex_vertex_distances_stand_npl_pl);

		fprintf(fp, "\nREMARK     1  Angles made by two atoms with another atom of Trigonal Bipyramid:");
		for (int i=0; i<3; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_stand_pl_pl[i][j]=ori_stand.get_tribipyr().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_stand.get_tribipyr_atoms_pl(), _vertex_vertex_angles_stand_pl_pl, 3);

		for (int i=0; i<18; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_stand_npl_pl[i][j]=ori_stand.get_tribipyr().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_stand.get_tribipyr_atoms_pl(), ori_stand.get_tribipyr_atoms_npl(), _vertex_vertex_angles_stand_npl_pl);

		fprintf(fp, "\nREMARK     2");
		fprintf(fp, "\nREMARK     2  Comparable Tetragonal Plane:");
		fprintf(fp, "\nREMARK     1  Molecule Type of Comparable Trigonal Bipyramid: %s", comparable_site_type.c_str());

		if (_structure_size=="larger"){
			fprintf(fp, "\nREMARK     2  Size of This Site is Larger than the Standard Site");
		}else if (_structure_size=="smaller"){
			fprintf(fp, "\nREMARK     2  Size of This Site is smaller than the Standard Site");
		}else if(_structure_size=="same"){
			fprintf(fp, "\nREMARK     2  Size of This Site is same as that of the Standard Site");
		}else if(_structure_size=="irregular"){
			fprintf(fp, "\nREMARK     2  Size of This Site is irregular in comparison to the Standard Site");
		}

		fprintf(fp, "\nREMARK     2  Distances of Atoms from  Metal Ion:");

		fprintf(fp, "\nREMARK     2  Distances of Atoms-Atoms:");

		for (int i=0; i<3; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_comp_pl_pl[i][j]=ori_comp.get_tribipyr().get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_comp.get_tribipyr_atoms_pl(), _vertex_vertex_distances_comp_pl_pl, 3);

		for (int i=0; i<6; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_comp_npl_pl[i][j]=ori_comp.get_tribipyr().get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_comp.get_tribipyr_atoms_pl(), ori_comp.get_tribipyr_atoms_npl(), _vertex_vertex_distances_comp_npl_pl);

		fprintf(fp, "\nREMARK     2  Angles made by two atoms with another atom of Tetragonal Plane:");
		for (int i=0; i<3; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_comp_pl_pl[i][j]=ori_comp.get_tribipyr().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_comp.get_tribipyr_atoms_pl(), _vertex_vertex_angles_comp_pl_pl, 3);

		for (int i=0; i<18; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_comp_npl_pl[i][j]=ori_comp.get_tribipyr().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_comp.get_tribipyr_atoms_pl(), ori_comp.get_tribipyr_atoms_npl(), _vertex_vertex_angles_comp_npl_pl);

		fprintf(fp, "\nREMARK     3");
		fprintf(fp, "\nREMARK     3  Deviation Results of Comparable Trigonal Bipyramid from Standard Trigonal Bipyramid:");

		fprintf(fp, "\nREMARK     3  Distance Differences of Atoms-Atoms: ");
		for (int i=0; i<_no_of_sides; i++){
			fprintf(fp, "\nREMARK     3  SIDE %d: %.3lf    ", (i+1), _devt_dist_vert_vert[i]);
		}

		fprintf(fp, "\nREMARK     3  RMSE for Atom-Atom Distance %s: %.3lf   ", mt.c_str(), _mse[0]);

		fprintf(fp, "\nREMARK     3  Differences of Angles Between Two Atoms with Another:");
		for (int i=0; i<_no_of_vertex_vertex_angles; i++){
			fprintf(fp, "\nREMARK     3  Angle %d: %.3lf   ", (i+1), _devt_angle_vertex_vertex[i]);
		}

		fprintf(fp, "\nREMARK     3  RMSE for Angle Differences AA %s: %.3lf   ", mt.c_str(), _mse[1]);
	}

	void tribipyrMetalicStructure::gen_report_vertex_vertex_distances(FILE* fp, atom atm[], double vert_vert_distances[][3], int no_of_sides){
		for (int i=0; i<no_of_sides; i++){
			fprintf(fp, "\nDISTRR %ld %ld %8.3lf", atm[(int) vert_vert_distances[i][0]].get_serial_no(), atm[(int) vert_vert_distances[i][1]].get_serial_no(), vert_vert_distances[i][2]);
		}
	}

	void tribipyrMetalicStructure::gen_report_vertex_vertex_distances(FILE* fp, atom atm_pl[], atom atm_npl[], double vert_vert_distances[][3]){
		int p=0;
		long int l, m;

		for (int p=0; p<6; p++){
			if (vert_vert_distances[p][0]==0){
				l=atm_npl[0].get_serial_no();
			}else if (vert_vert_distances[p][0]==1){
				l=atm_npl[1].get_serial_no();
			}

			if (vert_vert_distances[p][1]==0){
				m=atm_pl[0].get_serial_no();
			}else if (vert_vert_distances[p][1]==1){
				m=atm_pl[1].get_serial_no();
			}else if (vert_vert_distances[p][1]==2){
				m=atm_pl[2].get_serial_no();
			}

			fprintf(fp, "\nDISTRR %ld %ld %8.3lf", l, m, vert_vert_distances[p][2]);
		}
	}

	void tribipyrMetalicStructure::gen_report_vertex_vertex_angles(FILE* fp, atom atm[], double vert_vert_angles[][4], int no_of_angles){
		for (int i=0; i<no_of_angles; i++){
			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", atm[(int) vert_vert_angles[i][0]].get_serial_no(), atm[(int) vert_vert_angles[i][1]].get_serial_no(), atm[(int) vert_vert_angles[i][2]].get_serial_no(), vert_vert_angles[i][3]);
		}
	}

	void tribipyrMetalicStructure::gen_report_vertex_vertex_angles(FILE* fp, atom atm_pl[], atom atm_npl[], double vert_vert_angles[][4]){
		int p=0;

		for (int p=0; p<9; p++){
			long int l, m, n;
			if (vert_vert_angles[p][0]==0){
				l= atm_pl[0].get_serial_no();
			}else if (vert_vert_angles[p][0]==1){
				l= atm_pl[1].get_serial_no();
			}else if (vert_vert_angles[p][0]==2){
				l= atm_pl[2].get_serial_no();
			}

			if (vert_vert_angles[p][1]==0){
				m= atm_pl[0].get_serial_no();
			}else if (vert_vert_angles[p][1]==1){
				m= atm_pl[1].get_serial_no();
			}else if (vert_vert_angles[p][1]==2){
				m= atm_pl[2].get_serial_no();
			}else if (vert_vert_angles[p][1]==3){
				m= atm_npl[0].get_serial_no();
			}

			if (vert_vert_angles[p][2]==1){
				n= atm_pl[1].get_serial_no();
			}else if (vert_vert_angles[p][2]==2){
				n= atm_pl[2].get_serial_no();
			}else if (vert_vert_angles[p][2]==3){
				n= atm_npl[0].get_serial_no();
			}

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
		}

		for (int p=9; p<18; p++){
			long int l, m, n;
			if (vert_vert_angles[p][0]==0){
				l= atm_pl[0].get_serial_no();
			}else if (vert_vert_angles[p][0]==1){
				l= atm_pl[1].get_serial_no();
			}else if (vert_vert_angles[p][0]==2){
				l= atm_pl[2].get_serial_no();
			}

			if (vert_vert_angles[p][1]==0){
				m= atm_pl[0].get_serial_no();
			}else if (vert_vert_angles[p][1]==1){
				m= atm_pl[1].get_serial_no();
			}else if (vert_vert_angles[p][1]==2){
				m= atm_pl[2].get_serial_no();
			}else if (vert_vert_angles[p][1]==3){
				m= atm_npl[1].get_serial_no();
			}

			if (vert_vert_angles[p][2]==1){
				n= atm_pl[1].get_serial_no();
			}else if (vert_vert_angles[p][2]==2){
				n= atm_pl[2].get_serial_no();
			}else if (vert_vert_angles[p][2]==3){
				n= atm_npl[1].get_serial_no();
			}

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
		}
	}
