#include "tetragonal_pyramidal.h"

void tetrapyrMetalicStructure::print(atom mt, atom* non_mt, int no_of_atoms){
	mt.fprint(stdout);

	for (int i=0; i<no_of_atoms; i++){
		non_mt[i].fprint(stdout);
	}
}

void tetrapyrMetalicStructure::print(Point3D* vert, int no_of_atoms){
	for (int i=0; i<no_of_atoms; i++){
		cout<<"\nVertex Pl "<<i<<": "<<vert[i].X()<<"  "<<vert[i].Y()<<"  "<<vert[i].Z();
	}
}

void tetrapyrMetalicStructure::print(Point3D vert, int no_of_atoms){
	cout<<"\nVertex Npl: "<<vert.X()<<"  "<<vert.Y()<<"  "<<vert.Z();
}

tetrapyrMetalicStructure::tetrapyrMetalicStructure(){
	_no_of_binding_atoms=5;
	assert(_no_of_binding_atoms==5 && "For Tetragonal Pyramidal No. of Vertices should be 5!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=8;
	assert(_no_of_sides==8 && "For Tetragonal Pyramidal No. of Sides should be 8!!");

	_no_of_planar_atoms=4;
	assert(_no_of_planar_atoms==4 && "For Tetragonal Pyramidal No. of Planar Atoms should be 4!!");

	_no_of_non_planar_atom=1;
	assert(_no_of_non_planar_atom==1 && "For Tetragonal Pyramidal No. of Non Planar Atoms should be 1!!");

	_no_of_vertex_vertex_angles=16;
	assert(_no_of_vertex_vertex_angles==16 && "For Tetragonal Pyramidal No. of Angles of Two Vertices with Another Vertex should be 16!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Tetragonal Pyramidal No. of Features to Compare should be 2!!");

	_npl_vertex_index=0;
}

tetrapyrMetalicStructure::tetrapyrMetalicStructure(tetrahedron left, atom left_atoms[], tetrahedron right, atom right_atoms[], atom metal_atom, vector<long int> platom, vector<long int> nonplatom){
	_no_of_binding_atoms=5;
	assert(_no_of_binding_atoms==5 && "For Tetragonal Pyramidal No. of Vertices should be 5!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=8;
	assert(_no_of_sides==8 && "For Tetragonal Pyramidal No. of Sides should be 8!!");

	_no_of_planar_atoms=4;
	assert(_no_of_planar_atoms==4 && "For Tetragonal Pyramidal No. of Planar Atoms should be 4!!");

	_no_of_non_planar_atom=1;
	assert(_no_of_non_planar_atom==1 && "For Tetragonal Pyramidal No. of Non Planar Atoms should be 1!!");

	_no_of_vertex_vertex_angles=16;
	assert(_no_of_vertex_vertex_angles==16 && "For Tetragonal Pyramidal No. of Angles of Two Vertices with Another Vertex should be 16!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Tetragonal Pyramidal No. of Features to Compare should be 2!!");

	_metal_atom=metal_atom;

	_planar_nos=platom;
	_non_planar_no.push_back(nonplatom.at(0));

	_left_tetra=left;
	for (int i=0; i<4; i++){
		_left_tetra_atoms[i]=left_atoms[i];
	}

	_right_tetra=right;
	for (int i=0; i<4; i++){
		_right_tetra_atoms[i]=right_atoms[i];
	}
}

	tetrapyrMetalicStructure::tetrapyrMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom, vector<long int> platom, vector<long int> nonplatom){
		_no_of_binding_atoms=5;
		assert(_no_of_binding_atoms==5 && "For Tetragonal Pyramidal No. of Vertices should be 5!!");

		_no_of_coords=3;
		assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

		_no_of_sides=8;
		assert(_no_of_sides==8 && "For Tetragonal Pyramidal No. of Sides should be 8!!");

		_no_of_planar_atoms=4;
		assert(_no_of_planar_atoms==4 && "For Tetragonal Pyramidal No. of Planar Atoms should be 4!!");

		_no_of_non_planar_atom=1;
		assert(_no_of_non_planar_atom==1 && "For Tetragonal Pyramidal No. of Non Planar Atoms should be 1!!");

		_no_of_vertex_vertex_angles=16;
		assert(_no_of_vertex_vertex_angles==16 && "For Tetragonal Pyramidal No. of Angles of Two Vertices with Another Vertex should be 16!!");

		_no_of_features_to_compare=2;
		assert(_no_of_features_to_compare==2 && "For Tetragonal Pyramidal No. of Features to Compare should be 2!!");

		_metal_atom=metal_atom;

		atom all_tetrapyr_atoms[_no_of_binding_atoms];
		for (int i=0; i<new_atm.size(); i++){
			all_tetrapyr_atoms[i]=new_atm.at(i).get_atom_direct_contact();
		}

		_planar_nos=platom;
		_non_planar_no.push_back(nonplatom.at(0));

		int p=0;
		for (auto pl: _planar_nos){
			cout<<"\nPL: "<<pl;
			for (int i=0; i<_no_of_binding_atoms; i++){
				if (pl==all_tetrapyr_atoms[i].get_serial_no()){
					_tetrapyr_atoms_pl[p]=all_tetrapyr_atoms[i];
					_tetrapyr_vertices_pl[p]=Point3D(_tetrapyr_atoms_pl[p].X(), _tetrapyr_atoms_pl[p].Y(), _tetrapyr_atoms_pl[p].Z());
					cout<<"\nSr. No.: "<<_tetrapyr_atoms_pl[p].get_serial_no()<<"; Tetra Pyr Pl: "<<_tetrapyr_atoms_pl[p].X()<<", "<<_tetrapyr_atoms_pl[p].Y()<<", "<<_tetrapyr_atoms_pl[p].Z();
					p++;
				}
			}
		}

		for (int i=0; i<_no_of_binding_atoms; i++){
			if (_non_planar_no.at(0)==all_tetrapyr_atoms[i].get_serial_no()){
				_tetrapyr_atom_npl=all_tetrapyr_atoms[i];
				_tetrapyr_vertex_npl=Point3D(_tetrapyr_atom_npl.X(), _tetrapyr_atom_npl.Y(), _tetrapyr_atom_npl.Z());
				cout<<"\nSr. No.: "<<_tetrapyr_atom_npl.get_serial_no()<<"; Tetra Pyr Npl: "<<_tetrapyr_atom_npl.X()<<", "<<_tetrapyr_atom_npl.Y()<<", "<<_tetrapyr_atom_npl.Z();
				break;
			}
		}

		_tetrapyr=tetragonal_pyramidal(_tetrapyr_vertices_pl, _tetrapyr_vertex_npl);

		print(_metal_atom, _tetrapyr_atoms_pl, _no_of_planar_atoms);
		print(_tetrapyr_atom_npl, _no_of_non_planar_atom);

		_npl_vertex_index=0;
	}

	tetrapyrMetalicStructure::tetrapyrMetalicStructure(tetraMetalicStructure left_transformed, tetraMetalicStructure right_transformed, atom metal_atom){
		_no_of_binding_atoms=5;
		assert(_no_of_binding_atoms==5 && "For Tetragonal Pyramidal No. of Vertices should be 5!!");

		_no_of_coords=3;
		assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

		_no_of_sides=8;
		assert(_no_of_sides==8 && "For Tetragonal Pyramidal No. of Sides should be 8!!");

		_no_of_planar_atoms=4;
		assert(_no_of_planar_atoms==4 && "For Tetragonal Pyramidal No. of Planar Atoms should be 4!!");

		_no_of_non_planar_atom=1;
		assert(_no_of_non_planar_atom==1 && "For Tetragonal Pyramidal No. of Non Planar Atoms should be 1!!");

		_no_of_vertex_vertex_angles=16;
		assert(_no_of_vertex_vertex_angles==16 && "For Tetragonal Pyramidal No. of Angles of Two Vertices with Another Vertex should be 16!!");

		_no_of_features_to_compare=2;
		assert(_no_of_features_to_compare==2 && "For Tetragonal Pyramidal No. of Features to Compare should be 2!!");

		_metal_atom=metal_atom;

		_transformed_tetra_left=left_transformed;
		_left_tetra=left_transformed.get_tetrahedron();

		for (int i=0; i<4; i++){
			_left_tetra_atoms[i]=left_transformed.get_tetra_atoms()[i];
		}

		for (int i=0; i<4; i++){
			_left_vertex_order[i]=left_transformed.get_tetrahedron().get_vertex_order()[i];
		}

		_transformed_tetra_right=right_transformed;
		_right_tetra=right_transformed.get_tetrahedron();

		for (int i=0; i<4; i++){
			_right_tetra_atoms[i]=right_transformed.get_tetra_atoms()[i];
		}

		for (int i=0; i<4; i++){
			_right_vertex_order[i]=right_transformed.get_tetrahedron().get_vertex_order()[i];
		}
	}

	tetrapyrMetalicStructure tetrapyrMetalicStructure::apply_transformations(){
		cout<<"\nApp Trans Start";
		cout<<"\n Original Vertices Planar";
		print(_tetrapyr_vertices_pl, _no_of_planar_atoms);

		cout<<"\n Original Vertices Non-Planar";
		print(_tetrapyr_vertex_npl, _no_of_non_planar_atom);

		for (int i=0; i<4; i++){
			_pl_vertex_order[i]=_tetrapyr.get_pl_vertex_order()[i];
		}

		_npl_vertex_index=0;

		find_largest_side_and_update_vertex_order();
		tetrapyrMetalicStructure tms=make_transformation();
		cout<<"\nVERTEX ORDER Planar After Orientation: "<<_pl_vertex_order[0]<<"  "<<_pl_vertex_order[1]<<"  "<<_pl_vertex_order[2]<<"  "<<_pl_vertex_order[3];

		cout<<"\nApp Trans Ends";
		return tms;
	}

	atom tetrapyrMetalicStructure::get_metal(){
		return _metal_atom;
	}

	tetragonal_pyramidal tetrapyrMetalicStructure::get_tetrapyr(){
		return _tetrapyr;
	}

	Point3D* tetrapyrMetalicStructure::get_tetrapyr_vertices_pl(){
		return _tetrapyr_vertices_pl;
	}

	Point3D tetrapyrMetalicStructure::get_tetrapyr_vertex_npl(){
		return _tetrapyr_vertex_npl;
	}

	atom* tetrapyrMetalicStructure::get_tetrapyr_atoms_pl(){
		return _tetrapyr_atoms_pl;
	}

	atom tetrapyrMetalicStructure::get_tetrapyr_atom_npl(){
		return _tetrapyr_atom_npl;
	}

	tetrahedron tetrapyrMetalicStructure::get_left_tetra(){
		return _left_tetra;
	}

	atom* tetrapyrMetalicStructure::get_left_tetra_atoms(){
		return _left_tetra_atoms;
	}

	atom* tetrapyrMetalicStructure::get_right_tetra_atoms(){
		return _right_tetra_atoms;
	}

	tetrahedron tetrapyrMetalicStructure::get_right_tetra(){
		return _right_tetra;
	}

	int tetrapyrMetalicStructure::get_no_of_binding_atoms(){
		return _no_of_binding_atoms;
	}

	double* tetrapyrMetalicStructure::get_side_deviations(){
		return _devt_dist_vert_vert;
	}

	double* tetrapyrMetalicStructure::get_angle_deviations_vertex_vertex(){
		return _devt_angle_vertex_vertex;
	}

	double* tetrapyrMetalicStructure::get_mse(){
		return _mse;
	}

	string tetrapyrMetalicStructure::get_structure_size(){
		return _structure_size;
	}

	void tetrapyrMetalicStructure::comparison(tetrapyrMetalicStructure & comp, string out_comparable_name){
		double dist_stand_planar[4][3];
		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				dist_stand_planar[i][j]=_tetrapyr.get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		double dist_stand_non_planar[4][3];
		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				dist_stand_non_planar[i][j]=_tetrapyr.get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		Point3D* comp_vertices_pl=comp.get_tetrapyr_vertices_pl();
		Point3D comp_vertices_npl=comp.get_tetrapyr_vertex_npl();
		measure_direction_cosines();
		double proj_pl_pl[_no_of_planar_atoms];
		double proj_npl_pl[4];

		int k=0;
		for (int i=0; i<(_no_of_planar_atoms-1); i++){
			for (int j=i+1; j<_no_of_planar_atoms; j++){
				double prx=	_direction_cosines_pl_pl[k][0]*(comp_vertices_pl[j].X()-comp_vertices_pl[i].X());
				double pry= _direction_cosines_pl_pl[k][1]*(comp_vertices_pl[j].Y()-comp_vertices_pl[i].Y());
				double prz=	_direction_cosines_pl_pl[k][2]*(comp_vertices_pl[j].Z()-comp_vertices_pl[i].Z());

				proj_pl_pl[k]=(prx+pry+prz);

				k++;
			}
		}

		k=0;
		for (int j=0; j<_no_of_planar_atoms; j++){
			double prx=	_direction_cosines_npl_pl[k][0]*(comp_vertices_pl[j].X()-comp_vertices_npl.X());
			double pry= _direction_cosines_npl_pl[k][1]*(comp_vertices_pl[j].Y()-comp_vertices_npl.Y());
			double prz=	_direction_cosines_npl_pl[k][2]*(comp_vertices_pl[j].Z()-comp_vertices_npl.Z());

			proj_npl_pl[k]=(prx+pry+prz);

			k++;
		}

		double dist_dev_vert_vert_pl_pl[_no_of_planar_atoms];
		double dist_dev_vert_vert_npl_pl[4];

		FILE* outf=fopen(out_comparable_name.c_str(), "w");

		cout<<"\n\nComparison Dist:\n";
		for (int i=0; i<_no_of_planar_atoms; i++){
			double diff=dist_stand_planar[i][2]-proj_pl_pl[i];
			dist_dev_vert_vert_pl_pl[i]=diff;
			cout<<"Side "<<(i+1)<<": "<<dist_dev_vert_vert_pl_pl[i]<<endl;
		}

		for (int i=0; i<4; i++){
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
		mse_dist[1]=(sum2/4);
		double mse_distance=mse_dist[0]+mse_dist[1];
		cout<<"\n MSE Distance: "<<mse_distance;
		fprintf(outf, "\nOriginal MSE Dist: %8.3lf", mse_distance);

		cout<<"\n\nComparison Angle:\n";
		double angle_stand_planar[4][4];
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				angle_stand_planar[i][j]=_tetrapyr.get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		double angle_stand_planar_non_planar[12][4];
		for (int i=0; i<12; i++){
			for (int j=0; j<4; j++){
				angle_stand_planar_non_planar[i][j]=_tetrapyr.get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}


		double angle_comp_planar[4][4];
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				angle_comp_planar[i][j]=comp.get_tetrapyr().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		double angle_comp_planar_non_planar[12][4];
		for (int i=0; i<12; i++){
			for (int j=0; j<4; j++){
				angle_comp_planar_non_planar[i][j]=comp.get_tetrapyr().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		double angle_dev_vert_vert_pl_pl[4];
		double angle_dev_vert_vert_npl_pl[12];

		for (int p=0; p<_no_of_planar_atoms; p++){
			angle_dev_vert_vert_pl_pl[p]=angle_stand_planar[p][3]-angle_comp_planar[p][3];
			cout<<"\nAngle "<<(p+1)<<": "<<angle_dev_vert_vert_pl_pl[p];
		}

		for (int p=0; p<12; p++){
			angle_dev_vert_vert_npl_pl[p]=angle_stand_planar_non_planar[p][3]-angle_comp_planar_non_planar[p][3];
			cout<<"\nAngle "<<(p+1)<<": "<<angle_dev_vert_vert_npl_pl[p];
		}

		sum1=0.0;
		for (int p=0; p<_no_of_planar_atoms; p++){
			sum1+=(angle_dev_vert_vert_pl_pl[p] * angle_dev_vert_vert_pl_pl[p]);
		}

		sum2=0.0;
		for (int p=0; p<12; p++){
			sum2+=(angle_dev_vert_vert_npl_pl[p] * angle_dev_vert_vert_npl_pl[p]);
		}

		double mse_angle[2];
		mse_angle[0]=(sum1/_no_of_planar_atoms);
		mse_angle[1]=(sum2/12);
		double mse_ang=mse_angle[0]+mse_angle[1];
		cout<<"\n MSE Angle: "<<mse_ang;
		fprintf(outf, "\nMSE Angle: %8.3lf", mse_ang);

		fclose(outf);
	}

	void tetrapyrMetalicStructure:: find_largest_side_and_update_vertex_order(){
		double dist[4][3];
		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				dist[i][j]=_tetrapyr.get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		int* _initial_vertex_order=new int[4];

		for (int i=0; i<4; i++){
			_initial_vertex_order[i]=_tetrapyr.get_pl_vertex_order()[i];
		}

		cout<<"\nInitial V. Order: "<<_pl_vertex_order[0]<<", "<<_pl_vertex_order[1]<<", "<<_pl_vertex_order[2]<<", "<<_pl_vertex_order[3];
		if ((dist[0][2]>=dist[1][2]) && (dist[0][2]>=dist[2][2]) && (dist[0][2]>=dist[3][2])){
			_pl_vertex_order[0]=_initial_vertex_order[0];
			_pl_vertex_order[1]=_initial_vertex_order[1];

			if (dist[1][2]>=dist[3][2]){
				_pl_vertex_order[2]=_initial_vertex_order[2];
				_pl_vertex_order[3]=_initial_vertex_order[3];
			}else{
				_pl_vertex_order[0]=_initial_vertex_order[1];
				_pl_vertex_order[1]=_initial_vertex_order[0];
				_pl_vertex_order[2]=_initial_vertex_order[3];
				_pl_vertex_order[3]=_initial_vertex_order[2];
			}
		}else if ((dist[1][2]>=dist[0][2]) && (dist[1][2]>=dist[2][2]) && (dist[1][2]>=dist[3][2])){
			_pl_vertex_order[0]=_initial_vertex_order[1];
			_pl_vertex_order[1]=_initial_vertex_order[2];

			if (dist[2][2]>=dist[0][2]){
				_pl_vertex_order[2]=_initial_vertex_order[3];
				_pl_vertex_order[3]=_initial_vertex_order[0];
			}else{
				_pl_vertex_order[0]=_initial_vertex_order[2];
				_pl_vertex_order[1]=_initial_vertex_order[1];
				_pl_vertex_order[2]=_initial_vertex_order[0];
				_pl_vertex_order[3]=_initial_vertex_order[3];
			}
		}else if ((dist[2][2]>=dist[0][2]) && (dist[2][2]>=dist[1][2]) && (dist[2][2]>=dist[3][2])){
			_pl_vertex_order[0]=_initial_vertex_order[2];
			_pl_vertex_order[1]=_initial_vertex_order[3];

			if (dist[3][2]>=dist[1][2]){
				_pl_vertex_order[2]=_initial_vertex_order[0];
				_pl_vertex_order[3]=_initial_vertex_order[1];
			}else{
				_pl_vertex_order[0]=_initial_vertex_order[3];
				_pl_vertex_order[1]=_initial_vertex_order[2];
				_pl_vertex_order[2]=_initial_vertex_order[1];
				_pl_vertex_order[3]=_initial_vertex_order[0];
			}
		}else if ((dist[3][2]>=dist[0][2]) && (dist[3][2]>=dist[1][2]) && (dist[3][2]>=dist[2][2])){
			_pl_vertex_order[0]=_initial_vertex_order[3];
			_pl_vertex_order[1]=_initial_vertex_order[0];

			if (dist[0][2]>=dist[2][2]){
				_pl_vertex_order[2]=_initial_vertex_order[1];
				_pl_vertex_order[3]=_initial_vertex_order[2];
			}else{
				_pl_vertex_order[0]=_initial_vertex_order[0];
				_pl_vertex_order[1]=_initial_vertex_order[3];
				_pl_vertex_order[2]=_initial_vertex_order[2];
				_pl_vertex_order[3]=_initial_vertex_order[1];
			}
		}

		cout<<"\nLar Size Up V Order: "<<_pl_vertex_order[0]<<", "<<_pl_vertex_order[1]<<", "<<_pl_vertex_order[2]<<", "<<_pl_vertex_order[3];
	}

	tetrapyrMetalicStructure tetrapyrMetalicStructure::make_transformation(){
		Point3D* left_tetra_vertices=new Point3D[4];
		atom left_tetra_atoms[4];
		Point3D* right_tetra_vertices=new Point3D[4];
		atom right_tetra_atoms[4];

		left_tetra_vertices[0]=_tetrapyr_vertices_pl[_pl_vertex_order[0]];
		left_tetra_vertices[1]=_tetrapyr_vertices_pl[_pl_vertex_order[1]];
		left_tetra_vertices[2]=_tetrapyr_vertices_pl[_pl_vertex_order[3]];
		left_tetra_vertices[3]=_tetrapyr_vertex_npl;
		left_tetra_atoms[0]=_tetrapyr_atoms_pl[_pl_vertex_order[0]];
		left_tetra_atoms[1]=_tetrapyr_atoms_pl[_pl_vertex_order[1]];
		left_tetra_atoms[2]=_tetrapyr_atoms_pl[_pl_vertex_order[3]];
		left_tetra_atoms[3]=_tetrapyr_atom_npl;

		tetraMetalicStructure leftTstruct=tetraMetalicStructure(left_tetra_atoms, left_tetra_vertices, _metal_atom, 0, 1, 2, 3);
		delete [] left_tetra_vertices;
		tetraMetalicStructure  left_transformed=leftTstruct.transform_tetrahedron();

		right_tetra_vertices[0]=_tetrapyr_vertices_pl[_pl_vertex_order[3]];
		right_tetra_vertices[1]=_tetrapyr_vertices_pl[_pl_vertex_order[1]];
		right_tetra_vertices[2]=_tetrapyr_vertices_pl[_pl_vertex_order[2]];
		right_tetra_vertices[3]=_tetrapyr_vertex_npl;
		right_tetra_atoms[0]=_tetrapyr_atoms_pl[_pl_vertex_order[3]];
		right_tetra_atoms[1]=_tetrapyr_atoms_pl[_pl_vertex_order[1]];
		right_tetra_atoms[2]=_tetrapyr_atoms_pl[_pl_vertex_order[2]];
		right_tetra_atoms[3]=_tetrapyr_atom_npl;

		tetraMetalicStructure rightTstruct=tetraMetalicStructure(right_tetra_atoms, right_tetra_vertices, _metal_atom, 0, 1, 2, 3);
		delete [] right_tetra_vertices;
		tetraMetalicStructure  right_transformed=rightTstruct.transform_tetrahedron();

		tetrapyrMetalicStructure tps=tetrapyrMetalicStructure(left_transformed, right_transformed, _metal_atom);

		return tps;
	}

	void tetrapyrMetalicStructure::measure_direction_cosines(){
		int k=0;
		double x, y, z, l, m, n, r;
		for (int i=0; i<_no_of_planar_atoms-1; i++){
			for (int j=i+1; j<_no_of_planar_atoms; j++){
				x=_tetrapyr_vertices_pl[j].X() - _tetrapyr_vertices_pl[i].X();
				y=_tetrapyr_vertices_pl[j].Y() - _tetrapyr_vertices_pl[i].Y();
				z=_tetrapyr_vertices_pl[j].Z() - _tetrapyr_vertices_pl[i].Z();
				r=sqrt(x*x+y*y+z*z);

				l=x/r;
				m=y/r;
				n=z/r;

				_direction_cosines_pl_pl[k][0]=l;
				_direction_cosines_pl_pl[k][1]=m;
				_direction_cosines_pl_pl[k][2]=n;
				k++;
			}
		}

		k=0;
		for (int j=0; j<_no_of_planar_atoms; j++){
			x=_tetrapyr_vertices_pl[j].X() - _tetrapyr_vertex_npl.X();
			y=_tetrapyr_vertices_pl[j].Y() - _tetrapyr_vertex_npl.Y();
			z=_tetrapyr_vertices_pl[j].Z() - _tetrapyr_vertex_npl.Z();
			r=sqrt(x*x+y*y+z*z);

			l=x/r;
			m=y/r;
			n=z/r;

			_direction_cosines_npl_pl[k][0]=l;
			_direction_cosines_npl_pl[k][1]=m;
			_direction_cosines_npl_pl[k][2]=n;
			k++;
		}
	}

	void tetrapyrMetalicStructure::write_atoms_pdb_format(string fname){
		FILE* outf=fopen(fname.c_str(), "w");

		for (int i=0; i<_no_of_planar_atoms; i++){
			fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _tetrapyr_atoms_pl[i].get_atom_type().c_str(), _tetrapyr_atoms_pl[i].get_serial_no(), _tetrapyr_atoms_pl[i].get_atom_name().c_str(), _tetrapyr_atoms_pl[i].get_residue_name().c_str(), _tetrapyr_atoms_pl[i].get_chain_name().c_str(), _tetrapyr_atoms_pl[i].get_residue_no(), _tetrapyr_atoms_pl[i].X(), _tetrapyr_atoms_pl[i].Y(), _tetrapyr_atoms_pl[i].Z(), _tetrapyr_atoms_pl[i].get_occupancy(), _tetrapyr_atoms_pl[i].get_temp_factor(), _tetrapyr_atoms_pl[i].get_atom_symbol().c_str());
		}

		fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _tetrapyr_atom_npl.get_atom_type().c_str(), _tetrapyr_atom_npl.get_serial_no(), _tetrapyr_atom_npl.get_atom_name().c_str(), _tetrapyr_atom_npl.get_residue_name().c_str(), _tetrapyr_atom_npl.get_chain_name().c_str(), _tetrapyr_atom_npl.get_residue_no(), _tetrapyr_atom_npl.X(), _tetrapyr_atom_npl.Y(), _tetrapyr_atom_npl.Z(), _tetrapyr_atom_npl.get_occupancy(), _tetrapyr_atom_npl.get_temp_factor(), _tetrapyr_atom_npl.get_atom_symbol().c_str());

		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[0].get_serial_no(), _tetrapyr_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[1].get_serial_no(), _tetrapyr_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[1].get_serial_no(), _tetrapyr_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[2].get_serial_no(), _tetrapyr_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[2].get_serial_no(), _tetrapyr_atoms_pl[3].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[3].get_serial_no(), _tetrapyr_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[3].get_serial_no(), _tetrapyr_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[0].get_serial_no(), _tetrapyr_atoms_pl[3].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atom_npl.get_serial_no(), _tetrapyr_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[0].get_serial_no(), _tetrapyr_atom_npl.get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atom_npl.get_serial_no(), _tetrapyr_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[1].get_serial_no(), _tetrapyr_atom_npl.get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atom_npl.get_serial_no(), _tetrapyr_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[2].get_serial_no(), _tetrapyr_atom_npl.get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atom_npl.get_serial_no(), _tetrapyr_atoms_pl[3].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapyr_atoms_pl[3].get_serial_no(), _tetrapyr_atom_npl.get_serial_no());

		fclose(outf);
	}

	void tetrapyrMetalicStructure::write_atoms_combined_tex_format_left(tetrapyrMetalicStructure & comp, string fname){
		string out_file_name=fname+".tex";
		FILE* outf=fopen(out_file_name.c_str(), "w");
		fprintf(outf, "\\documentclass[a4paper]{article}\n");
		fprintf(outf, "\\usepackage[utf8]{inputenc}\n");
		fprintf(outf, "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}");
		fprintf(outf, "\\usepackage{chemfig}\n");
		fprintf(outf, "\\begin{document}\n");
		fprintf(outf, "\\setbondoffset{3pt}\n");
		string st_left[4], st1, str;
		double bond_angle_standard_left[4];
		double bond_angle_comparable_left[4];
		//atom** _comp_tetrapyr_atoms=comp.get_tetrapyr_atoms();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		tetrahedron left_comp=comp.get_left_tetra();

		Point3D* left_vertices_standard=_left_tetra.get_vertices();
		Point3D* left_vertices_comparable=left_comp.get_vertices();

		int* stand_vertex_order_left=_left_tetra.get_vertex_order();
		int* comp_vertex_order_left=left_comp.get_vertex_order();

		double bond_len_standard_left=_left_tetra.get_circum_radius();
		double bond_len_comparable_left=left_comp.get_circum_radius();

		string adjust_stand_left="";
		string adjust_comp_left="";

		if (bond_len_standard_left>=8.0){
			bond_len_standard_left=3.5;
			adjust_stand_left="adjusted";
		}

		if (bond_len_comparable_left>=8.0){
			bond_len_comparable_left=3.5;
			adjust_comp_left="adjusted";
		}

		str="\\center\\chemname[30pt]{\\chemfig{\\chembelow[4pt]{Circumcentre}{Origin}";
		fprintf(outf, "%s", str.c_str());

		st_left[0]="-";
		st_left[1]="-";
		st_left[2]="-";
		st_left[3]="<";

		//Compute bond angles for Standard left Tetrahedron
		double ang1=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_standard[stand_vertex_order_left[0]], &left_vertices_standard[stand_vertex_order_left[1]]);
		bond_angle_standard_left[0]=180.0+ang1;
		double ang2=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_standard[stand_vertex_order_left[1]], &left_vertices_standard[stand_vertex_order_left[0]]);
		bond_angle_standard_left[1]=-ang2;

		if (left_vertices_standard[stand_vertex_order_left[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&left_vertices_standard[stand_vertex_order_left[1]], &_centre, &left_vertices_standard[stand_vertex_order_left[2]]);
			bond_angle_standard_left[2]=ang-ang2;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&left_vertices_standard[stand_vertex_order_left[0]], &_centre, &left_vertices_standard[stand_vertex_order_left[2]]);
			bond_angle_standard_left[2]=ang-ang1;
		}

		if ((left_vertices_standard[stand_vertex_order_left[3]].X()>=0.0) && (left_vertices_standard[stand_vertex_order_left[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_standard[stand_vertex_order_left[3]], &left_vertices_standard[stand_vertex_order_left[2]]);
			bond_angle_standard_left[3]=90-ang/2.0;
		}else if ((left_vertices_standard[stand_vertex_order_left[3]].X()<0.0) && (left_vertices_standard[stand_vertex_order_left[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_standard[stand_vertex_order_left[3]], &left_vertices_standard[stand_vertex_order_left[2]]);
			bond_angle_standard_left[3]=90+ang/2.0;
		}else if ((left_vertices_standard[stand_vertex_order_left[3]].X()<0.0) && (left_vertices_standard[stand_vertex_order_left[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_standard[stand_vertex_order_left[3]], &left_vertices_standard[stand_vertex_order_left[0]]);
			bond_angle_standard_left[3]=180+ang/2.0;
		}else  if ((left_vertices_standard[stand_vertex_order_left[3]].X()>=0.0) && (left_vertices_standard[stand_vertex_order_left[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_standard[stand_vertex_order_left[1]], &left_vertices_standard[stand_vertex_order_left[3]]);
			bond_angle_standard_left[3]=360-ang/2.0;
		}

		//Compute bond angles for Comparable left Tetrahedron
		double ang3=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_comparable[comp_vertex_order_left[0]], &left_vertices_comparable[comp_vertex_order_left[1]]);
		bond_angle_comparable_left[0]=180.0+ang3;
		double ang4=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_comparable[comp_vertex_order_left[1]], &left_vertices_comparable[comp_vertex_order_left[0]]);
		bond_angle_comparable_left[1]=-ang4;

		if (left_vertices_comparable[comp_vertex_order_left[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&left_vertices_comparable[comp_vertex_order_left[1]], &_centre, &left_vertices_comparable[comp_vertex_order_left[2]]);
			bond_angle_comparable_left[2]=ang-ang4;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&left_vertices_comparable[comp_vertex_order_left[0]], &_centre, &left_vertices_comparable[comp_vertex_order_left[2]]);
			bond_angle_comparable_left[2]=ang-ang3;
		}

		if ((left_vertices_comparable[comp_vertex_order_left[3]].X()>=0.0) && (left_vertices_comparable[comp_vertex_order_left[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_comparable[comp_vertex_order_left[3]], &left_vertices_comparable[comp_vertex_order_left[2]]);
			bond_angle_comparable_left[3]=90-ang/2.0;
		}else if ((left_vertices_comparable[comp_vertex_order_left[3]].X()<0.0) && (left_vertices_comparable[comp_vertex_order_left[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_comparable[comp_vertex_order_left[3]], &left_vertices_comparable[comp_vertex_order_left[2]]);
			bond_angle_comparable_left[3]=90+ang/2.0;
		}else if ((left_vertices_comparable[comp_vertex_order_left[3]].X()<0.0) && (left_vertices_comparable[comp_vertex_order_left[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_comparable[comp_vertex_order_left[3]], &left_vertices_comparable[comp_vertex_order_left[0]]);
			bond_angle_comparable_left[3]=180+ang/2.0;
		}else  if ((left_vertices_comparable[comp_vertex_order_left[3]].X()>=0.0) && (left_vertices_comparable[comp_vertex_order_left[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &left_vertices_comparable[comp_vertex_order_left[1]], &left_vertices_comparable[comp_vertex_order_left[3]]);
			bond_angle_comparable_left[3]=360-ang/2.0;
		}

		//string main_string="";
		for (int p=0; p<4; p++){
			//Write Standard Tetrahedron left Vertices
			str="";

			if (_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()=="HOH"){
				st1=_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()=="O"){
					st1=_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()=="OXT"){
					st1=_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_left[p].c_str(), bond_angle_standard_left[p], bond_len_standard_left, str.c_str());

			//Write Comparable Tetrahedron left Vertices
			str="";
			atom* comp_left_tetra_atoms=comp.get_left_tetra_atoms();
			if (comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_name()=="HOH"){
				st1=comp_left_tetra_atoms[comp_vertex_order_left[p]].get_chain_name()+":"+comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_name()+"_{"+to_string(comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{blue}H_2O}{\\color{blue}"+st1+"}";
			}else{
				if (comp_left_tetra_atoms[comp_vertex_order_left[p]].get_atom_name()=="O"){
					st1=comp_left_tetra_atoms[comp_vertex_order_left[p]].get_chain_name()+":"+comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_name()+"_{"+to_string(comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_left_tetra_atoms[comp_vertex_order_left[p]].get_atom_symbol()+"1}{\\color{blue}"+st1+"}";
				}else if (comp_left_tetra_atoms[comp_vertex_order_left[p]].get_atom_name()=="OXT"){
					st1=comp_left_tetra_atoms[comp_vertex_order_left[p]].get_chain_name()+":"+comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_name()+"_{"+to_string(comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_left_tetra_atoms[comp_vertex_order_left[p]].get_atom_symbol()+"2}{\\color{blue}"+st1+"}";
				}else{
					st1=comp_left_tetra_atoms[comp_vertex_order_left[p]].get_chain_name()+":"+comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_name()+"_{"+to_string(comp_left_tetra_atoms[comp_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_left_tetra_atoms[comp_vertex_order_left[p]].get_atom_name()+"}{\\color{blue}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,blue]%s)", st_left[p].c_str(), bond_angle_comparable_left[p], bond_len_comparable_left, str.c_str());
		}

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms left}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Tetragonal Pyramidal Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Tetragonal Pyramidal Site to Any Atom\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); If Bond Length Exceeds Page Width, It Is Set Equal to the Page Width\\end{flushleft}";

		fprintf(outf, "%s", str.c_str());

		if (adjust_stand_left!=""){
			fprintf(outf, "\n\\vspace{2ex}Here Bond Length for Standard Site is Adjusted to the Page Width");
		}

		if (adjust_comp_left!=""){
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

	void tetrapyrMetalicStructure::write_atoms_combined_tex_format_right(tetrapyrMetalicStructure & comp, string fname){
		string out_file_name=fname+".tex";
		FILE* outf=fopen(out_file_name.c_str(), "w");
		fprintf(outf, "\\documentclass[a4paper]{article}\n");
		fprintf(outf, "\\usepackage[utf8]{inputenc}\n");
		fprintf(outf, "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}");
		fprintf(outf, "\\usepackage{chemfig}\n");
		fprintf(outf, "\\begin{document}\n");
		fprintf(outf, "\\setbondoffset{3pt}\n");
		string st_right[4], st1, str;
		double bond_angle_standard_right[4];
		double bond_angle_comparable_right[4];
		//atom** _comp_tetrapyr_atoms=comp.get_tetrapyr_atoms();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		tetrahedron right_comp=comp.get_right_tetra();

		Point3D* right_vertices_standard=_right_tetra.get_vertices();
		Point3D* right_vertices_comparable=right_comp.get_vertices();

		int* stand_vertex_order_right=_right_tetra.get_vertex_order();
		int* comp_vertex_order_right=right_comp.get_vertex_order();

		double bond_len_standard_right=_right_tetra.get_circum_radius();
		double bond_len_comparable_right=right_comp.get_circum_radius();

		string adjust_stand_right="";
		string adjust_comp_right="";

		if (bond_len_standard_right>=8.0){
			bond_len_standard_right=3.5;
			adjust_stand_right="adjusted";
		}

		if (bond_len_comparable_right>=8.0){
			bond_len_comparable_right=3.5;
			adjust_comp_right="adjusted";
		}

		str="\\center\\chemname[30pt]{\\chemfig{\\chembelow[4pt]{Circumcentre}{Origin}";
		fprintf(outf, "%s", str.c_str());

		st_right[0]="-";
		st_right[1]="-";
		st_right[2]="-";
		st_right[3]="<:";

		//Compute bond angles for Standard right Tetrahedron
		double ang11=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_standard[stand_vertex_order_right[0]], &right_vertices_standard[stand_vertex_order_right[1]]);
		bond_angle_standard_right[0]=180.0+ang11;
		double ang12=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_standard[stand_vertex_order_right[1]], &right_vertices_standard[stand_vertex_order_right[0]]);
		bond_angle_standard_right[1]=-ang12;

		if (right_vertices_standard[stand_vertex_order_right[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&right_vertices_standard[stand_vertex_order_right[1]], &_centre, &right_vertices_standard[stand_vertex_order_right[2]]);
			bond_angle_standard_right[2]=ang-ang12;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&right_vertices_standard[stand_vertex_order_right[0]], &_centre, &right_vertices_standard[stand_vertex_order_right[2]]);
			bond_angle_standard_right[2]=ang-ang11;
		}

		if ((right_vertices_standard[stand_vertex_order_right[3]].X()>=0.0) && (right_vertices_standard[stand_vertex_order_right[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_standard[stand_vertex_order_right[3]], &right_vertices_standard[stand_vertex_order_right[2]]);
			bond_angle_standard_right[3]=90-ang/2.0;
		}else if ((right_vertices_standard[stand_vertex_order_right[3]].X()<0.0) && (right_vertices_standard[stand_vertex_order_right[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_standard[stand_vertex_order_right[3]], &right_vertices_standard[stand_vertex_order_right[2]]);
			bond_angle_standard_right[3]=90+ang/2.0;
		}else if ((right_vertices_standard[stand_vertex_order_right[3]].X()<0.0) && (right_vertices_standard[stand_vertex_order_right[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_standard[stand_vertex_order_right[3]], &right_vertices_standard[stand_vertex_order_right[0]]);
			bond_angle_standard_right[3]=180+ang/2.0;
		}else  if ((right_vertices_standard[stand_vertex_order_right[3]].X()>=0.0) && (right_vertices_standard[stand_vertex_order_right[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_standard[stand_vertex_order_right[1]], &right_vertices_standard[stand_vertex_order_right[3]]);
			bond_angle_standard_right[3]=360-ang/2.0;
		}

		//Compute bond angles for Comparable right Tetrahedron
		double ang13=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_comparable[comp_vertex_order_right[0]], &right_vertices_comparable[comp_vertex_order_right[1]]);
		bond_angle_comparable_right[0]=180.0+ang13;
		double ang14=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_comparable[comp_vertex_order_right[1]], &right_vertices_comparable[comp_vertex_order_right[0]]);
		bond_angle_comparable_right[1]=-ang14;

		if (right_vertices_comparable[comp_vertex_order_right[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&right_vertices_comparable[comp_vertex_order_right[1]], &_centre, &right_vertices_comparable[comp_vertex_order_right[2]]);
			bond_angle_comparable_right[2]=ang-ang14;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&right_vertices_comparable[comp_vertex_order_right[0]], &_centre, &right_vertices_comparable[comp_vertex_order_right[2]]);
			bond_angle_comparable_right[2]=ang-ang13;
		}

		if ((right_vertices_comparable[comp_vertex_order_right[3]].X()>=0.0) && (right_vertices_comparable[comp_vertex_order_right[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_comparable[comp_vertex_order_right[3]], &right_vertices_comparable[comp_vertex_order_right[2]]);
			bond_angle_comparable_right[3]=90-ang/2.0;
		}else if ((right_vertices_comparable[comp_vertex_order_right[3]].X()<0.0) && (right_vertices_comparable[comp_vertex_order_right[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_comparable[comp_vertex_order_right[3]], &right_vertices_comparable[comp_vertex_order_right[2]]);
			bond_angle_comparable_right[3]=90+ang/2.0;
		}else if ((right_vertices_comparable[comp_vertex_order_right[3]].X()<0.0) && (right_vertices_comparable[comp_vertex_order_right[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_comparable[comp_vertex_order_right[3]], &right_vertices_comparable[comp_vertex_order_right[0]]);
			bond_angle_comparable_right[3]=180+ang/2.0;
		}else  if ((right_vertices_comparable[comp_vertex_order_right[3]].X()>=0.0) && (right_vertices_comparable[comp_vertex_order_right[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_comparable[comp_vertex_order_right[1]], &right_vertices_comparable[comp_vertex_order_right[3]]);
			bond_angle_comparable_right[3]=360-ang/2.0;
		}

		for (int p=0; p<4; p++){
			//Write Standard Tetrahedron right Vertices
			str="";

			if (_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()=="HOH"){
				st1=_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()=="O"){
					st1=_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()=="OXT"){
					st1=_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_right[p].c_str(), bond_angle_standard_right[p], bond_len_standard_right, str.c_str());

			//Write Comparable Tetrahedron right Vertices
			str="";
			atom* comp_right_tetra_atoms=comp.get_right_tetra_atoms();
			if (comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_name()=="HOH"){
				st1=comp_right_tetra_atoms[comp_vertex_order_right[p]].get_chain_name()+":"+comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_name()+"_{"+to_string(comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{blue}H_2O}{\\color{blue}"+st1+"}";
			}else{
				if (comp_right_tetra_atoms[comp_vertex_order_right[p]].get_atom_name()=="O"){
					st1=comp_right_tetra_atoms[comp_vertex_order_right[p]].get_chain_name()+":"+comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_name()+"_{"+to_string(comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_right_tetra_atoms[comp_vertex_order_right[p]].get_atom_symbol()+"1}{\\color{blue}"+st1+"}";
				}else if (comp_right_tetra_atoms[comp_vertex_order_right[p]].get_atom_name()=="OXT"){
					st1=comp_right_tetra_atoms[comp_vertex_order_right[p]].get_chain_name()+":"+comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_name()+"_{"+to_string(comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_right_tetra_atoms[comp_vertex_order_right[p]].get_atom_symbol()+"2}{\\color{blue}"+st1+"}";
				}else{
					st1=comp_right_tetra_atoms[comp_vertex_order_right[p]].get_chain_name()+":"+comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_name()+"_{"+to_string(comp_right_tetra_atoms[comp_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+comp_right_tetra_atoms[comp_vertex_order_right[p]].get_atom_name()+"}{\\color{blue}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,blue]%s)", st_right[p].c_str(), bond_angle_comparable_right[p], bond_len_comparable_right, str.c_str());
		}

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms right}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Tetragonal Pyramidal Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Tetragonal Pyramidal Site to Any Atom\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); If Bond Length Exceeds Page Width, It Is Set Equal to the Page Width\\end{flushleft}";

		fprintf(outf, "%s", str.c_str());

		if (adjust_stand_right!=""){
			fprintf(outf, "\n\\vspace{2ex}Here Bond Length for Standard Site is Adjusted to the Page Width");
		}

		if (adjust_comp_right!=""){
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

	void tetrapyrMetalicStructure::compare(tetrapyrMetalicStructure & comp){
		double deviate_dist_left[6];
		double deviate_dist_right[6];

		_left_tetra.deviation_distances(comp.get_left_tetra(), _left_tetra.measure_direction_cosines_of_sides(), deviate_dist_left);
		_right_tetra.deviation_distances(comp.get_right_tetra(), _right_tetra.measure_direction_cosines_of_sides(), deviate_dist_right);

		int k=0;
		for (int p=0; p<6; p++){
			cout<<"\nleft Dist:"<<deviate_dist_left[p]<<", right Dist: "<<deviate_dist_right[p];
			_devt_dist_vert_vert_left[p]=deviate_dist_left[p];
			_devt_dist_vert_vert_right[p]=deviate_dist_right[p];

			if (p!=3){
				_devt_dist_vert_vert[k]=_devt_dist_vert_vert_left[p];
				k++;
			}
		}

		k=5;
		for (int p=0; p<6; p++){
			if ((p!=0) || (p!=2) || (p!=4)){
				_devt_dist_vert_vert[k]=_devt_dist_vert_vert_right[p];
				k++;
			}
		}

		double deviate_angle_left[12];
		double deviate_angle_right[12];

		cout<<"\nDev. Angle Left: ";
		_left_tetra.deviation_angles_vertex_vertex(comp.get_left_tetra(), deviate_angle_left);

		cout<<"\nDev. Angle Right: ";
		_right_tetra.deviation_angles_vertex_vertex(comp.get_right_tetra(), deviate_angle_right);

		k=0;
		for (int p=0; p<12; p++){
			_devt_angle_vertex_vertex_left[p]=deviate_angle_left[p];
			_devt_angle_vertex_vertex_right[p]=deviate_angle_right[p];

			if ((p==0) || (p==1) || (p==2) || (p==4) || (p==7) || (p==9) || (p==10)){
				_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_left[p];
				k++;
			}
			//cout<<"\nVert-Vert Angle IN COMPARE: "<<(p+1)<<": "<<_devt_angle_vertex_vertex_left[p]<<", "<<_devt_angle_vertex_vertex_right[p];
		}

		k=7;
		for (int p=0; p<12; p++){
			if ((p==2) || (p==5) || (p==6) || (p==7) || (p==8) || (p==10) || (p==11)){
				_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_right[p];
				k++;
			}
		}

		_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_left[4]+_devt_angle_vertex_vertex_right[3];
		k++;
		_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_left[6]+_devt_angle_vertex_vertex_right[0];

		compute_mse();

		_dist_devt_from_cc_left=(_left_tetra.get_circum_radius()-comp.get_left_tetra().get_circum_radius());
		_dist_devt_from_cc_right=(_right_tetra.get_circum_radius()-comp.get_right_tetra().get_circum_radius());

		if ((_dist_devt_from_cc_left==0.0) && (_dist_devt_from_cc_right==0.0)){
			_structure_size="same";
		}else if ((_dist_devt_from_cc_left>0.0) && (_dist_devt_from_cc_right>0.0)){
			_structure_size="smaller";
		}else if ((_dist_devt_from_cc_left<0.0) && (_dist_devt_from_cc_right<0.0)){
			_structure_size="larger";
		}else{
			_structure_size="irregular";
		}

		cout<<"\nCompare End: "<<_structure_size<<", left CCRad Stand: "<<_left_tetra.get_circum_radius()<<", left CCRad Comp: "<<comp.get_left_tetra().get_circum_radius()<<", Diff: "<<_dist_devt_from_cc_left;
		cout<<"\nCompare End: "<<_structure_size<<", right CCRad Stand: "<<_right_tetra.get_circum_radius()<<", right CCRad Comp: "<<comp.get_right_tetra().get_circum_radius()<<", Diff: "<<_dist_devt_from_cc_right;
	}

	void tetrapyrMetalicStructure::compute_mse(){
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

	void tetrapyrMetalicStructure::gen_report_metal_binding_sites(tetrapyrMetalicStructure & ori_stand, tetrapyrMetalicStructure & ori_comp, tetrapyrMetalicStructure & trans_comp, string mt, string standard_site_type, string comparable_site_type, FILE* fp){
		fprintf(fp, "HEADER     PROGRAM NAME: StructureDeviation  VERSION: 1.1");
		string str="TITLE      Deviation Report for Computation of Structural Deviations of Two Tetragonal Pyramidal "+mt+" ION Binding Sites";
		fprintf(fp, "\n%s", str.c_str());

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Standard Tetragonal Pyramid:                   ");
		fprintf(fp, "\nREMARK     1  Molecule Type of Standard Tetragonal Pyramid: %s", standard_site_type.c_str());

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Distances of Atoms-Atoms:                         ");

		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_stand_pl_pl[i][j]=ori_stand.get_tetrapyr().get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_stand.get_tetrapyr_atoms_pl(), _vertex_vertex_distances_stand_pl_pl, 4);

		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_stand_npl_pl[i][j]=ori_stand.get_tetrapyr().get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_stand.get_tetrapyr_atoms_pl(), ori_stand.get_tetrapyr_atom_npl(), _vertex_vertex_distances_stand_npl_pl);

		fprintf(fp, "\nREMARK     1  Angles made by two atoms with another atom of Tetragonal Pyramid:");
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_stand_pl_pl[i][j]=ori_stand.get_tetrapyr().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_stand.get_tetrapyr_atoms_pl(), _vertex_vertex_angles_stand_pl_pl, 4);

		for (int i=0; i<12; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_stand_npl_pl[i][j]=ori_stand.get_tetrapyr().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_stand.get_tetrapyr_atoms_pl(), ori_stand.get_tetrapyr_atom_npl(), _vertex_vertex_angles_stand_npl_pl);

		fprintf(fp, "\nREMARK     2");
		fprintf(fp, "\nREMARK     2  Comparable Tetrahedrons:");
		fprintf(fp, "\nREMARK     1  Molecule Type of Comparable Tetragonal Pyramid: %s", comparable_site_type.c_str());

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

		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_comp_pl_pl[i][j]=ori_comp.get_tetrapyr().get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_comp.get_tetrapyr_atoms_pl(), _vertex_vertex_distances_comp_pl_pl, 4);

		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_comp_npl_pl[i][j]=ori_comp.get_tetrapyr().get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_comp.get_tetrapyr_atoms_pl(), ori_comp.get_tetrapyr_atom_npl(), _vertex_vertex_distances_comp_npl_pl);

		fprintf(fp, "\nREMARK     2  Angles made by two atoms with another atom of Tetragonal Pyramid:");
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_comp_pl_pl[i][j]=ori_comp.get_tetrapyr().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_comp.get_tetrapyr_atoms_pl(), _vertex_vertex_angles_comp_pl_pl, 4);

		for (int i=0; i<12; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_comp_npl_pl[i][j]=ori_comp.get_tetrapyr().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_comp.get_tetrapyr_atoms_pl(), ori_comp.get_tetrapyr_atom_npl(), _vertex_vertex_angles_comp_npl_pl);

		fprintf(fp, "\nREMARK     3");
		fprintf(fp, "\nREMARK     3  Deviation Results of Comparable Tetragonal Pyramid from Standard Tetragonal Pyramid:");

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

	void tetrapyrMetalicStructure::gen_report_vertex_vertex_distances(FILE* fp, atom atm[], double vert_vert_distances[][3], int no_of_sides){
		for (int i=0; i<no_of_sides; i++){
			fprintf(fp, "\nDISTRR %ld %ld %8.3lf", atm[(int) vert_vert_distances[i][0]].get_serial_no(), atm[(int) vert_vert_distances[i][1]].get_serial_no(), vert_vert_distances[i][2]);
		}
	}

	void tetrapyrMetalicStructure::gen_report_vertex_vertex_distances(FILE* fp, atom atm_pl[], atom atm_npl, double vert_vert_distances[][3]){
		int p=0;
		long int l, m;

		l=atm_npl.get_serial_no();
		for (int j=0; j<_no_of_planar_atoms; j++){
			m=atm_pl[j].get_serial_no();

			fprintf(fp, "\nDISTRR %ld %ld %8.3lf", l, m, vert_vert_distances[p][2]);
			p++;
		}
	}

	void tetrapyrMetalicStructure::gen_report_vertex_vertex_angles(FILE* fp, atom atm[], double vert_vert_angles[][4], int no_of_angles){
		for (int i=0; i<no_of_angles; i++){
			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", atm[(int) vert_vert_angles[i][0]].get_serial_no(), atm[(int) vert_vert_angles[i][1]].get_serial_no(), atm[(int) vert_vert_angles[i][2]].get_serial_no(), vert_vert_angles[i][3]);
		}
	}

	void tetrapyrMetalicStructure::gen_report_vertex_vertex_angles(FILE* fp, atom atm_pl[], atom atm_npl, double vert_vert_angles[][4]){
		int p=0, i;
		long int l, m, n;
		for (i=0; i<(_no_of_planar_atoms-1); i++){
			l=atm_pl[i].get_serial_no();
			m=atm_npl.get_serial_no();
			n=atm_pl[i+1].get_serial_no();

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
			p++;

			l=atm_npl.get_serial_no();
			m=atm_pl[i].get_serial_no();
			n=atm_pl[i+1].get_serial_no();

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
			p++;

			l=atm_npl.get_serial_no();
			m=atm_pl[i+1].get_serial_no();
			n=atm_pl[i].get_serial_no();

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
			p++;
		}

		l=atm_pl[i].get_serial_no();
		m=atm_npl.get_serial_no();
		n=atm_pl[0].get_serial_no();

		fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);

		p++;
		l=atm_npl.get_serial_no();
		m=atm_pl[i].get_serial_no();
		n=atm_pl[0].get_serial_no();

		fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);

		p++;
		l=atm_npl.get_serial_no();
		m=atm_pl[0].get_serial_no();
		n=atm_pl[i].get_serial_no();

		fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
	}
