#include "octahedral.h"

void octaMetalicStructure::print(atom mt, atom* non_mt, int no_of_atoms){
	mt.fprint(stdout);

	for (int i=0; i<no_of_atoms; i++){
		non_mt[i].fprint(stdout);
	}
}

void octaMetalicStructure::print(Point3D* vert, int no_of_atoms){
	for (int i=0; i<no_of_atoms; i++){
		cout<<"\nVertex "<<i<<": "<<vert[i].X()<<"  "<<vert[i].Y()<<"  "<<vert[i].Z();
	}
}

void octaMetalicStructure::print(Point3D vert, int no_of_atoms){
	cout<<"\nVertex : "<<vert.X()<<"  "<<vert.Y()<<"  "<<vert.Z();
}

octaMetalicStructure::octaMetalicStructure(){
	_no_of_binding_atoms=6;
	assert(_no_of_binding_atoms==6 && "For Octahedral No. of Vertices should be 6!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=12;
	assert(_no_of_sides==12 && "For Octahedral No. of Sides should be 12!!");

	_no_of_planar_atoms=4;
	assert(_no_of_planar_atoms==4 && "For Octahedral No. of Planar Atoms should be 4!!");

	_no_of_non_planar_atoms=2;
	assert(_no_of_non_planar_atoms==2 && "For Octahedral No. of Non Planar Atoms should be 2!!");

	_no_of_vertex_vertex_angles=28;
	assert(_no_of_vertex_vertex_angles==28 && "For Octahedral No. of Angles of Two Vertices with Another Vertex should be 28!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Octahedral No. of Features to Compare should be 2!!");
}

octaMetalicStructure::octaMetalicStructure(tetrahedron upper_left, atom upper_left_atoms[], tetrahedron upper_right, atom upper_right_atoms[], tetrahedron lower_left, atom lower_left_atoms[], tetrahedron lower_right, atom lower_right_atoms[], atom metal_atom, vector<long int> platom, vector<long int> nonplatom){
	_no_of_binding_atoms=6;
	assert(_no_of_binding_atoms==6 && "For Octahedral No. of Vertices should be 6!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=12;
	assert(_no_of_sides==12 && "For Octahedral No. of Sides should be 12!!");

	_no_of_planar_atoms=4;
	assert(_no_of_planar_atoms==4 && "For Octahedral No. of Planar Atoms should be 4!!");

	_no_of_non_planar_atoms=2;
	assert(_no_of_non_planar_atoms==2 && "For Octahedral No. of Non Planar Atoms should be 2!!");

	_no_of_vertex_vertex_angles=28;
	assert(_no_of_vertex_vertex_angles==28 && "For Octahedral No. of Angles of Two Vertices with Another Vertex should be 28!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Octahedral No. of Features to Compare should be 2!!");

	_metal_atom=metal_atom;

	_planar_nos=platom;
	_non_planar_nos=nonplatom;

	_upper_left_tetra=upper_left;
	for (int i=0; i<4; i++){
		_upper_left_tetra_atoms[i]=upper_left_atoms[i];
	}

	_upper_right_tetra=upper_right;
	for (int i=0; i<4; i++){
		_upper_right_tetra_atoms[i]=upper_right_atoms[i];
	}

	_lower_left_tetra=lower_left;
	for (int i=0; i<4; i++){
		_lower_left_tetra_atoms[i]=lower_left_atoms[i];
	}

	_lower_right_tetra=lower_right;
	for (int i=0; i<4; i++){
		_lower_right_tetra_atoms[i]=lower_right_atoms[i];
	}
}

	octaMetalicStructure::octaMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom, vector<long int> platom, vector<long int> nonplatom){
		_no_of_binding_atoms=6;
		assert(_no_of_binding_atoms==6 && "For Octahedral No. of Vertices should be 6!!");

		_no_of_coords=3;
		assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

		_no_of_sides=12;
		assert(_no_of_sides==12 && "For Octahedral No. of Sides should be 12!!");

		_no_of_planar_atoms=4;
		assert(_no_of_planar_atoms==4 && "For Octahedral No. of Planar Atoms should be 4!!");

		_no_of_non_planar_atoms=2;
		assert(_no_of_non_planar_atoms==2 && "For Octahedral No. of Non Planar Atoms should be 2!!");

		_no_of_vertex_vertex_angles=28;
		assert(_no_of_vertex_vertex_angles==28 && "For Octahedral No. of Angles of Two Vertices with Another Vertex should be 28!!");

		_no_of_features_to_compare=2;
		assert(_no_of_features_to_compare==2 && "For Octahedral No. of Features to Compare should be 2!!");

		_metal_atom=metal_atom;

		atom all_octa_atoms[_no_of_binding_atoms];
		for (int i=0; i<new_atm.size(); i++){
			all_octa_atoms[i]=new_atm.at(i).get_atom_direct_contact();
		}

		_planar_nos=platom;
		_non_planar_nos=nonplatom;

		int p=0;
		for (auto pl: _planar_nos){
			cout<<"\nPL: "<<pl;
			for (int i=0; i<_no_of_binding_atoms; i++){
				if (pl==all_octa_atoms[i].get_serial_no()){
					_octa_atoms_pl[p]=all_octa_atoms[i];
					_octa_vertices_pl[p]=Point3D(_octa_atoms_pl[p].X(), _octa_atoms_pl[p].Y(), _octa_atoms_pl[p].Z());
					cout<<"\nSr. No.: "<<_octa_atoms_pl[p].get_serial_no()<<"; Octa Pl: "<<_octa_atoms_pl[p].X()<<", "<<_octa_atoms_pl[p].Y()<<", "<<_octa_atoms_pl[p].Z();
					p++;
				}
			}
		}

		p=0;
		for (auto npl: _non_planar_nos){
			cout<<"\nNPL: "<<npl;
			for (int i=0; i<_no_of_binding_atoms; i++){
				if (npl==all_octa_atoms[i].get_serial_no()){
					_octa_atoms_npl[p]=all_octa_atoms[i];
					_octa_vertices_npl[p]=Point3D(_octa_atoms_npl[p].X(), _octa_atoms_npl[p].Y(), _octa_atoms_npl[p].Z());
					cout<<"\nSr. No.: "<<_octa_atoms_npl[p].get_serial_no()<<"; Octa Npl: "<<_octa_atoms_npl[p].X()<<", "<<_octa_atoms_npl[p].Y()<<", "<<_octa_atoms_npl[p].Z();
					p++;
				}
			}
		}

		_octa=octahedral(_octa_vertices_pl, _octa_vertices_npl);

		print(_metal_atom, _octa_atoms_pl, _no_of_planar_atoms);
		print(_octa_atoms_npl, _no_of_non_planar_atoms);

		for (int i=0; i<4; i++){
			_pl_vertex_order[i]=_octa.get_pl_vertex_order()[i];
		}

		for (int i=0; i<2; i++){
			_npl_vertex_order[i]=_octa.get_npl_vertex_order()[i];
		}
	}

	octaMetalicStructure::octaMetalicStructure(tetraMetalicStructure upper_left_transformed, tetraMetalicStructure upper_right_transformed, tetraMetalicStructure lower_left_transformed, tetraMetalicStructure lower_right_transformed, atom metal_atom){
		_no_of_binding_atoms=6;
		assert(_no_of_binding_atoms==6 && "For Octahedral No. of Vertices should be 6!!");

		_no_of_coords=3;
		assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

		_no_of_sides=12;
		assert(_no_of_sides==12 && "For Octahedral No. of Sides should be 12!!");

		_no_of_planar_atoms=4;
		assert(_no_of_planar_atoms==4 && "For Octahedral No. of Planar Atoms should be 4!!");

		_no_of_non_planar_atoms=2;
		assert(_no_of_non_planar_atoms==2 && "For Octahedral No. of Non Planar Atoms should be 2!!");

		_no_of_vertex_vertex_angles=28;
		assert(_no_of_vertex_vertex_angles==28 && "For Octahedral No. of Angles of Two Vertices with Another Vertex should be 28!!");

		_no_of_features_to_compare=2;
		assert(_no_of_features_to_compare==2 && "For Octahedral No. of Features to Compare should be 2!!");

		_metal_atom=metal_atom;

		_transformed_tetra_upper_left=upper_left_transformed;
		_upper_left_tetra=upper_left_transformed.get_tetrahedron();

		for (int i=0; i<4; i++){
			_upper_left_tetra_atoms[i]=upper_left_transformed.get_tetra_atoms()[i];
		}

		for (int i=0; i<4; i++){
			_upper_left_vertex_order[i]=upper_left_transformed.get_tetrahedron().get_vertex_order()[i];
		}

		_transformed_tetra_upper_right=upper_right_transformed;
		_upper_right_tetra=upper_right_transformed.get_tetrahedron();

		for (int i=0; i<4; i++){
			_upper_right_tetra_atoms[i]=upper_right_transformed.get_tetra_atoms()[i];
		}

		for (int i=0; i<4; i++){
			_upper_right_vertex_order[i]=upper_right_transformed.get_tetrahedron().get_vertex_order()[i];
		}

		_transformed_tetra_lower_left=lower_left_transformed;
		_lower_left_tetra=lower_left_transformed.get_tetrahedron();

		for (int i=0; i<4; i++){
			_lower_left_tetra_atoms[i]=lower_left_transformed.get_tetra_atoms()[i];
		}

		for (int i=0; i<4; i++){
			_lower_left_vertex_order[i]=lower_left_transformed.get_tetrahedron().get_vertex_order()[i];
		}

		_transformed_tetra_lower_right=lower_right_transformed;
		_lower_right_tetra=lower_right_transformed.get_tetrahedron();

		for (int i=0; i<4; i++){
			_lower_right_tetra_atoms[i]=lower_right_transformed.get_tetra_atoms()[i];
		}

		for (int i=0; i<4; i++){
			_lower_right_vertex_order[i]=lower_right_transformed.get_tetrahedron().get_vertex_order()[i];
		}
	}

	octaMetalicStructure octaMetalicStructure::apply_transformations(){
		cout<<"\n Original Vertices Planar";
		print(_octa_vertices_pl, _no_of_planar_atoms);

		cout<<"\n Original Vertices Non-Planar";
		print(_octa_vertices_npl, _no_of_non_planar_atoms);

		for (int i=0; i<4; i++){
			_pl_vertex_order[i]=_octa.get_pl_vertex_order()[i];
		}

		for (int i=0; i<2; i++){
			_npl_vertex_order[i]=_octa.get_npl_vertex_order()[i];
		}

		find_largest_side_and_update_vertex_order();
		octaMetalicStructure tms=make_transformation();
		cout<<"\nVERTEX ORDER Planar After Orientation: "<<_pl_vertex_order[0]<<"  "<<_pl_vertex_order[1]<<"  "<<_pl_vertex_order[2]<<"  "<<_pl_vertex_order[3];
		cout<<"\nVERTEX ORDER Non Planar After Orientation: "<<_npl_vertex_order[0]<<"  "<<_npl_vertex_order[1];

		return tms;
	}

	atom octaMetalicStructure::get_metal(){
		return _metal_atom;
	}

	octahedral octaMetalicStructure::get_octa(){
		return _octa;
	}

	Point3D* octaMetalicStructure::get_octa_vertices_pl(){
		return _octa_vertices_pl;
	}

	Point3D* octaMetalicStructure::get_octa_vertices_npl(){
		return _octa_vertices_npl;
	}

	atom* octaMetalicStructure::get_octa_atoms_pl(){
		return _octa_atoms_pl;
	}

	atom* octaMetalicStructure::get_octa_atoms_npl(){
		return _octa_atoms_npl;
	}

	int* octaMetalicStructure::get_pl_vertex_order(){
		return _pl_vertex_order;
	}

	int* octaMetalicStructure::get_npl_vertex_order(){
		return _npl_vertex_order;
	}

	tetrahedron octaMetalicStructure::get_upper_left_tetra(){
		return _upper_left_tetra;
	}

	atom* octaMetalicStructure::get_upper_left_tetra_atoms(){
		return _upper_left_tetra_atoms;
	}

	atom* octaMetalicStructure::get_upper_right_tetra_atoms(){
		return _upper_right_tetra_atoms;
	}

	tetrahedron octaMetalicStructure::get_upper_right_tetra(){
		return _upper_right_tetra;
	}

	tetrahedron octaMetalicStructure::get_lower_left_tetra(){
		return _lower_left_tetra;
	}

	atom* octaMetalicStructure::get_lower_left_tetra_atoms(){
		return _lower_left_tetra_atoms;
	}

	atom* octaMetalicStructure::get_lower_right_tetra_atoms(){
		return _lower_right_tetra_atoms;
	}

	tetrahedron octaMetalicStructure::get_lower_right_tetra(){
		return _lower_right_tetra;
	}

	int octaMetalicStructure::get_no_of_binding_atoms(){
		return _no_of_binding_atoms;
	}

	double* octaMetalicStructure::get_side_deviations(){
		return _devt_dist_vert_vert;
	}

	double* octaMetalicStructure::get_angle_deviations_vertex_vertex(){
		return _devt_angle_vertex_vertex;
	}

	double* octaMetalicStructure::get_mse(){
		return _mse;
	}

	string octaMetalicStructure::get_structure_size(){
		return _structure_size;
	}

	void octaMetalicStructure::comparison(octaMetalicStructure & comp, string out_comparable_name){
		double dist_stand_planar[4][3];
		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				dist_stand_planar[i][j]=_octa.get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		double dist_stand_non_planar[4][3];
		for (int i=0; i<8; i++){
			for (int j=0; j<3; j++){
				dist_stand_non_planar[i][j]=_octa.get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		Point3D comp_vertices_pl[4];
		for (int i=0; i<4; i++){
			comp_vertices_pl[i]=comp.get_octa_vertices_pl()[i];
		}

		Point3D comp_vertices_npl[2];
		for (int i=0; i<2; i++){
			comp_vertices_npl[i]=comp.get_octa_vertices_npl()[i];
		}

		measure_direction_cosines();
		double proj_pl_pl[_no_of_planar_atoms];
		double proj_npl_pl[8];

		int k=0;
		for (int i=0; i<(_no_of_planar_atoms-1); i++){
			for (int j=i+1; j<_no_of_planar_atoms; j++){
				if (_octa_vertices_pl[i].dist(&comp_vertices_pl[i])==0 && _octa_vertices_pl[j].dist(&comp_vertices_pl[j])==0){
					proj_pl_pl[k]=comp_vertices_pl[i].dist(&comp_vertices_pl[j]);
					cout<<"\nProj 1st Pl";
				}else{
					cout<<"\nProj 2nd Pl";
					double prx=	_direction_cosines_pl_pl[k][0]*(comp_vertices_pl[j].X()-comp_vertices_pl[i].X());
					double pry= _direction_cosines_pl_pl[k][1]*(comp_vertices_pl[j].Y()-comp_vertices_pl[i].Y());
					double prz=	_direction_cosines_pl_pl[k][2]*(comp_vertices_pl[j].Z()-comp_vertices_pl[i].Z());

					proj_pl_pl[k]=(prx+pry+prz);

					k++;
				}
			}
		}

		k=0;
		for (int i=0; i<8; i++){
			for (int j=0; j<_no_of_planar_atoms; j++){
				if (_octa_vertices_npl[i].dist(&comp_vertices_pl[i])==0 && _octa_vertices_npl[j].dist(&comp_vertices_pl[j])==0){
					proj_npl_pl[k]=comp_vertices_pl[i].dist(&comp_vertices_pl[j]);
					cout<<"\nProj 1st Npl";
				}else{
					cout<<"\nProj 2nd Npl";
					double prx=	_direction_cosines_npl_pl[k][0]*(comp_vertices_pl[j].X()-comp_vertices_npl[i].X());
					double pry= _direction_cosines_npl_pl[k][1]*(comp_vertices_pl[j].Y()-comp_vertices_npl[i].Y());
					double prz=	_direction_cosines_npl_pl[k][2]*(comp_vertices_pl[j].Z()-comp_vertices_npl[i].Z());

					proj_npl_pl[k]=(prx+pry+prz);

					k++;
				}
			}
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

		for (int i=0; i<8; i++){
			double diff=dist_stand_non_planar[i][2]-proj_npl_pl[i];
			dist_dev_vert_vert_npl_pl[i]=diff;
			cout<<"Side "<<(i+1)<<": "<<dist_dev_vert_vert_npl_pl[i]<<endl;
		}

		double sum1=0.0;
		for (int p=0; p<_no_of_planar_atoms; p++){
			sum1+=(dist_dev_vert_vert_pl_pl[p] * dist_dev_vert_vert_pl_pl[p]);
		}

		double sum2=0.0;
		for (int p=0; p<8; p++){
			sum2+=(dist_dev_vert_vert_npl_pl[p] * dist_dev_vert_vert_npl_pl[p]);
		}

		double mse_dist[2];
		mse_dist[0]=(sum1/_no_of_planar_atoms);
		mse_dist[1]=(sum2/8);
		double mse_distance=mse_dist[0]+mse_dist[1];
		cout<<"\n MSE Distance: "<<mse_distance;
		fprintf(outf, "\nOriginal MSE Dist: %8.3lf", mse_distance);

		cout<<"\n\nComparison Angle:\n";
		double angle_stand_planar[4][4];
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				angle_stand_planar[i][j]=_octa.get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		double angle_stand_planar_non_planar[24][4];
		for (int i=0; i<24; i++){
			for (int j=0; j<4; j++){
				angle_stand_planar_non_planar[i][j]=_octa.get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}


		double angle_comp_planar[4][4];
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				angle_comp_planar[i][j]=comp.get_octa().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		double angle_comp_planar_non_planar[24][4];
		for (int i=0; i<24; i++){
			for (int j=0; j<4; j++){
				angle_comp_planar_non_planar[i][j]=comp.get_octa().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		double angle_dev_vert_vert_pl_pl[4];
		double angle_dev_vert_vert_npl_pl[24];

		for (int p=0; p<_no_of_planar_atoms; p++){
			angle_dev_vert_vert_pl_pl[p]=angle_stand_planar[p][3]-angle_comp_planar[p][3];
			cout<<"\nAngle "<<(p+1)<<": "<<angle_dev_vert_vert_pl_pl[p];
		}

		for (int p=0; p<24; p++){
			angle_dev_vert_vert_npl_pl[p]=angle_stand_planar_non_planar[p][3]-angle_comp_planar_non_planar[p][3];
			cout<<"\nAngle "<<(p+1)<<": "<<angle_dev_vert_vert_npl_pl[p];
		}

		sum1=0.0;
		for (int p=0; p<_no_of_planar_atoms; p++){
			sum1+=(angle_dev_vert_vert_pl_pl[p] * angle_dev_vert_vert_pl_pl[p]);
		}

		sum2=0.0;
		for (int p=0; p<24; p++){
			sum2+=(angle_dev_vert_vert_npl_pl[p] * angle_dev_vert_vert_npl_pl[p]);
		}

		double mse_angle[2];
		mse_angle[0]=(sum1/_no_of_planar_atoms);
		mse_angle[1]=(sum2/24);
		double mse_ang=mse_angle[0]+mse_angle[1];
		cout<<"\n MSE Angle: "<<mse_ang;
		fprintf(outf, "\nMSE Angle: %8.3lf", mse_ang);

		fclose(outf);
	}

	void octaMetalicStructure:: find_largest_side_and_update_vertex_order(){
		double dist[4][3];
		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				dist[i][j]=_octa.get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		int _initial_vertex_order[4];

		for (int i=0; i<4; i++){
			_initial_vertex_order[i]=_octa.get_pl_vertex_order()[i];
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

	octaMetalicStructure octaMetalicStructure::make_transformation(){
		Point3D* upper_left_tetra_vertices=new Point3D[4];
		atom upper_left_tetra_atoms[4];
		Point3D* upper_right_tetra_vertices=new Point3D[4];
		atom upper_right_tetra_atoms[4];
		Point3D* lower_left_tetra_vertices=new Point3D[4];
		atom lower_left_tetra_atoms[4];
		Point3D* lower_right_tetra_vertices=new Point3D[4];
		atom lower_right_tetra_atoms[4];

		upper_left_tetra_vertices[0]=_octa_vertices_pl[_pl_vertex_order[0]];
		upper_left_tetra_vertices[1]=_octa_vertices_pl[_pl_vertex_order[1]];
		upper_left_tetra_vertices[2]=_octa_vertices_pl[_pl_vertex_order[3]];
		upper_left_tetra_vertices[3]=_octa_vertices_npl[0];
		upper_left_tetra_atoms[0]=_octa_atoms_pl[_pl_vertex_order[0]];
		upper_left_tetra_atoms[1]=_octa_atoms_pl[_pl_vertex_order[1]];
		upper_left_tetra_atoms[2]=_octa_atoms_pl[_pl_vertex_order[3]];
		upper_left_tetra_atoms[3]=_octa_atoms_npl[0];

		tetraMetalicStructure upper_leftTstruct=tetraMetalicStructure(upper_left_tetra_atoms, upper_left_tetra_vertices, _metal_atom, 0, 1, 2, 3);
		delete [] upper_left_tetra_vertices;

		tetraMetalicStructure  upper_left_transformed=upper_leftTstruct.transform_tetrahedron();

		upper_right_tetra_vertices[0]=_octa_vertices_pl[_pl_vertex_order[3]];
		upper_right_tetra_vertices[1]=_octa_vertices_pl[_pl_vertex_order[1]];
		upper_right_tetra_vertices[2]=_octa_vertices_pl[_pl_vertex_order[2]];
		upper_right_tetra_vertices[3]=_octa_vertices_npl[0];
		upper_right_tetra_atoms[0]=_octa_atoms_pl[_pl_vertex_order[3]];
		upper_right_tetra_atoms[1]=_octa_atoms_pl[_pl_vertex_order[1]];
		upper_right_tetra_atoms[2]=_octa_atoms_pl[_pl_vertex_order[2]];
		upper_right_tetra_atoms[3]=_octa_atoms_npl[0];

		tetraMetalicStructure upper_rightTstruct=tetraMetalicStructure(upper_right_tetra_atoms, upper_right_tetra_vertices, _metal_atom, 0, 1, 2, 3);
		delete [] upper_right_tetra_vertices;

		tetraMetalicStructure  upper_right_transformed=upper_rightTstruct.transform_tetrahedron();

		lower_left_tetra_vertices[0]=_octa_vertices_pl[_pl_vertex_order[0]];
		lower_left_tetra_vertices[1]=_octa_vertices_pl[_pl_vertex_order[1]];
		lower_left_tetra_vertices[2]=_octa_vertices_pl[_pl_vertex_order[3]];
		lower_left_tetra_vertices[3]=_octa_vertices_npl[1];
		lower_left_tetra_atoms[0]=_octa_atoms_pl[_pl_vertex_order[0]];
		lower_left_tetra_atoms[1]=_octa_atoms_pl[_pl_vertex_order[1]];
		lower_left_tetra_atoms[2]=_octa_atoms_pl[_pl_vertex_order[3]];
		lower_left_tetra_atoms[3]=_octa_atoms_npl[1];

		tetraMetalicStructure lower_leftTstruct=tetraMetalicStructure(lower_left_tetra_atoms, lower_left_tetra_vertices, _metal_atom, 0, 1, 2, 3);
		delete [] lower_left_tetra_vertices;

		tetraMetalicStructure  lower_left_transformed=lower_leftTstruct.transform_tetrahedron();

		lower_right_tetra_vertices[0]=_octa_vertices_pl[_pl_vertex_order[3]];
		lower_right_tetra_vertices[1]=_octa_vertices_pl[_pl_vertex_order[1]];
		lower_right_tetra_vertices[2]=_octa_vertices_pl[_pl_vertex_order[2]];
		lower_right_tetra_vertices[3]=_octa_vertices_npl[1];
		lower_right_tetra_atoms[0]=_octa_atoms_pl[_pl_vertex_order[3]];
		lower_right_tetra_atoms[1]=_octa_atoms_pl[_pl_vertex_order[1]];
		lower_right_tetra_atoms[2]=_octa_atoms_pl[_pl_vertex_order[2]];
		lower_right_tetra_atoms[3]=_octa_atoms_npl[1];

		tetraMetalicStructure lower_rightTstruct=tetraMetalicStructure(lower_right_tetra_atoms, lower_right_tetra_vertices, _metal_atom, 0, 1, 2, 3);
		delete [] lower_right_tetra_vertices;

		tetraMetalicStructure  lower_right_transformed=lower_rightTstruct.transform_tetrahedron();

		octaMetalicStructure tps=octaMetalicStructure(upper_left_transformed, upper_right_transformed, lower_left_transformed, lower_right_transformed, _metal_atom);

		return tps;
	}

	void octaMetalicStructure::measure_direction_cosines(){
		int k=0;
		double x, y, z, l, m, n, r;
		for (int i=0; i<_no_of_planar_atoms-1; i++){
			for (int j=i+1; j<_no_of_planar_atoms; j++){
				x=_octa_vertices_pl[j].X() - _octa_vertices_pl[i].X();
				y=_octa_vertices_pl[j].Y() - _octa_vertices_pl[i].Y();
				z=_octa_vertices_pl[j].Z() - _octa_vertices_pl[i].Z();
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
		for (int i=0; i<_no_of_non_planar_atoms; i++){
			for (int j=0; j<_no_of_planar_atoms; j++){

				x=_octa_vertices_pl[j].X() - _octa_vertices_npl[i].X();
				y=_octa_vertices_pl[j].Y() - _octa_vertices_npl[i].Y();
				z=_octa_vertices_pl[j].Z() - _octa_vertices_npl[i].Z();
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
	}

	void octaMetalicStructure::write_atoms_pdb_format(string fname){
		FILE* outf=fopen(fname.c_str(), "w");

		for (int i=0; i<_no_of_planar_atoms; i++){
			fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _octa_atoms_pl[i].get_atom_type().c_str(), _octa_atoms_pl[i].get_serial_no(), _octa_atoms_pl[i].get_atom_name().c_str(), _octa_atoms_pl[i].get_residue_name().c_str(), _octa_atoms_pl[i].get_chain_name().c_str(), _octa_atoms_pl[i].get_residue_no(), _octa_atoms_pl[i].X(), _octa_atoms_pl[i].Y(), _octa_atoms_pl[i].Z(), _octa_atoms_pl[i].get_occupancy(), _octa_atoms_pl[i].get_temp_factor(), _octa_atoms_pl[i].get_atom_symbol().c_str());
		}

		for (int i=0; i<_no_of_non_planar_atoms; i++){
			fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _octa_atoms_npl[i].get_atom_type().c_str(), _octa_atoms_npl[i].get_serial_no(), _octa_atoms_npl[i].get_atom_name().c_str(), _octa_atoms_npl[i].get_residue_name().c_str(), _octa_atoms_npl[i].get_chain_name().c_str(), _octa_atoms_npl[i].get_residue_no(), _octa_atoms_npl[i].X(), _octa_atoms_npl[i].Y(), _octa_atoms_npl[i].Z(), _octa_atoms_npl[i].get_occupancy(), _octa_atoms_npl[i].get_temp_factor(), _octa_atoms_npl[i].get_atom_symbol().c_str());
		}

		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[0].get_serial_no(), _octa_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[1].get_serial_no(), _octa_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[1].get_serial_no(), _octa_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[2].get_serial_no(), _octa_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[2].get_serial_no(), _octa_atoms_pl[3].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[3].get_serial_no(), _octa_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[3].get_serial_no(), _octa_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[0].get_serial_no(), _octa_atoms_pl[3].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_npl[0].get_serial_no(), _octa_atoms_pl[0].get_serial_no());
//		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[0].get_serial_no(), _octa_atoms_npl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_npl[0].get_serial_no(), _octa_atoms_pl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[1].get_serial_no(), _octa_atoms_npl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_npl[0].get_serial_no(), _octa_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[2].get_serial_no(), _octa_atoms_npl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_npl[0].get_serial_no(), _octa_atoms_pl[3].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[3].get_serial_no(), _octa_atoms_npl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_npl[1].get_serial_no(), _octa_atoms_pl[0].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[0].get_serial_no(), _octa_atoms_npl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_npl[1].get_serial_no(), _octa_atoms_pl[1].get_serial_no());
//		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[1].get_serial_no(), _octa_atoms_npl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_npl[1].get_serial_no(), _octa_atoms_pl[2].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[2].get_serial_no(), _octa_atoms_npl[1].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_npl[1].get_serial_no(), _octa_atoms_pl[3].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _octa_atoms_pl[3].get_serial_no(), _octa_atoms_npl[1].get_serial_no());
		fclose(outf);
	}

	void octaMetalicStructure::write_atoms_combined_tex_format_upper_left(octaMetalicStructure & comp, string fname){
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
		//atom** _comp_octa_atoms=comp.get_octa_atoms();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		tetrahedron left_comp=comp.get_upper_left_tetra();

		Point3D* left_vertices_standard=_upper_left_tetra.get_vertices();
		Point3D* left_vertices_comparable=left_comp.get_vertices();

		int* stand_vertex_order_left=_upper_left_tetra.get_vertex_order();
		int* comp_vertex_order_left=left_comp.get_vertex_order();

		double bond_len_standard_left=_upper_left_tetra.get_circum_radius();
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

		//Compute bond angles for Standard upper left Tetrahedron
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

		//Compute bond angles for Comparable upper left Tetrahedron
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
			//Write Standard Tetrahedron Upper left Vertices
			str="";

			if (_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()=="HOH"){
				st1=_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()=="O"){
					st1=_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()=="OXT"){
					st1=_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_left[p].c_str(), bond_angle_standard_left[p], bond_len_standard_left, str.c_str());

			//Write Comparable Tetrahedron Uper left Vertices
			str="";
			atom* comp_left_tetra_atoms=comp.get_upper_left_tetra_atoms();
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

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms Upper Left}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Octahedral Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Octahedral Site to Any Atom\\par\n";
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

	void octaMetalicStructure::write_atoms_combined_tex_format_upper_right(octaMetalicStructure & comp, string fname){
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
		//atom** _comp_octa_atoms=comp.get_octa_atoms();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		tetrahedron right_comp=comp.get_upper_right_tetra();

		Point3D* right_vertices_standard=_upper_right_tetra.get_vertices();
		Point3D* right_vertices_comparable=right_comp.get_vertices();

		int* stand_vertex_order_right=_upper_right_tetra.get_vertex_order();
		int* comp_vertex_order_right=right_comp.get_vertex_order();

		double bond_len_standard_right=_upper_right_tetra.get_circum_radius();
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
		st_right[3]="<";

		//Compute bond angles for Standard Upper right Tetrahedron
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

		//Compute bond angles for Comparable Upper right Tetrahedron
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
			//Write Standard Tetrahedron upper right Vertices
			str="";

			if (_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()=="HOH"){
				st1=_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()=="O"){
					st1=_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()=="OXT"){
					st1=_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_upper_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_right[p].c_str(), bond_angle_standard_right[p], bond_len_standard_right, str.c_str());

			//Write Comparable Tetrahedron upper right Vertices
			str="";
			atom* comp_right_tetra_atoms=comp.get_upper_right_tetra_atoms();
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

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms Upper Right}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Octahedral Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Octahedral Site to Any Atom\\par\n";
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

	void octaMetalicStructure::write_atoms_combined_tex_format_lower_left(octaMetalicStructure & comp, string fname){
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
		//atom** _comp_octa_atoms=comp.get_octa_atoms();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		tetrahedron left_comp=comp.get_lower_left_tetra();

		Point3D* left_vertices_standard=_lower_left_tetra.get_vertices();
		Point3D* left_vertices_comparable=left_comp.get_vertices();

		int* stand_vertex_order_left=_lower_left_tetra.get_vertex_order();
		int* comp_vertex_order_left=left_comp.get_vertex_order();

		double bond_len_standard_left=_lower_left_tetra.get_circum_radius();
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
		st_left[3]="<:";

		//Compute bond angles for Standard lower left Tetrahedron
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

		//Compute bond angles for Comparable lower left Tetrahedron
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
			//Write Standard Tetrahedron Upper left Vertices
			str="";

			if (_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()=="HOH"){
				st1=_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()=="O"){
					st1=_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()=="OXT"){
					st1=_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_left_tetra_atoms[stand_vertex_order_left[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_left[p].c_str(), bond_angle_standard_left[p], bond_len_standard_left, str.c_str());

			//Write Comparable Tetrahedron Uper left Vertices
			str="";
			atom* comp_left_tetra_atoms=comp.get_lower_left_tetra_atoms();
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

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms Lower Left}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Octahedral Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Octahedral Site to Any Atom\\par\n";
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

	void octaMetalicStructure::write_atoms_combined_tex_format_lower_right(octaMetalicStructure & comp, string fname){
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
		//atom** _comp_octa_atoms=comp.get_octa_atoms();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		tetrahedron right_comp=comp.get_lower_right_tetra();

		Point3D* right_vertices_standard=_lower_right_tetra.get_vertices();
		Point3D* right_vertices_comparable=right_comp.get_vertices();

		int* stand_vertex_order_right=_lower_right_tetra.get_vertex_order();
		int* comp_vertex_order_right=right_comp.get_vertex_order();

		double bond_len_standard_right=_lower_right_tetra.get_circum_radius();
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

		//Compute bond angles for Standard Upper right Tetrahedron
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

		//Compute bond angles for Comparable Upper right Tetrahedron
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
			//Write Standard Tetrahedron lower right Vertices
			str="";

			if (_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()=="HOH"){
				st1=_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()=="O"){
					st1=_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()=="OXT"){
					st1=_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_lower_right_tetra_atoms[stand_vertex_order_right[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_right[p].c_str(), bond_angle_standard_right[p], bond_len_standard_right, str.c_str());

			//Write Comparable Tetrahedron lower right Vertices
			str="";
			atom* comp_right_tetra_atoms=comp.get_lower_right_tetra_atoms();
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

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms Lower Right}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Octahedral Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Octahedral Site to Any Atom\\par\n";
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

	void octaMetalicStructure::compare(octaMetalicStructure & comp){
		double deviate_dist_upper_left[6];
		double deviate_dist_upper_right[6];
		double deviate_dist_lower_left[6];
		double deviate_dist_lower_right[6];

		cout<<"\nTetra V Order U L: "<<_upper_left_tetra.get_vertex_order()[0]<<", "<<_upper_left_tetra.get_vertex_order()[1]<<", "<<_upper_left_tetra.get_vertex_order()[2]<<", "<<_upper_left_tetra.get_vertex_order()[3]<<"; Comp: "<<comp.get_upper_left_tetra().get_vertex_order()[0]<<", "<<comp.get_upper_left_tetra().get_vertex_order()[1]<<", "<<comp.get_upper_left_tetra().get_vertex_order()[2]<<", "<<comp.get_upper_left_tetra().get_vertex_order()[3];

		_upper_left_tetra.deviation_distances(comp.get_upper_left_tetra(), _upper_left_tetra.measure_direction_cosines_of_sides(), deviate_dist_upper_left);
		_upper_right_tetra.deviation_distances(comp.get_upper_right_tetra(), _upper_right_tetra.measure_direction_cosines_of_sides(), deviate_dist_upper_right);
		_lower_left_tetra.deviation_distances(comp.get_lower_left_tetra(), _lower_left_tetra.measure_direction_cosines_of_sides(), deviate_dist_lower_left);
		_lower_right_tetra.deviation_distances(comp.get_lower_right_tetra(), _lower_right_tetra.measure_direction_cosines_of_sides(), deviate_dist_lower_right);

		int k=0;
		for (int p=0; p<6; p++){
			cout<<"\nUpper Left Dist:"<<deviate_dist_upper_left[p]<<", Upper Right Dist: "<<deviate_dist_upper_right[p]<<", Lower Left Dist: "<<deviate_dist_lower_left[p]<<", Lower Right Dist: "<<deviate_dist_lower_right[p];
			_devt_dist_vert_vert_upper_left[p]=deviate_dist_upper_left[p];
			_devt_dist_vert_vert_upper_right[p]=deviate_dist_upper_right[p];
			_devt_dist_vert_vert_lower_left[p]=deviate_dist_lower_left[p];
			_devt_dist_vert_vert_lower_right[p]=deviate_dist_lower_right[p];

			if (p!=3){
				_devt_dist_vert_vert[k]=_devt_dist_vert_vert_upper_left[p];
				k++;
			}
		}

		k=5;
		for (int p=0; p<6; p++){
			if ((p!=0) || (p!=2) || (p!=4)){
				_devt_dist_vert_vert[k]=_devt_dist_vert_vert_upper_right[p];
				k++;
			}
		}

		k=8;
		for (int p=0; p<6; p++){
			if ((p!=0) || (p!=1) || (p!=3)){
				_devt_dist_vert_vert[k]=_devt_dist_vert_vert_lower_right[p];
				k++;
			}
		}

		k=11;
		for (int p=0; p<6; p++){
			if (p==5){
				_devt_dist_vert_vert[k]=_devt_dist_vert_vert_lower_right[p];
				break;
			}
		}

		cout<<"\nCompare Angle Started";
		double deviate_angle_upper_left[12];
		double deviate_angle_upper_right[12];
		double deviate_angle_lower_left[12];
		double deviate_angle_lower_right[12];

		cout<<"\nDev. Angle Upper Left: ";
		_upper_left_tetra.deviation_angles_vertex_vertex(comp.get_upper_left_tetra(), deviate_angle_upper_left);

		cout<<"\nDev. Angle Upper Right: ";
		_upper_right_tetra.deviation_angles_vertex_vertex(comp.get_upper_right_tetra(), deviate_angle_upper_right);

		cout<<"\nDev. Angle Lower Left: ";
		_lower_left_tetra.deviation_angles_vertex_vertex(comp.get_lower_left_tetra(), deviate_angle_lower_left);

		cout<<"\nDev. Angle Lower Right: ";
		_lower_right_tetra.deviation_angles_vertex_vertex(comp.get_lower_right_tetra(), deviate_angle_lower_right);

		k=0;
		for (int p=0; p<12; p++){
			_devt_angle_vertex_vertex_upper_left[p]=deviate_angle_upper_left[p];
			_devt_angle_vertex_vertex_upper_right[p]=deviate_angle_upper_right[p];
			_devt_angle_vertex_vertex_lower_left[p]=deviate_angle_lower_left[p];
			_devt_angle_vertex_vertex_lower_right[p]=deviate_angle_lower_right[p];

			if ((p==0) || (p==1) || (p==2) || (p==4) || (p==7) || (p==9) || (p==10)){
				_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_upper_left[p];
				k++;
			}
			//cout<<"\nVert-Vert Angle IN COMPARE: "<<(p+1)<<": "<<_devt_angle_vertex_vertex_left[p]<<", "<<_devt_angle_vertex_vertex_right[p];
		}

		for (int p=0; p<12; p++){
			if ((p==2) || (p==5) || (p==6) || (p==7) || (p==8) || (p==10) || (p==11)){
				_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_upper_right[p];
				k++;
			}
		}

		_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_upper_left[4]+_devt_angle_vertex_vertex_upper_right[3];
		k++;
		_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_upper_left[6]+_devt_angle_vertex_vertex_upper_right[0];
		k++;

		for (int p=0; p<12; p++){
			if ((p==1) || (p==2) || (p==4) || (p==7) || (p==9) || (p==10)){
				_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_lower_left[p];
				k++;
			}
		}

		for (int p=0; p<12; p++){
			if ((p==2) || (p==5) || (p==7) || (p==8) || (p==10) || (p==11)){
				_devt_angle_vertex_vertex[k]=_devt_angle_vertex_vertex_lower_right[p];
				k++;
			}
		}

		cout<<"\nCompare MSE Started";
		compute_mse();

		_dist_devt_from_cc_upper_left=(_upper_left_tetra.get_circum_radius()-comp.get_upper_left_tetra().get_circum_radius());
		_dist_devt_from_cc_upper_right=(_upper_right_tetra.get_circum_radius()-comp.get_upper_right_tetra().get_circum_radius());
		_dist_devt_from_cc_lower_left=(_lower_left_tetra.get_circum_radius()-comp.get_lower_left_tetra().get_circum_radius());
		_dist_devt_from_cc_lower_right=(_lower_right_tetra.get_circum_radius()-comp.get_lower_right_tetra().get_circum_radius());

		if ((_dist_devt_from_cc_upper_left==0.0) && (_dist_devt_from_cc_upper_right==0.0) && (_dist_devt_from_cc_lower_left==0.0) && (_dist_devt_from_cc_lower_right==0.0)){
			_structure_size="same";
		}else if ((_dist_devt_from_cc_upper_left>0.0) && (_dist_devt_from_cc_upper_right>0.0) && (_dist_devt_from_cc_lower_left>0.0) && (_dist_devt_from_cc_lower_right>0.0)){
			_structure_size="smaller";
		}else if ((_dist_devt_from_cc_upper_left<0.0) && (_dist_devt_from_cc_upper_right<0.0) && (_dist_devt_from_cc_lower_left<0.0) && (_dist_devt_from_cc_lower_right<0.0)){
			_structure_size="larger";
		}else{
			_structure_size="irregular";
		}
	}

	void octaMetalicStructure::compute_mse(){
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

	void octaMetalicStructure::gen_report_metal_binding_sites(octaMetalicStructure & ori_stand, octaMetalicStructure & ori_comp, octaMetalicStructure & trans_comp, string mt, string standard_site_type, string comparable_site_type, FILE* fp){
		fprintf(fp, "HEADER     PROGRAM NAME: StructureDeviation  VERSION: 1.1");
		string str="TITLE      Deviation Report for Computation of Structural Deviations of Two Octahedral "+mt+" ION Binding Sites";
		fprintf(fp, "\n%s", str.c_str());

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Standard Octahedron:                   ");
		fprintf(fp, "\nREMARK     1  Molecule Type of Standard Octahedron: %s", standard_site_type.c_str());

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Distances of Atoms-Atoms:   s                      ");

		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_stand_pl_pl[i][j]=ori_stand.get_octa().get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_stand.get_octa_atoms_pl(), _vertex_vertex_distances_stand_pl_pl, 4);

		for (int i=0; i<8; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_stand_npl_pl[i][j]=ori_stand.get_octa().get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_stand.get_octa_atoms_pl(), ori_stand.get_octa_atoms_npl(), ori_comp.get_pl_vertex_order(), ori_comp.get_npl_vertex_order(), _vertex_vertex_distances_stand_npl_pl);

		fprintf(fp, "\nREMARK     1  Angles made by two atoms with another atom of Octahedron:");
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_stand_pl_pl[i][j]=ori_stand.get_octa().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_stand.get_octa_atoms_pl(), _vertex_vertex_angles_stand_pl_pl, 4);

		for (int i=0; i<24; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_stand_npl_pl[i][j]=ori_stand.get_octa().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_stand.get_octa_atoms_pl(), ori_stand.get_octa_atoms_npl(), ori_stand.get_pl_vertex_order(), ori_stand.get_npl_vertex_order(), _vertex_vertex_angles_stand_npl_pl);

		fprintf(fp, "\nREMARK     2");
		fprintf(fp, "\nREMARK     2  Comparable Tetrahedrons:");
		fprintf(fp, "\nREMARK     1  Molecule Type of Comparable Octahedron: %s", comparable_site_type.c_str());

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
				_vertex_vertex_distances_comp_pl_pl[i][j]=ori_comp.get_octa().get_vertex_vertex_distances_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_comp.get_octa_atoms_pl(), _vertex_vertex_distances_comp_pl_pl, 4);

		for (int i=0; i<8; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_comp_npl_pl[i][j]=ori_comp.get_octa().get_vertex_vertex_distances_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_comp.get_octa_atoms_pl(), ori_comp.get_octa_atoms_npl(), ori_comp.get_pl_vertex_order(), ori_comp.get_npl_vertex_order(), _vertex_vertex_distances_comp_npl_pl);

		fprintf(fp, "\nREMARK     2  Angles made by two atoms with another atom of Octahedron:");
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_comp_pl_pl[i][j]=ori_comp.get_octa().get_vertex_vertex_angles_pl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_comp.get_octa_atoms_pl(), _vertex_vertex_angles_comp_pl_pl, 4);

		for (int i=0; i<24; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_comp_npl_pl[i][j]=ori_comp.get_octa().get_vertex_vertex_angles_npl_pl()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_comp.get_octa_atoms_pl(), ori_comp.get_octa_atoms_npl(), ori_comp.get_pl_vertex_order(), ori_comp.get_npl_vertex_order(), _vertex_vertex_angles_comp_npl_pl);

		fprintf(fp, "\nREMARK     3");
		fprintf(fp, "\nREMARK     3  Deviation Results of Comparable Octahedron from Standard Octahedron:");

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

	void octaMetalicStructure::gen_report_vertex_vertex_distances(FILE* fp, atom atm[], double vert_vert_distances[][3], int no_of_sides){
		for (int i=0; i<no_of_sides; i++){
			fprintf(fp, "\nDISTRR %ld %ld %8.3lf", atm[(int) vert_vert_distances[i][0]].get_serial_no(), atm[(int) vert_vert_distances[i][1]].get_serial_no(), vert_vert_distances[i][2]);
		}
	}

	void octaMetalicStructure::gen_report_vertex_vertex_distances(FILE* fp, atom atm_pl[], atom atm_npl[], int* plv, int* npv, double vert_vert_distances[][3]){
		long int l, m;
		int p=0;
		for (int i=0; i<_no_of_non_planar_atoms; i++){
			for (int j=0; j<_no_of_planar_atoms; j++){
				l=atm_npl[npv[i]].get_serial_no();
				m=atm_pl[plv[j]].get_serial_no();
				fprintf(fp, "\nDISTRR %ld %ld %8.3lf", l, m, vert_vert_distances[p][2]);
				p++;
			}
		}
	}

	void octaMetalicStructure::gen_report_vertex_vertex_angles(FILE* fp, atom atm[], double vert_vert_angles[][4], int no_of_angles){
		for (int i=0; i<no_of_angles; i++){
			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", atm[(int) vert_vert_angles[i][0]].get_serial_no(), atm[(int) vert_vert_angles[i][1]].get_serial_no(), atm[(int) vert_vert_angles[i][2]].get_serial_no(), vert_vert_angles[i][3]);
		}
	}

	void octaMetalicStructure::gen_report_vertex_vertex_angles(FILE* fp, atom atm_pl[], atom atm_npl[], int* plv, int* npv, double vert_vert_angles[][4]){
		for (int p=0; p<6; p++){
			long int l, m, n;
			if (vert_vert_angles[p][0]==0){
				l= atm_pl[plv[0]].get_serial_no();
			}else if (vert_vert_angles[p][0]==1){
				l= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][0]==2){
				l= atm_pl[plv[3]].get_serial_no();
			}

			if (vert_vert_angles[p][1]==0){
				m= atm_pl[plv[0]].get_serial_no();
			}else if (vert_vert_angles[p][1]==1){
				m= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][1]==2){
				m= atm_pl[plv[3]].get_serial_no();
			}else if (vert_vert_angles[p][1]==3){
				m= atm_npl[npv[0]].get_serial_no();
			}

			if (vert_vert_angles[p][2]==1){
				n= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][2]==2){
				n= atm_pl[plv[3]].get_serial_no();
			}else if (vert_vert_angles[p][2]==3){
				n= atm_npl[npv[0]].get_serial_no();
			}

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
		}

		for (int p=6; p<12; p++){
			long int l, m, n;
			if (vert_vert_angles[p][0]==0){
				l= atm_pl[plv[3]].get_serial_no();
			}else if (vert_vert_angles[p][0]==1){
				l= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][0]==2){
				l= atm_pl[plv[2]].get_serial_no();
			}

			if (vert_vert_angles[p][1]==0){
				m= atm_pl[plv[3]].get_serial_no();
			}else if (vert_vert_angles[p][1]==1){
				m= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][1]==2){
				m= atm_pl[plv[2]].get_serial_no();
			}else if (vert_vert_angles[p][1]==3){
				m= atm_npl[npv[0]].get_serial_no();
			}

			if (vert_vert_angles[p][2]==1){
				n= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][2]==2){
				n= atm_pl[plv[2]].get_serial_no();
			}else if (vert_vert_angles[p][2]==3){
				n= atm_npl[npv[0]].get_serial_no();
			}

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
		}

		for (int p=12; p<18; p++){
			long int l, m, n;
			if (vert_vert_angles[p][0]==0){
				l= atm_pl[plv[0]].get_serial_no();
			}else if (vert_vert_angles[p][0]==1){
				l= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][0]==2){
				l= atm_pl[plv[3]].get_serial_no();
			}

			if (vert_vert_angles[p][1]==0){
				m= atm_pl[plv[0]].get_serial_no();
			}else if (vert_vert_angles[p][1]==1){
				m= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][1]==2){
				m= atm_pl[plv[3]].get_serial_no();
			}else if (vert_vert_angles[p][1]==3){
				m= atm_npl[npv[1]].get_serial_no();
			}

			if (vert_vert_angles[p][2]==1){
				n= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][2]==2){
				n= atm_pl[plv[3]].get_serial_no();
			}else if (vert_vert_angles[p][2]==3){
				n= atm_npl[npv[1]].get_serial_no();
			}

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
		}

		for (int p=18; p<24; p++){
			long int l, m, n;
			if (vert_vert_angles[p][0]==0){
				l= atm_pl[plv[3]].get_serial_no();
			}else if (vert_vert_angles[p][0]==1){
				l= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][0]==2){
				l= atm_pl[plv[2]].get_serial_no();
			}

			if (vert_vert_angles[p][1]==0){
				m= atm_pl[plv[3]].get_serial_no();
			}else if (vert_vert_angles[p][1]==1){
				m= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][1]==2){
				m= atm_pl[plv[2]].get_serial_no();
			}else if (vert_vert_angles[p][1]==3){
				m= atm_npl[npv[1]].get_serial_no();
			}

			if (vert_vert_angles[p][2]==1){
				n= atm_pl[plv[1]].get_serial_no();
			}else if (vert_vert_angles[p][2]==2){
				n= atm_pl[plv[2]].get_serial_no();
			}else if (vert_vert_angles[p][2]==3){
				n= atm_npl[npv[1]].get_serial_no();
			}

			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", l, m, n, vert_vert_angles[p][3]);
		}
	}
