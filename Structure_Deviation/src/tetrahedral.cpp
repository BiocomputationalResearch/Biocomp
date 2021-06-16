#include "tetrahedral.h"

void tetraMetalicStructure::print(atom mt, atom* non_mt){
	mt.fprint(stdout);

	for (int i=0; i<_no_of_binding_atoms; i++){
		non_mt[i].fprint(stdout);
	}
}

void tetraMetalicStructure::print(Point3D* vert){
	for (int i=0; i<_no_of_binding_atoms; i++){
		cout<<"\nVertex "<<i<<": "<<vert[i].X()<<"  "<<vert[i].Y()<<"  "<<vert[i].Z();
	}
}

tetraMetalicStructure::tetraMetalicStructure(){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetrahedron No. of Vertices should be 4!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=6;
	assert(_no_of_sides==6 && "For Tetrahedron No. of Sides should be 6!!");

	_no_of_angles_with_circum_centre=6;
	assert(_no_of_angles_with_circum_centre==6 && "For Tetrahedron No. of Angles of Two Vertices with Circumcentre should be 6!!");

	_no_of_vertex_vertex_angles=12;
	assert(_no_of_vertex_vertex_angles==12 && "For Tetrahedron No. of Angles of Two Vertices with Another Vertex should be 12!!");

	_no_of_features_to_compare=4;
	assert(_no_of_features_to_compare==4 && "For Tetrahedron No. of Features to Compare should be 4!!");

	_dist_devt_from_cc=0.0;
}

tetraMetalicStructure::tetraMetalicStructure(atom new_atm[], atom metal_atom){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetrahedron No. of Vertices should be 4!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=6;
	assert(_no_of_sides==6 && "For Tetrahedron No. of Sides should be 6!!");

	_no_of_angles_with_circum_centre=6;
	assert(_no_of_angles_with_circum_centre==6 && "For Tetrahedron No. of Angles of Two Vertices with Circumcentre should be 6!!");

	_no_of_vertex_vertex_angles=12;
	assert(_no_of_vertex_vertex_angles==12 && "For Tetrahedron No. of Angles of Two Vertices with Another Vertex should be 12!!");

	_no_of_features_to_compare=4;
	assert(_no_of_features_to_compare==4 && "For Tetrahedron No. of Features to Compare should be 4!!");

	_metal_atom=metal_atom;

	for (int i=0; i<_no_of_binding_atoms; i++){
		_tetra_atoms[i]=new_atm[i];
	}

	for (int i=0; i<_no_of_binding_atoms; i++){
		_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
	}

	_tetra=tetrahedron(_tetra_vertices);
	_tetra.compute_centre_of_mass();
	_tetra.measure_circum_centre_to_vertex_distances();
	_tetra.measure_angles_vertices_with_circum_centre();
	_tetra.measure_surface_areas_and_order_of_vertices();

	for (int i=0; i<_no_of_binding_atoms; i++){
		_vertex_order[i]=_tetra.get_vertex_order()[i];
	}


	print(_metal_atom, _tetra_atoms);

	_dist_devt_from_cc=0.0;
}

tetraMetalicStructure::tetraMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetrahedron No. of Vertices should be 4!!");
	assert(new_atm.size()==_no_of_binding_atoms && "Size of Atom Vector Should be 4");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=6;
	assert(_no_of_sides==6 && "For Tetrahedron No. of Sides should be 6!!");

	_no_of_angles_with_circum_centre=6;
	assert(_no_of_angles_with_circum_centre==6 && "For Tetrahedron No. of Angles of Two Vertices with Circumcentre should be 6!!");

	_no_of_vertex_vertex_angles=12;
	assert(_no_of_vertex_vertex_angles==12 && "For Tetrahedron No. of Angles of Two Vertices with Another Vertex should be 12!!");

	_no_of_features_to_compare=4;
	assert(_no_of_features_to_compare==4 && "For Tetrahedron No. of Features to Compare should be 4!!");

	_metal_atom=metal_atom;

	for (int i=0; i<new_atm.size(); i++){
		_tetra_atoms[i]=new_atm.at(i).get_atom_direct_contact();
	}

	for (int i=0; i<_no_of_binding_atoms; i++){
		_tetra_vertices[i]=Point3D(_tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z());
	}

	_tetra=tetrahedron(_tetra_vertices);
	_tetra.compute_centre_of_mass();
	_tetra.measure_circum_centre_to_vertex_distances();
	_tetra.measure_angles_vertices_with_circum_centre();
	_tetra.measure_surface_areas_and_order_of_vertices();

	for (int i=0; i<_no_of_binding_atoms; i++){
		_vertex_order[i]=_tetra.get_vertex_order()[i];
	}

	print(_metal_atom, _tetra_atoms);

	_dist_devt_from_cc=0.0;
}

tetraMetalicStructure::tetraMetalicStructure(atom new_atm[], Point3D new_vert[], atom metal_atom, int v1, int v2, int v3, int v4){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetrahedron No. of Vertices should be 4!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=6;
	assert(_no_of_sides==6 && "For Tetrahedron No. of Sides should be 6!!");

	_no_of_angles_with_circum_centre=6;
	assert(_no_of_angles_with_circum_centre==6 && "For Tetrahedron No. of Angles of Two Vertices with Circumcentre should be 6!!");

	_no_of_vertex_vertex_angles=12;
	assert(_no_of_vertex_vertex_angles==12 && "For Tetrahedron No. of Angles of Two Vertices with Another Vertex should be 12!!");

	_no_of_features_to_compare=4;
	assert(_no_of_features_to_compare==4 && "For Tetrahedron No. of Features to Compare should be 4!!");

	_metal_atom=metal_atom;

	for (int i=0; i<_no_of_binding_atoms; i++){
		_tetra_atoms[i]=new_atm[i];
	}

	for (int i=0; i<_no_of_binding_atoms; i++){
		_tetra_vertices[i]=new_vert[i];
	}

	_tetra=tetrahedron(_tetra_vertices);
	_tetra.compute_centre_of_mass();
	_tetra.measure_circum_centre_to_vertex_distances();
	_tetra.measure_angles_vertices_with_circum_centre();
	_tetra.set_vertex_order(v1, v2, v3, v4);

	for (int i=0; i<_no_of_binding_atoms; i++){
		_vertex_order[i]=_tetra.get_vertex_order()[i];
	}

	cout<<"\nReturn V Order: "<<_vertex_order[0]<<", "<<_vertex_order[1]<<", "<<_vertex_order[2]<<", "<<_vertex_order[3];
	print(_metal_atom, _tetra_atoms);

	_dist_devt_from_cc=0.0;
}

tetraMetalicStructure tetraMetalicStructure::apply_transformations(){
	make_base_parallel_xy();
	check_and_change_orientation();
	tetraMetalicStructure tms=shift_circum_centre_to_origin();
	return tms;
}

tetraMetalicStructure tetraMetalicStructure::transform_tetrahedron(){
	cout<<"\n Original Vertices";
	print(_tetra_vertices);

	make_base_parallel_xy();
	change_orientation();
	cout<<"\nVERTEX ORDER After Orientation: "<<_vertex_order[0]<<"  "<<_vertex_order[1]<<"  "<<_vertex_order[2]<<"  "<<_vertex_order[3];
	tetraMetalicStructure tms=shift_circum_centre_to_origin(_vertex_order);
	return tms;
}

atom tetraMetalicStructure::get_metal(){
	return _metal_atom;
}

atom* tetraMetalicStructure::get_tetra_atoms(){
	return _tetra_atoms;
}

tetrahedron tetraMetalicStructure::get_tetrahedron(){
	return _tetra;
}

int tetraMetalicStructure::get_no_of_binding_atoms(){
	return _no_of_binding_atoms;
}

double* tetraMetalicStructure::get_side_deviations(){
	return _devt_dist_vert_vert;
}

double tetraMetalicStructure::get_distance_deviations_from_cc(){
	return _dist_devt_from_cc;
}

double* tetraMetalicStructure::get_angle_deviations_with_cc(){
	return _devt_angle_with_cc;
}

double* tetraMetalicStructure::get_angle_deviations_vertex_vertex(){
	return _devt_angle_vertex_vertex;
}

double* tetraMetalicStructure::get_mse(){
	return _mse;
}

bool tetraMetalicStructure::get_nearer(){
	return _nearer;
}

string tetraMetalicStructure::get_structure_size(){
	return _structure_size;
}

	void tetraMetalicStructure::comparison(tetraMetalicStructure & ctdrn, string out_comparable_name){
		double dist_stand[6][3];

		for (int i=0; i<6; i++){
			for (int j=0; j<3; j++){
				dist_stand[i][j]=_tetra.get_vertex_vertex_distances()[i][j];
			}
		}

		measure_direction_cosines();

		int k=0;
		double proj[_no_of_sides];

		for (int i=0; i<(_no_of_binding_atoms-1); i++){
			for (int j=i+1; j<_no_of_binding_atoms; j++){
				double prx=	_direction_cosines[k][0]*(ctdrn.get_tetrahedron().get_vertices()[j].X()-ctdrn.get_tetrahedron().get_vertices()[i].X());
				double pry= _direction_cosines[k][1]*(ctdrn.get_tetrahedron().get_vertices()[j].Y()-ctdrn.get_tetrahedron().get_vertices()[i].Y());
				double prz=	_direction_cosines[k][2]*(ctdrn.get_tetrahedron().get_vertices()[j].Z()-ctdrn.get_tetrahedron().get_vertices()[i].Z());

				proj[k]=(prx+pry+prz);

				k++;
			}
		}

		double dist_dev_vert_vert[_no_of_sides];
		FILE* outf=fopen(out_comparable_name.c_str(), "w");

		cout<<"\n\nComparison Dist:\n";
		for (int i=0; i<_no_of_sides; i++){
			double diff=dist_stand[i][2]-proj[i];
			dist_dev_vert_vert[i]=diff;
			cout<<"Side "<<(i+1)<<": "<<dist_dev_vert_vert[i]<<endl;
		}

		double sum=0.0;
		for (int p=0; p<_no_of_sides; p++){
			sum+=(dist_dev_vert_vert[p] * dist_dev_vert_vert[p]);
		}

		double mse_dist=(sum/_no_of_sides);
		cout<<"\n MSE Distance: "<<mse_dist;
		fprintf(outf, "\nOriginal MSE Dist: %8.3lf", mse_dist);

		cout<<"\n\nComparison Angle:\n";
		double angle_stand[12][4];
		double angle_comp[12][4];
		double angle_dev_vert_vert[_no_of_vertex_vertex_angles];

		for (int i=0; i<12; i++){
			for (int j=0; j<4; j++){
				angle_stand[i][j]=_tetra.get_vertex_vertex_angles()[i][j];
				angle_comp[i][j]=ctdrn.get_tetrahedron().get_vertex_vertex_angles()[i][j];
			}
		}

		for (int p=0; p<_no_of_vertex_vertex_angles; p++){
			angle_dev_vert_vert[p]=angle_stand[p][3]-angle_comp[p][3];
			cout<<"\nAngle "<<(p+1)<<": "<<angle_dev_vert_vert[p];
		}

		sum=0.0;
		for (int p=0; p<_no_of_vertex_vertex_angles; p++){
			sum+=(angle_dev_vert_vert[p] * angle_dev_vert_vert[p]);
		}

		double mse_angle=(sum/_no_of_vertex_vertex_angles);
		cout<<"\n MSE Angle: "<<mse_angle;
		fprintf(outf, "\nMSE Angle: %8.3lf", mse_angle);

		fclose(outf);
	}

	void tetraMetalicStructure::make_base_parallel_xy(){
		for (int i=0; i<_no_of_binding_atoms; i++){
			_vertex_order[i]=_tetra.get_vertex_order()[i];
		}

		cout<<"\nVERTEX ORDER: "<<_vertex_order[0]<<"  "<<_vertex_order[1]<<"  "<<_vertex_order[2]<<"  "<<_vertex_order[3];
		shift_1st_vertex_to_origin();

		cout<<"\n\nShifted Vertices to Origin (1st Vertex):";
		print(_shifted_vertices);
		cout<<"\n\nShifted Metal Parallel:";
		_shifted_metal.fprint(stdout);

		//Rotate largest side of base (2nd vertex) to make it parallel to x-axis
		rotate_largest_side();
		cout<<"\n\nRotated Vertices (Largest Side):";
		print(_rotated_vertices_largest_side);
		cout<<"\nRotated Metal Largest Side: ";
		_rotated_metal_yaxis.fprint(stdout);

		//Rotate 3rd vertex of tetrahedron to make it parallel to x-axis
		rotate_3rd_vertex();
		cout<<"\n\nParallel Vertices (After Rotating 3rd Vertex):";
		print(_parallel_vertices);
		cout<<"\n\nParallel Atoms:";
		for (int i=0; i<_no_of_binding_atoms; i++){
			_parallel_tetra_atoms[i]=atom(_tetra_atoms[i].get_atom_type(), _tetra_atoms[i].get_serial_no(), _tetra_atoms[i].get_atom_name(), _tetra_atoms[i].get_alt_conf_indicator(), _tetra_atoms[i].get_residue_name(), _tetra_atoms[i].get_chain_name(), _tetra_atoms[i].get_residue_no(), _parallel_vertices[i].X(), _parallel_vertices[i].Y(), _parallel_vertices[i].Z(), _tetra_atoms[i].get_occupancy(), _tetra_atoms[i].get_temp_factor(), 0.0, _tetra_atoms[i].get_atom_symbol(), _tetra_atoms[i].get_partial_charge());
		}

		//Print Atoms
		print(_parallel_metal_atom, _parallel_tetra_atoms);
	}

	void tetraMetalicStructure::shift_1st_vertex_to_origin(){
		cout<<"\nVERTEX ORDER: "<<_vertex_order[0]<<"  "<<_vertex_order[1]<<"  "<<_vertex_order[2]<<"  "<<_vertex_order[3];
		Point3D shift_amount=_tetra_vertices[_vertex_order[0]].shift_to_origin();
		cout<<"\nShift Amount: "<<shift_amount.X()<<"  "<<shift_amount.Y()<<"  "<<shift_amount.Z();
		//Shift Metal Atom and Vertices
		_shifted_metal=atom(_metal_atom.get_atom_type(), _metal_atom.get_serial_no(), _metal_atom.get_atom_name(), _metal_atom.get_alt_conf_indicator(), _metal_atom.get_residue_name(), _metal_atom.get_chain_name(), _metal_atom.get_residue_no(), (_metal_atom.X()+shift_amount.X()), (_metal_atom.Y()+shift_amount.Y()), (_metal_atom.Z()+shift_amount.Z()), _metal_atom.get_occupancy(), _metal_atom.get_temp_factor(), 0.0, _metal_atom.get_atom_symbol(), _metal_atom.get_partial_charge());
		for (int i=0; i<_no_of_binding_atoms; i++){
			_shifted_vertices[i]=Point3D((_tetra_vertices[i].X()+shift_amount.X()), (_tetra_vertices[i].Y()+shift_amount.Y()), (_tetra_vertices[i].Z()+shift_amount.Z()));
		}
	}

	void tetraMetalicStructure::rotate_largest_side(){
		double rotate_angle=_shifted_vertices[_vertex_order[1]].rotate_about_z_axis();

		//Rotate Metal and Vertices about z-axis
		double _new_x_coord=_shifted_metal.X()*cos(rotate_angle) - _shifted_metal.Y()*sin(rotate_angle);
		double _new_y_coord=_shifted_metal.X()*sin(rotate_angle) + _shifted_metal.Y()*cos(rotate_angle);
		double _new_z_coord=_shifted_metal.Z();
		_rotated_metal_zaxis=atom(_shifted_metal.get_atom_type(), _shifted_metal.get_serial_no(), _shifted_metal.get_atom_name(), _shifted_metal.get_alt_conf_indicator(), _shifted_metal.get_residue_name(), _shifted_metal.get_chain_name(), _shifted_metal.get_residue_no(), _new_x_coord, _new_y_coord, _new_z_coord, _shifted_metal.get_occupancy(), _shifted_metal.get_temp_factor(), 0.0, _shifted_metal.get_atom_symbol(), _shifted_metal.get_partial_charge());

		Point3D rotated_vertex[_no_of_binding_atoms];
		for (int i=0; i<_no_of_binding_atoms; i++){	//Rotate about z-axis
			_new_x_coord=_shifted_vertices[i].X()*cos(rotate_angle) - _shifted_vertices[i].Y()*sin(rotate_angle);
			_new_y_coord=_shifted_vertices[i].X()*sin(rotate_angle) + _shifted_vertices[i].Y()*cos(rotate_angle);
			_new_z_coord=_shifted_vertices[i].Z();

			rotated_vertex[i]=Point3D(_new_x_coord, _new_y_coord, _new_z_coord);
		}

		rotate_angle=rotated_vertex[_vertex_order[1]].rotate_about_y_axis();
		//Rotate Metal about y-axis
		_new_x_coord=_rotated_metal_zaxis.X()*cos(rotate_angle) - _rotated_metal_zaxis.Z()*sin(rotate_angle);
		_new_y_coord=_rotated_metal_zaxis.Y();
		_new_z_coord=_rotated_metal_zaxis.X()*sin(rotate_angle) + _rotated_metal_zaxis.Z()*cos(rotate_angle);
		_rotated_metal_yaxis=atom(_rotated_metal_zaxis.get_atom_type(), _rotated_metal_zaxis.get_serial_no(), _rotated_metal_zaxis.get_atom_name(), _rotated_metal_zaxis.get_alt_conf_indicator(), _rotated_metal_zaxis.get_residue_name(), _rotated_metal_zaxis.get_chain_name(), _rotated_metal_zaxis.get_residue_no(), _new_x_coord, _new_y_coord, _new_z_coord, _rotated_metal_zaxis.get_occupancy(), _rotated_metal_zaxis.get_temp_factor(), 0.0, _rotated_metal_zaxis.get_atom_symbol(), _rotated_metal_zaxis.get_partial_charge());

		for (int i=0; i<_no_of_binding_atoms; i++){	//Rotate about y-axis
			_new_x_coord=rotated_vertex[i].X()*cos(rotate_angle) - rotated_vertex[i].Z()*sin(rotate_angle);
			_new_y_coord=rotated_vertex[i].Y();
			_new_z_coord=rotated_vertex[i].X()*sin(rotate_angle) + rotated_vertex[i].Z()*cos(rotate_angle);

			_rotated_vertices_largest_side[i]=Point3D(_new_x_coord, _new_y_coord, _new_z_coord);
		}
	}

	void tetraMetalicStructure::rotate_3rd_vertex(){
		double rotate_angle=_rotated_vertices_largest_side[_vertex_order[2]].rotate_about_x_axis();
		//Rotate Metal about x-axis
		double _new_x_coord=_rotated_metal_yaxis.X();
		double _new_y_coord=_rotated_metal_yaxis.Y()*cos(rotate_angle) - _rotated_metal_yaxis.Z()*sin(rotate_angle);
		double _new_z_coord=_rotated_metal_yaxis.Y()*sin(rotate_angle) + _rotated_metal_yaxis.Z()*cos(rotate_angle);
		_parallel_metal_atom=atom(_rotated_metal_yaxis.get_atom_type(), _rotated_metal_yaxis.get_serial_no(), _rotated_metal_yaxis.get_atom_name(), _rotated_metal_yaxis.get_alt_conf_indicator(), _rotated_metal_yaxis.get_residue_name(), _rotated_metal_yaxis.get_chain_name(), _rotated_metal_yaxis.get_residue_no(), _new_x_coord, _new_y_coord, _new_z_coord, _rotated_metal_yaxis.get_occupancy(), _rotated_metal_yaxis.get_temp_factor(), 0.0, _rotated_metal_yaxis.get_atom_symbol(), _rotated_metal_yaxis.get_partial_charge());

		for (int i=0; i<_no_of_binding_atoms; i++){
			_new_x_coord=_rotated_vertices_largest_side[i].X();
			_new_y_coord=_rotated_vertices_largest_side[i].Y()*cos(rotate_angle) - _rotated_vertices_largest_side[i].Z()*sin(rotate_angle);
			_new_z_coord=_rotated_vertices_largest_side[i].Y()*sin(rotate_angle) + _rotated_vertices_largest_side[i].Z()*cos(rotate_angle);

			_parallel_vertices[i]=Point3D(_new_x_coord, _new_y_coord, _new_z_coord);
		}
	}

	bool tetraMetalicStructure::check_distance(Point3D vert[], int _vertex_order[]){
		if (vert[_vertex_order[0]].dist(&vert[_vertex_order[2]]) >= vert[_vertex_order[1]].dist(&vert[_vertex_order[2]])){
			return true;
		}else{
			return false;
		}
	}

	void tetraMetalicStructure::check_and_change_orientation(){
		bool flag=false, flag1=false;
		Point3D* temp_vertices;

		if (_parallel_vertices[_vertex_order[2]].Y()>0){	//check for vertex3
			flag=(_parallel_vertices[_vertex_order[3]].Z()>=0)?true:false;	//Check Fourth Vertex

			if (flag==true){
				flag1=check_distance(_parallel_vertices, _vertex_order);	//check_distance_of_larger_and_large_base_edges
				if (flag1==true){
					_oriented_metal_atom=atom(_parallel_metal_atom.get_atom_type(), _parallel_metal_atom.get_serial_no(), _parallel_metal_atom.get_atom_name(), _parallel_metal_atom.get_alt_conf_indicator(), _parallel_metal_atom.get_residue_name(), _parallel_metal_atom.get_chain_name(), _parallel_metal_atom.get_residue_no(), _parallel_metal_atom.X(), _parallel_metal_atom.Y(), _parallel_metal_atom.Z(), _parallel_metal_atom.get_occupancy(), _parallel_metal_atom.get_temp_factor(), 0.0, _parallel_metal_atom.get_atom_symbol(), _parallel_metal_atom.get_partial_charge());

					for (int p=0; p<_no_of_binding_atoms; p++){
						_oriented_vertices[p]=Point3D(_parallel_vertices[p].X(), _parallel_vertices[p].Y(), _parallel_vertices[p].Z());
					}
				}else{
					//Reflect about y-axis and Shift 2nd vertex to origin and update vertex order
					reflect_and_shift_2nd_vertex_update_vertex_order(_parallel_vertices, _parallel_metal_atom);
				}
			}else{
				//Reflect 4th Vertex about the XY-plane
				atom _reflected_metal_atom=reflect_4th_vertex_about_xy_plane(_parallel_vertices, _parallel_metal_atom);

				flag1=check_distance(_reflected_vertices_4th_vertex, _vertex_order);
				if (flag1==true){
					_oriented_metal_atom=atom(_reflected_metal_atom.get_atom_type(), _reflected_metal_atom.get_serial_no(), _reflected_metal_atom.get_atom_name(), _reflected_metal_atom.get_alt_conf_indicator(), _reflected_metal_atom.get_residue_name(), _reflected_metal_atom.get_chain_name(), _reflected_metal_atom.get_residue_no(), _reflected_metal_atom.X(), _reflected_metal_atom.Y(), _reflected_metal_atom.Z(), _reflected_metal_atom.get_occupancy(), _reflected_metal_atom.get_temp_factor(), 0.0, _reflected_metal_atom.get_atom_symbol(), _reflected_metal_atom.get_partial_charge());
					for (int i=0; i<_no_of_binding_atoms; i++){
						_oriented_vertices[i]=Point3D(_reflected_vertices_4th_vertex[i].X(), _reflected_vertices_4th_vertex[i].Y(), _reflected_vertices_4th_vertex[i].Z());
					}
				}else{
					//Reflect about y-axis and Shift 2nd vertex to origin and update vertex order
					reflect_and_shift_2nd_vertex_update_vertex_order(_reflected_vertices_4th_vertex, _reflected_metal_atom);
				}
			}
		}else{
			//Reflect 3rd vertex
			Point3D metal_point=Point3D(_parallel_metal_atom.X(), _parallel_metal_atom.Y(), _parallel_metal_atom.Z());
			Point3D reflected_metal_point=metal_point.reflect_about_x_axis();
			atom _reflected_metal_atom_3rd_vertex=atom(_parallel_metal_atom.get_atom_type(), _parallel_metal_atom.get_serial_no(), _parallel_metal_atom.get_atom_name(), _parallel_metal_atom.get_alt_conf_indicator(), _parallel_metal_atom.get_residue_name(), _parallel_metal_atom.get_chain_name(), _parallel_metal_atom.get_residue_no(), reflected_metal_point.X(), reflected_metal_point.Y(), reflected_metal_point.Z(), _parallel_metal_atom.get_occupancy(), _parallel_metal_atom.get_temp_factor(), 0.0, _parallel_metal_atom.get_atom_symbol(), _parallel_metal_atom.get_partial_charge());
			for (int p=0; p<_no_of_binding_atoms; p++){
				_reflected_vertices_3rd_vertex[p]=_parallel_vertices[p].reflect_about_x_axis();
			}

			flag=(_reflected_vertices_3rd_vertex[_vertex_order[3]].Z()>=0)?true:false;	//Check Fourth Vertex

			if (flag==true){
				flag1=check_distance(_reflected_vertices_3rd_vertex, _vertex_order);
				if (flag1==true){
					_oriented_metal_atom=atom(_reflected_metal_atom_3rd_vertex.get_atom_type(), _reflected_metal_atom_3rd_vertex.get_serial_no(), _reflected_metal_atom_3rd_vertex.get_atom_name(), _reflected_metal_atom_3rd_vertex.get_alt_conf_indicator(), _reflected_metal_atom_3rd_vertex.get_residue_name(), _reflected_metal_atom_3rd_vertex.get_chain_name(), _reflected_metal_atom_3rd_vertex.get_residue_no(), _reflected_metal_atom_3rd_vertex.X(), _reflected_metal_atom_3rd_vertex.Y(), _reflected_metal_atom_3rd_vertex.Z(), _reflected_metal_atom_3rd_vertex.get_occupancy(), _reflected_metal_atom_3rd_vertex.get_temp_factor(), 0.0, _reflected_metal_atom_3rd_vertex.get_atom_symbol(), _reflected_metal_atom_3rd_vertex.get_partial_charge());
					for (int p=0; p<_no_of_binding_atoms; p++){
						_oriented_vertices[p]=Point3D(_reflected_vertices_3rd_vertex[p].X(), _reflected_vertices_3rd_vertex[p].Y(), _reflected_vertices_3rd_vertex[p].Z());
					}
				}else{
					//Reflect about y-axis and Shift 2nd vertex to origin and update vertex order
					reflect_and_shift_2nd_vertex_update_vertex_order(_reflected_vertices_3rd_vertex, _reflected_metal_atom_3rd_vertex);
				}
			}else{
				//Reflect 4th Vertex about the XY-plane
				atom _reflected_metal_atom=reflect_4th_vertex_about_xy_plane(_reflected_vertices_3rd_vertex, _reflected_metal_atom_3rd_vertex);

				flag1=check_distance(_reflected_vertices_4th_vertex, _vertex_order);
				if (flag1==true){
					_oriented_metal_atom=atom(_reflected_metal_atom.get_atom_type(), _reflected_metal_atom.get_serial_no(), _reflected_metal_atom.get_atom_name(), _reflected_metal_atom.get_alt_conf_indicator(), _reflected_metal_atom.get_residue_name(), _reflected_metal_atom.get_chain_name(), _reflected_metal_atom.get_residue_no(), _reflected_metal_atom.X(), _reflected_metal_atom.Y(), _reflected_metal_atom.Z(), _reflected_metal_atom.get_occupancy(), _reflected_metal_atom.get_temp_factor(), 0.0, _reflected_metal_atom.get_atom_symbol(), _reflected_metal_atom.get_partial_charge());
					for (int p=0; p<_no_of_binding_atoms; p++){
						_oriented_vertices[p]=Point3D(_reflected_vertices_4th_vertex[p].X(), _reflected_vertices_4th_vertex[p].Y(), _reflected_vertices_4th_vertex[p].Z());
					}
				}else{
					//Reflect about y-axis and Shift 2nd vertex to origin and update vertex order
					reflect_and_shift_2nd_vertex_update_vertex_order(_reflected_vertices_4th_vertex, _reflected_metal_atom);
				}
			}
		}

		cout<<"\n\nOriented Vertices:";
		print(_oriented_vertices);
		cout<<"\n\nOriented Metal:";
		_oriented_metal_atom.fprint(stdout);

		//Generate atoms corresponding to oriented vertices
		for (int i=0; i<_no_of_binding_atoms; i++){
			_oriented_tetra_atoms[i]=atom(_parallel_tetra_atoms[i].get_atom_type(), _parallel_tetra_atoms[i].get_serial_no(), _parallel_tetra_atoms[i].get_atom_name(), _parallel_tetra_atoms[i].get_alt_conf_indicator(), _parallel_tetra_atoms[i].get_residue_name(), _parallel_tetra_atoms[i].get_chain_name(), _parallel_tetra_atoms[i].get_residue_no(), _oriented_vertices[i].X(), _oriented_vertices[i].Y(), _oriented_vertices[i].Z(), _parallel_tetra_atoms[i].get_occupancy(), _parallel_tetra_atoms[i].get_temp_factor(), 0.0, _parallel_tetra_atoms[i].get_atom_symbol(), _parallel_tetra_atoms[i].get_partial_charge());
		}
	}

	void tetraMetalicStructure::change_orientation(){
		bool flag=false, flag1=false;
		Point3D* temp_vertices;

		if (_parallel_vertices[_vertex_order[2]].Y()>0){	//check for vertex3
			flag=(_parallel_vertices[_vertex_order[3]].Z()>=0)?true:false;	//Check Fourth Vertex

			if (flag==true){
				_oriented_metal_atom=atom(_parallel_metal_atom.get_atom_type(), _parallel_metal_atom.get_serial_no(), _parallel_metal_atom.get_atom_name(), _parallel_metal_atom.get_alt_conf_indicator(), _parallel_metal_atom.get_residue_name(), _parallel_metal_atom.get_chain_name(), _parallel_metal_atom.get_residue_no(), _parallel_metal_atom.X(), _parallel_metal_atom.Y(), _parallel_metal_atom.Z(), _parallel_metal_atom.get_occupancy(), _parallel_metal_atom.get_temp_factor(), 0.0, _parallel_metal_atom.get_atom_symbol(), _parallel_metal_atom.get_partial_charge());

				for (int p=0; p<_no_of_binding_atoms; p++){
					_oriented_vertices[p]=Point3D(_parallel_vertices[p].X(), _parallel_vertices[p].Y(), _parallel_vertices[p].Z());
				}
			}else{
				//Reflect 4th Vertex about the XY-plane
				atom _reflected_metal_atom=reflect_4th_vertex_about_xy_plane(_parallel_vertices, _parallel_metal_atom);

				_oriented_metal_atom=atom(_reflected_metal_atom.get_atom_type(), _reflected_metal_atom.get_serial_no(), _reflected_metal_atom.get_atom_name(), _reflected_metal_atom.get_alt_conf_indicator(), _reflected_metal_atom.get_residue_name(), _reflected_metal_atom.get_chain_name(), _reflected_metal_atom.get_residue_no(), _reflected_metal_atom.X(), _reflected_metal_atom.Y(), _reflected_metal_atom.Z(), _reflected_metal_atom.get_occupancy(), _reflected_metal_atom.get_temp_factor(), 0.0, _reflected_metal_atom.get_atom_symbol(), _reflected_metal_atom.get_partial_charge());
				for (int i=0; i<_no_of_binding_atoms; i++){
					_oriented_vertices[i]=Point3D(_reflected_vertices_4th_vertex[i].X(), _reflected_vertices_4th_vertex[i].Y(), _reflected_vertices_4th_vertex[i].Z());
				}
			}
		}else{
			//Reflect 3rd vertex
			Point3D metal_point=Point3D(_parallel_metal_atom.X(), _parallel_metal_atom.Y(), _parallel_metal_atom.Z());
			Point3D reflected_metal_point=metal_point.reflect_about_x_axis();
			atom _reflected_metal_atom_3rd_vertex=atom(_parallel_metal_atom.get_atom_type(), _parallel_metal_atom.get_serial_no(), _parallel_metal_atom.get_atom_name(), _parallel_metal_atom.get_alt_conf_indicator(), _parallel_metal_atom.get_residue_name(), _parallel_metal_atom.get_chain_name(), _parallel_metal_atom.get_residue_no(), reflected_metal_point.X(), reflected_metal_point.Y(), reflected_metal_point.Z(), _parallel_metal_atom.get_occupancy(), _parallel_metal_atom.get_temp_factor(), 0.0, _parallel_metal_atom.get_atom_symbol(), _parallel_metal_atom.get_partial_charge());
			for (int p=0; p<_no_of_binding_atoms; p++){
				_reflected_vertices_3rd_vertex[p]=_parallel_vertices[p].reflect_about_x_axis();
			}

			flag=(_reflected_vertices_3rd_vertex[_vertex_order[3]].Z()>=0)?true:false;	//Check Fourth Vertex

			if (flag==true){
				_oriented_metal_atom=atom(_reflected_metal_atom_3rd_vertex.get_atom_type(), _reflected_metal_atom_3rd_vertex.get_serial_no(), _reflected_metal_atom_3rd_vertex.get_atom_name(), _reflected_metal_atom_3rd_vertex.get_alt_conf_indicator(), _reflected_metal_atom_3rd_vertex.get_residue_name(), _reflected_metal_atom_3rd_vertex.get_chain_name(), _reflected_metal_atom_3rd_vertex.get_residue_no(), _reflected_metal_atom_3rd_vertex.X(), _reflected_metal_atom_3rd_vertex.Y(), _reflected_metal_atom_3rd_vertex.Z(), _reflected_metal_atom_3rd_vertex.get_occupancy(), _reflected_metal_atom_3rd_vertex.get_temp_factor(), 0.0, _reflected_metal_atom_3rd_vertex.get_atom_symbol(), _reflected_metal_atom_3rd_vertex.get_partial_charge());
				for (int p=0; p<_no_of_binding_atoms; p++){
					_oriented_vertices[p]=Point3D(_reflected_vertices_3rd_vertex[p].X(), _reflected_vertices_3rd_vertex[p].Y(), _reflected_vertices_3rd_vertex[p].Z());
				}
			}else{
				//Reflect 4th Vertex about the XY-plane
				atom _reflected_metal_atom=reflect_4th_vertex_about_xy_plane(_reflected_vertices_3rd_vertex, _reflected_metal_atom_3rd_vertex);

				_oriented_metal_atom=atom(_reflected_metal_atom.get_atom_type(), _reflected_metal_atom.get_serial_no(), _reflected_metal_atom.get_atom_name(), _reflected_metal_atom.get_alt_conf_indicator(), _reflected_metal_atom.get_residue_name(), _reflected_metal_atom.get_chain_name(), _reflected_metal_atom.get_residue_no(), _reflected_metal_atom.X(), _reflected_metal_atom.Y(), _reflected_metal_atom.Z(), _reflected_metal_atom.get_occupancy(), _reflected_metal_atom.get_temp_factor(), 0.0, _reflected_metal_atom.get_atom_symbol(), _reflected_metal_atom.get_partial_charge());
				for (int p=0; p<_no_of_binding_atoms; p++){
					_oriented_vertices[p]=Point3D(_reflected_vertices_4th_vertex[p].X(), _reflected_vertices_4th_vertex[p].Y(), _reflected_vertices_4th_vertex[p].Z());
				}
			}
		}

		cout<<"\n\nOriented Vertices:";
		print(_oriented_vertices);
		cout<<"\n\nOriented Metal:";
		_oriented_metal_atom.fprint(stdout);

		//Generate atoms corresponding to oriented vertices
		for (int i=0; i<_no_of_binding_atoms; i++){
			_oriented_tetra_atoms[i]=atom(_parallel_tetra_atoms[i].get_atom_type(), _parallel_tetra_atoms[i].get_serial_no(), _parallel_tetra_atoms[i].get_atom_name(), _parallel_tetra_atoms[i].get_alt_conf_indicator(), _parallel_tetra_atoms[i].get_residue_name(), _parallel_tetra_atoms[i].get_chain_name(), _parallel_tetra_atoms[i].get_residue_no(), _oriented_vertices[i].X(), _oriented_vertices[i].Y(), _oriented_vertices[i].Z(), _parallel_tetra_atoms[i].get_occupancy(), _parallel_tetra_atoms[i].get_temp_factor(), 0.0, _parallel_tetra_atoms[i].get_atom_symbol(), _parallel_tetra_atoms[i].get_partial_charge());
		}
	}

	void tetraMetalicStructure::reflect_and_shift_2nd_vertex_update_vertex_order(Point3D vertices[], atom metal_atom){
		Point3D reflected_vertices[_no_of_binding_atoms];
		Point3D metal_point=Point3D(metal_atom.X(), metal_atom.Y(), metal_atom.Z());
		Point3D reflected_metal_point=metal_point.reflect_about_y_axis();
		atom _reflected_metal_atom=atom(metal_atom.get_atom_type(), metal_atom.get_serial_no(), metal_atom.get_atom_name(), metal_atom.get_alt_conf_indicator(), metal_atom.get_residue_name(), metal_atom.get_chain_name(), metal_atom.get_residue_no(), reflected_metal_point.X(), reflected_metal_point.Y(), reflected_metal_point.Z(), metal_atom.get_occupancy(), metal_atom.get_temp_factor(), 0.0, metal_atom.get_atom_symbol(), metal_atom.get_partial_charge());
		for (int i=0; i<_no_of_binding_atoms; i++){
			reflected_vertices[i]=vertices[i].reflect_about_y_axis();
		}

		Point3D shift_amount=reflected_vertices[_vertex_order[1]].shift_to_origin();
		_oriented_metal_atom=atom(_reflected_metal_atom.get_atom_type(), _reflected_metal_atom.get_serial_no(), _reflected_metal_atom.get_atom_name(), _reflected_metal_atom.get_alt_conf_indicator(), _reflected_metal_atom.get_residue_name(), _reflected_metal_atom.get_chain_name(), _reflected_metal_atom.get_residue_no(), (_reflected_metal_atom.X()+shift_amount.X()), (_reflected_metal_atom.Y()+shift_amount.Y()), (_reflected_metal_atom.Z()+shift_amount.Z()), _reflected_metal_atom.get_occupancy(), _reflected_metal_atom.get_temp_factor(), 0.0, _reflected_metal_atom.get_atom_symbol(), _reflected_metal_atom.get_partial_charge());
		for (int i=0; i<_no_of_binding_atoms; i++){
			_oriented_vertices[i]=Point3D((reflected_vertices[i].X()+shift_amount.X()), (reflected_vertices[i].Y()+shift_amount.Y()), (reflected_vertices[i].Z()+shift_amount.Z()));
		}

		int v=_vertex_order[0];
		_vertex_order[0]=_vertex_order[1];
		_vertex_order[1]=v;

		cout<<"\n\nReflect-Shift 2nd Vertex:";
		print(_oriented_vertices);
		cout<<"\n\nReflect-Shift_2nd Metal:";
		_oriented_metal_atom.fprint(stdout);
		cout<<"\n\nVertex Order after Reflect-Shift 2nd Vertex: "<<_vertex_order[0]<<"  "<<_vertex_order[1]<<"  "<<_vertex_order[2]<<"  "<<_vertex_order[3]<<endl;
	}

	atom tetraMetalicStructure::reflect_4th_vertex_about_xy_plane(Point3D vertices[], atom metal_atom){
		Point3D reflected_vertices[_no_of_binding_atoms];
		Point3D metal_point=Point3D(metal_atom.X(), metal_atom.Y(), metal_atom.Z());
		Point3D reflected_metal_point=metal_point.reflect_about_xy_plane();
		atom _reflected_metal_atom=atom(metal_atom.get_atom_type(), metal_atom.get_serial_no(), metal_atom.get_atom_name(), metal_atom.get_alt_conf_indicator(), metal_atom.get_residue_name(), metal_atom.get_chain_name(), metal_atom.get_residue_no(), reflected_metal_point.X(), reflected_metal_point.Y(), reflected_metal_point.Z(), metal_atom.get_occupancy(), metal_atom.get_temp_factor(), 0.0, metal_atom.get_atom_symbol(), metal_atom.get_partial_charge());

		for (int i=0; i<_no_of_binding_atoms; i++){
			_reflected_vertices_4th_vertex[i]=vertices[i].reflect_about_xy_plane();
		}

		return _reflected_metal_atom;
	}

	tetraMetalicStructure tetraMetalicStructure::shift_circum_centre_to_origin(){
		tetrahedron oriented_tetra=tetrahedron(_oriented_vertices);
		oriented_tetra.compute_circum_centre_radius();
		oriented_tetra.compute_centre_of_mass();
		oriented_tetra.measure_circum_centre_to_vertex_distances();
		oriented_tetra.measure_angles_vertices_with_circum_centre();
		oriented_tetra.set_vertex_order(_vertex_order[0], _vertex_order[1], _vertex_order[2], _vertex_order[3]);

		Point3D shift_amount=oriented_tetra.get_circum_centre().shift_to_origin();
		_cc_shifted_metal_atom=atom(_oriented_metal_atom.get_atom_type(), _oriented_metal_atom.get_serial_no(), _oriented_metal_atom.get_atom_name(), _oriented_metal_atom.get_alt_conf_indicator(), _oriented_metal_atom.get_residue_name(), _oriented_metal_atom.get_chain_name(), _oriented_metal_atom.get_residue_no(), (_oriented_metal_atom.X()+shift_amount.X()), (_oriented_metal_atom.Y()+shift_amount.Y()), (_oriented_metal_atom.Z()+shift_amount.Z()), _oriented_metal_atom.get_occupancy(), _oriented_metal_atom.get_temp_factor(), 0.0, _oriented_metal_atom.get_atom_symbol(), _oriented_metal_atom.get_partial_charge());

		Point3D* vertices=oriented_tetra.get_vertices();
		for (int i=0; i<_no_of_binding_atoms; i++){
			_cc_shifted_vertices[i]=Point3D((vertices[i].X()+shift_amount.X()), (vertices[i].Y()+shift_amount.Y()), (vertices[i].Z()+shift_amount.Z()));
		}

		//Generate atoms corresponding to circumcentre shifted vertices
		for (int i=0; i<_no_of_binding_atoms; i++){
			_cc_shifted_tetra_atoms[i]=atom(_oriented_tetra_atoms[i].get_atom_type(), _oriented_tetra_atoms[i].get_serial_no(), _oriented_tetra_atoms[i].get_atom_name(), _oriented_tetra_atoms[i].get_alt_conf_indicator(), _oriented_tetra_atoms[i].get_residue_name(), _oriented_tetra_atoms[i].get_chain_name(), _oriented_tetra_atoms[i].get_residue_no(), _cc_shifted_vertices[i].X(), _cc_shifted_vertices[i].Y(), _cc_shifted_vertices[i].Z(), _oriented_tetra_atoms[i].get_occupancy(), _oriented_tetra_atoms[i].get_temp_factor(), 0.0, _oriented_tetra_atoms[i].get_atom_symbol(), _oriented_tetra_atoms[i].get_partial_charge());
		}

		cout<<"\n\nCircumcentre Shifted Vertices:";
		print(_cc_shifted_vertices);

		cout<<"\n\nCC-Shift_2nd Metal:";
		_cc_shifted_metal_atom.fprint(stdout);

		cout<<"\n\nVertex Order after CC-Shift: "<<_vertex_order[0]<<"  "<<_vertex_order[1]<<"  "<<_vertex_order[2]<<"  "<<_vertex_order[3]<<endl;
		return tetraMetalicStructure(_cc_shifted_tetra_atoms, _cc_shifted_metal_atom);
	}

	tetraMetalicStructure tetraMetalicStructure::shift_circum_centre_to_origin(int v_order[]){
		tetrahedron oriented_tetra=tetrahedron(_oriented_vertices);
		oriented_tetra.compute_circum_centre_radius();
		oriented_tetra.compute_centre_of_mass();
		oriented_tetra.measure_circum_centre_to_vertex_distances();
		oriented_tetra.measure_angles_vertices_with_circum_centre();
		oriented_tetra.set_vertex_order(v_order[0], v_order[1], v_order[2], v_order[3]);

		Point3D shift_amount=oriented_tetra.get_circum_centre().shift_to_origin();
		_cc_shifted_metal_atom=atom(_oriented_metal_atom.get_atom_type(), _oriented_metal_atom.get_serial_no(), _oriented_metal_atom.get_atom_name(), _oriented_metal_atom.get_alt_conf_indicator(), _oriented_metal_atom.get_residue_name(), _oriented_metal_atom.get_chain_name(), _oriented_metal_atom.get_residue_no(), (_oriented_metal_atom.X()+shift_amount.X()), (_oriented_metal_atom.Y()+shift_amount.Y()), (_oriented_metal_atom.Z()+shift_amount.Z()), _oriented_metal_atom.get_occupancy(), _oriented_metal_atom.get_temp_factor(), 0.0, _oriented_metal_atom.get_atom_symbol(), _oriented_metal_atom.get_partial_charge());

		Point3D* vertices=oriented_tetra.get_vertices();
		for (int i=0; i<_no_of_binding_atoms; i++){
			_cc_shifted_vertices[i]=Point3D((vertices[i].X()+shift_amount.X()), (vertices[i].Y()+shift_amount.Y()), (vertices[i].Z()+shift_amount.Z()));
		}

		//Generate atoms corresponding to circumcentre shifted vertices
		for (int i=0; i<_no_of_binding_atoms; i++){
			_cc_shifted_tetra_atoms[i]=atom(_oriented_tetra_atoms[i].get_atom_type(), _oriented_tetra_atoms[i].get_serial_no(), _oriented_tetra_atoms[i].get_atom_name(), _oriented_tetra_atoms[i].get_alt_conf_indicator(), _oriented_tetra_atoms[i].get_residue_name(), _oriented_tetra_atoms[i].get_chain_name(), _oriented_tetra_atoms[i].get_residue_no(), _cc_shifted_vertices[i].X(), _cc_shifted_vertices[i].Y(), _cc_shifted_vertices[i].Z(), _oriented_tetra_atoms[i].get_occupancy(), _oriented_tetra_atoms[i].get_temp_factor(), 0.0, _oriented_tetra_atoms[i].get_atom_symbol(), _oriented_tetra_atoms[i].get_partial_charge());
		}

		cout<<"\n\nCircumcentre Shifted Vertices:";
		print(_cc_shifted_vertices);

		cout<<"\n\nCC-Shift_2nd Metal:";
		_cc_shifted_metal_atom.fprint(stdout);

		cout<<"\n\nVertex Order after CC-Shift: "<<_vertex_order[0]<<"  "<<_vertex_order[1]<<"  "<<_vertex_order[2]<<"  "<<_vertex_order[3]<<endl;
		return tetraMetalicStructure(_cc_shifted_tetra_atoms, _cc_shifted_vertices, _cc_shifted_metal_atom, v_order[0], v_order[1], v_order[2], v_order[3]);
	}

	void tetraMetalicStructure::measure_direction_cosines(){
		for (int p=0; p<_no_of_sides; p++){
			_direction_cosines[p][0]=_tetra.measure_direction_cosines_of_sides()[p][0];
			_direction_cosines[p][1]=_tetra.measure_direction_cosines_of_sides()[p][1];
			_direction_cosines[p][2]=_tetra.measure_direction_cosines_of_sides()[p][2];
		}
	}

	void tetraMetalicStructure::write_atoms_pdb_format(string fname){
		FILE* outf=fopen(fname.c_str(), "w");
		fprintf(outf, "%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _metal_atom.get_atom_type().c_str(), _metal_atom.get_serial_no(), _metal_atom.get_atom_name().c_str(), _metal_atom.get_residue_name().c_str(), _metal_atom.get_chain_name().c_str(), _metal_atom.get_residue_no(), _metal_atom.X(), _metal_atom.Y(), _metal_atom.Z(), _metal_atom.get_occupancy(), _metal_atom.get_temp_factor(), _metal_atom.get_atom_symbol().c_str());
		for (int i=0; i<_no_of_binding_atoms; i++){
			fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _tetra_atoms[i].get_atom_type().c_str(), _tetra_atoms[i].get_serial_no(), _tetra_atoms[i].get_atom_name().c_str(), _tetra_atoms[i].get_residue_name().c_str(), _tetra_atoms[i].get_chain_name().c_str(), _tetra_atoms[i].get_residue_no(), _tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z(), _tetra_atoms[i].get_occupancy(), _tetra_atoms[i].get_temp_factor(), _tetra_atoms[i].get_atom_symbol().c_str());
		}

		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetra_atoms[0].get_serial_no(), _metal_atom.get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetra_atoms[1].get_serial_no(), _metal_atom.get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetra_atoms[2].get_serial_no(), _metal_atom.get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetra_atoms[3].get_serial_no(), _metal_atom.get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld%5ld%5ld%5ld                                                 ", _metal_atom.get_serial_no(), _tetra_atoms[0].get_serial_no(), _tetra_atoms[1].get_serial_no(), _tetra_atoms[2].get_serial_no(), _tetra_atoms[3].get_serial_no());
		fclose(outf);
	}

	void tetraMetalicStructure::write_atoms_single_tex_format(string fname){
		string out_file_name=fname+".tex";
		FILE* outf=fopen(out_file_name.c_str(), "w");
		fprintf(outf, "\\documentclass[a4paper]{article}\n");
		fprintf(outf, "\\usepackage[utf8]{inputenc}\n");
		fprintf(outf, "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}");
		fprintf(outf, "\\usepackage{chemfig}\n");
		fprintf(outf, "\\begin{document}\n");
		fprintf(outf, "\\setbondoffset{3pt}\n");
		fprintf(outf, "\\setbondstyle{line width=0.5pt,red}\n");
		fprintf(outf, "\\setnodestyle{draw,inner sep=2pt, blue}\n");

		int _no_of_amino_acids=20;
		string amino_acid_O1[_no_of_amino_acids][2]={{"ALA", "=[::45](-[::45]{O2})(-[::-90](-[::90]{{CH_3}})(-[::-45]{{NH_2}}))"},
					{"ARG", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]-[::90]-[::-90]N(-[::-45]H)-[::90](=[::45]{{NH}})(-[::-90]{{NH_2}}))(-[::-45]{{NH_2}}))"},
					{"ASN", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]-[::90](=[::45]O)(-[::-90]{{NH_2}}))(-[::-45]{{NH_2}}))"},
					{"ASP", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90](=[::90]O)(-[::-45]O^{-}))(-[::-45]{{NH_2}}))"},
					{"CYS", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]{{SH}})(-[::-45]{{NH_2}}))"},
					{"GLN", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]-[::90](=[::45]O)(-[::-90]{{NH_2}}))(-[::-45]{{NH_2}}))"},
					{"GLU", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]-[::90](=[::45]O)(-[::-90]O^{-}))(-[::-45]{{NH_2}}))"},
					{"GLY", "=[::45](-[::45]{O2})(-[::-90](-[::-45]{{NH_2}}))"},
					{"HIS", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]([::45](*5(=-{NH}-=N-))))(-[::-45]{{NH_2}}))"},
					{"ILE", "=[::45](-[::45]{O2})(-[::-90](-[::90]{{CH_3}}-[::-90]-[::90]{{CH_3}})(-[::-45]{{NH_2}}))"},
					{"LEU", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90](-[::45]{{CH_3}})(-[::-45]{{CH_3}}))(-[::-45]{{NH_2}}))"},
					{"LYS", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]-[::90]-[::-90]-[::90]{{NH_2}})(-[::-45]{{NH_2}}))"},
					{"MET", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]-[::90]S-[::-90]{{CH_3}})(-[::-45]{{NH_2}}))"},
					{"PHE", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]([::45](*6(=-=-=-))))(-[::-45]{{NH_2}}))"},
					{"PRO", "=[::45](-[::45]{O2})(-[::-90](-[::135]H)([::0](*5(-{NH}----))))"},
					{"SER", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]O^{-})(-[::-45]{{NH_2}}))"},
					{"THR", "=[::45](-[::45]{O2})(-[::-90](-[::90](-[::45]O^{-})(-[::-90]{{CH_3}}))(-[::-45]{{NH_2}}))"},
					{"TRP", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]([::45]*5(=-{HN}-(*6(-=-=-=))--)))(-[::-45]{{NH_2}}))"},
					{"TYR", "=[::45](-[::45]{O2})(-[::-90](-[::90]-[::-90]([::0]*6(-=-(-[::-45]O^{-})=-=)))(-[::-45]{{NH_2}}))"},
					{"VAL", "=[::45](-[::45]{O2})(-[::-90](-[::90](-[::45]{{CH_3}})(-[::-45]{{CH_3}}))(-[::-45]{{NH_2}}))"}};

		string amino_acid_O2[_no_of_amino_acids][2]={{"ALA", "-[::45](=[::45]{O1})(-[::-90](-[::90]{{CH_3}})(-[::-45]{{NH_2}}))"},
					{"ARG", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]-[::90]-[::-90]N(-[::-45]H)-[::90](=[::45]{{NH}})(-[::-90]{{NH_2}}))(-[::-45]{{NH_2}}))"},
					{"ASN", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]-[::90](=[::45]O)(-[::-90]{{NH_2}}))(-[::-45]{{NH_2}}))"},
					{"ASP", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90](=[::90]O)(-[::-45]O^{-}))(-[::-45]{{NH_2}}))"},
					{"CYS", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]{{SH}})(-[::-45]{{NH_2}}))"},
					{"GLN", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]-[::90](=[::45]O)(-[::-90]{{NH_2}}))(-[::-45]{{NH_2}}))"},
					{"GLU", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]-[::90](=[::45]O)(-[::-90]O^{-}))(-[::-45]{{NH_2}}))"},
					{"GLY", "-[::45](=[::45]{O1})(-[::-90](-[::-45]{{NH_2}}))"},
					{"HIS", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]([::45](*5(=-{NH}-=N-))))(-[::-45]{{NH_2}}))"},
					{"ILE", "-[::45](=[::45]{O1})(-[::-90](-[::90]{{CH_3}}-[::-90]-[::90]{{CH_3}})(-[::-45]{{NH_2}}))"},
					{"LEU", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90](-[::45]{{CH_3}})(-[::-45]{{CH_3}}))(-[::-45]{{NH_2}}))"},
					{"LYS", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]-[::90]-[::-90]-[::90]{{NH_2}})(-[::-45]{{NH_2}}))"},
					{"MET", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]-[::90]S-[::-90]{{CH_3}})(-[::-45]{{NH_2}}))"},
					{"PHE", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]([::45](*6(=-=-=-))))(-[::-45]{{NH_2}}))"},
					{"PRO", "-[::45](=[::45]{O1})(-[::-90](-[::135]H)([::0](*5(-{NH}----))))"},
					{"SER", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]O^{-})(-[::-45]{{NH_2}}))"},
					{"THR", "-[::45](=[::45]{O1})(-[::-90](-[::90](-[::45]O^{-})(-[::-90]{{CH_3}}))(-[::-45]{{NH_2}}))"},
					{"TRP", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]([::45]*5(=-{HN}-(*6(-=-=-=))--)))(-[::-45]{{NH_2}}))"},
					{"TYR", "-[::45](=[::45]{O1})(-[::-90](-[::90]-[::-90]([::0]*6(-=-(-[::-45]O^{-})=-=)))(-[::-45]{{NH_2}}))"},
					{"VAL", "-[::45](=[::45]{O1})(-[::-90](-[::90](-[::45]{{CH_3}})(-[::-45]{{CH_3}}))(-[::-45]{{NH_2}}))"}};

		string st[_no_of_binding_atoms], st1, str;
		double bond_angle[_no_of_binding_atoms];
		double bond_len;
		string adjust="";

		if (_tetra.get_circum_radius()>8.0){
			bond_len=3.5;
			adjust="adjusted";
		}else{
			bond_len=_tetra.get_circum_radius();
		}

		Point3D _centre=_tetra.get_circum_centre();

		for (int i=0; i<_no_of_binding_atoms; i++){
			_vertex_order[i]=_tetra.get_vertex_order()[i];
		}

	//cout<<"\nVERTORDER in tex: "<<_vertex_order[0]<<"  "<<_vertex_order[1]<<"  "<<_vertex_order[2]<<"  "<<_vertex_order[3]<<endl;
		str="\\center\\chemname[30pt]{\\chemfig{\\chembelow[4pt]{"+_metal_atom.get_atom_symbol()+"}{"+_metal_atom.get_chain_name()+":"+_metal_atom.get_residue_name()+"_{"+to_string(_metal_atom.get_residue_no())+"}}";

		fprintf(outf, "%s", str.c_str());
		//fprintf(outf, "%s[:%lf,%lf]\\chembelow[4pt]{%s}{%s}", str.c_str(), bond_angle[0], bond_len, str.c_str(), st1.c_str());

		st[0]="<:";
		st[1]="<:";
		st[2]="<:";
		st[3]="<";

		double ang1=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[0]], &_tetra_atoms[_vertex_order[1]]);
		bond_angle[0]=180.0+ang1;
		double ang2=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[1]], &_tetra_atoms[_vertex_order[0]]);
		bond_angle[1]=-ang2;

		if (_tetra_atoms[_vertex_order[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&_tetra_atoms[_vertex_order[1]], &_centre, &_tetra_atoms[_vertex_order[2]]);
			bond_angle[2]=ang-ang2;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&_tetra_atoms[_vertex_order[0]], &_centre, &_tetra_atoms[_vertex_order[2]]);
			bond_angle[2]=ang-ang1;
		}

		if ((_tetra_atoms[_vertex_order[3]].X()>=0.0) && (_tetra_atoms[_vertex_order[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[3]], &_tetra_atoms[_vertex_order[2]]);
			bond_angle[3]=90-ang/2.0;
		}else if ((_tetra_atoms[_vertex_order[3]].X()<0.0) && (_tetra_atoms[_vertex_order[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[3]], &_tetra_atoms[_vertex_order[2]]);
			bond_angle[3]=90+ang/2.0;
		}else if ((_tetra_atoms[_vertex_order[3]].X()<0.0) && (_tetra_atoms[_vertex_order[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[3]], &_tetra_atoms[_vertex_order[0]]);
			bond_angle[3]=180+ang/2.0;
		}else  if ((_tetra_atoms[_vertex_order[3]].X()>=0.0) && (_tetra_atoms[_vertex_order[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[1]], &_tetra_atoms[_vertex_order[3]]);
			bond_angle[3]=360-ang/2.0;
		}

		for (int p=0; p<_no_of_binding_atoms; p++){
			str="";

			if (_tetra_atoms[_vertex_order[p]].get_residue_name()=="HOH"){
				st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{H_2O}{"+st1+"}";
			}else{
				if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="O"){
					for (int x=0; x<_no_of_amino_acids; x++){
						if (_tetra_atoms[_vertex_order[p]].get_residue_name()==amino_acid_O1[x][0]){
							st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
							str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_symbol()+"1}{"+st1+"}"+amino_acid_O1[x][1];
							break;
						}
					}
				}else if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="OXT"){
					for (int x=0; x<_no_of_amino_acids; x++){
						if (_tetra_atoms[_vertex_order[p]].get_residue_name()==amino_acid_O2[x][0]){
							st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
							str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_symbol()+"2}{"+st1+"}"+amino_acid_O2[x][1];
							break;
						}
					}
				}else if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="OD1"){
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_name()+"}{"+st1+"}"+"=[::45](-[::45]{{OD2}})(-[::-45,,,,dashed])";
				}else if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="OD2"){
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_name()+"}{"+st1+"}"+"-[::45](=[::45]{{OD1}})(-[::-45,,,,dashed])";
				}else if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="OE1"){
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_name()+"}{"+st1+"}"+"=[::45](-[::45]{{OE2}})(-[::-45,,,,dashed])";
				}else if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="OE2"){
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_name()+"}{"+st1+"}"+"-[::45](=[::45]{{OE1}})(-[::-45,,,,dashed])";
				}else if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="OG1"){
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_name()+"}{"+st1+"}"+"=[::45](-[::45]{{OG2}})(-[::-45,,,,dashed])";
				}else if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="OG2"){
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_name()+"}{"+st1+"}"+"-[::45](=[::45]{{OG1}})(-[::-45,,,,dashed])";
				}else{
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{"+_tetra_atoms[_vertex_order[p]].get_atom_name()+"}{"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf]%s)", st[p].c_str(), bond_angle[p], bond_len, str.c_str());
		}

		str="}}{"+_metal_atom.get_atom_symbol()+" Binding Atoms}\n";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Tetrahedral Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Tetrahedral Site to Any Atom\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); If Bond Length Exceeds Page Width, It Is Set Equal to the Page Width\\end{flushleft}";

		fprintf(outf, "%s", str.c_str());

		if (adjust!=""){
			fprintf(outf, "\n\\vspace{2ex}Here Bond Length is Adjusted to the Page Width");
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

	void tetraMetalicStructure::write_atoms_combined_tex_format(tetraMetalicStructure & ctdrn, string fname){
		string out_file_name=fname+".tex";
		FILE* outf=fopen(out_file_name.c_str(), "w");
		fprintf(outf, "\\documentclass[a4paper]{article}\n");
		fprintf(outf, "\\usepackage[utf8]{inputenc}\n");
		fprintf(outf, "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}");
		fprintf(outf, "\\usepackage{chemfig}\n");
		fprintf(outf, "\\begin{document}\n");
		fprintf(outf, "\\setbondoffset{3pt}\n");
		string st[_no_of_binding_atoms], st1, str;
		double bond_angle_standard[_no_of_binding_atoms];
		double bond_angle_comparable[_no_of_binding_atoms];
		atom* _comp_tetra_atoms=ctdrn.get_tetra_atoms();
		atom _comp_metal_atom=ctdrn.get_metal();
		Point3D _centre=Point3D(0.0, 0.0, 0.0);
		int* _comp_vertex_order=ctdrn.get_tetrahedron().get_vertex_order();

		double bond_len_standard;
		string adjust_stand="";
		string adjust_comp="";

		if (_tetra.get_circum_radius()>=8.0){
			bond_len_standard=3.5;
			adjust_stand="adjusted";

		}else{
			bond_len_standard=_tetra.get_circum_radius();
		}

		double bond_len_comparable;

		if (ctdrn.get_tetrahedron().get_circum_radius()>=8.0){
			bond_len_comparable=3.5;
			adjust_comp="adjusted";
		}else{
			bond_len_comparable=ctdrn.get_tetrahedron().get_circum_radius();
		}

		str="\\center\\chemname[30pt]{\\chemfig{\\chembelow[4pt]{Circumcentre}{Origin}";
		fprintf(outf, "%s", str.c_str());

		st[0]="<:";
		st[1]="<:";
		st[2]="<:";
		st[3]="<";

		//Compute bond angles for Standard Tetrahedron
		double ang1=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[0]], &_tetra_atoms[_vertex_order[1]]);
		bond_angle_standard[0]=180.0+ang1;
		double ang2=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[1]], &_tetra_atoms[_vertex_order[0]]);
		bond_angle_standard[1]=-ang2;

		if (_tetra_atoms[_vertex_order[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&_tetra_atoms[_vertex_order[1]], &_centre, &_tetra_atoms[_vertex_order[2]]);
			bond_angle_standard[2]=ang-ang2;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&_tetra_atoms[_vertex_order[0]], &_centre, &_tetra_atoms[_vertex_order[2]]);
			bond_angle_standard[2]=ang-ang1;
		}

		if ((_tetra_atoms[_vertex_order[3]].X()>=0.0) && (_tetra_atoms[_vertex_order[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[3]], &_tetra_atoms[_vertex_order[2]]);
			bond_angle_standard[3]=90-ang/2.0;
		}else if ((_tetra_atoms[_vertex_order[3]].X()<0.0) && (_tetra_atoms[_vertex_order[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[3]], &_tetra_atoms[_vertex_order[2]]);
			bond_angle_standard[3]=90+ang/2.0;
		}else if ((_tetra_atoms[_vertex_order[3]].X()<0.0) && (_tetra_atoms[_vertex_order[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[3]], &_tetra_atoms[_vertex_order[0]]);
			bond_angle_standard[3]=180+ang/2.0;
		}else  if ((_tetra_atoms[_vertex_order[3]].X()>=0.0) && (_tetra_atoms[_vertex_order[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_tetra_atoms[_vertex_order[1]], &_tetra_atoms[_vertex_order[3]]);
			bond_angle_standard[3]=360-ang/2.0;
		}

		//Compute bond angles for Comparable Tetrahedron
		double ang3=Point3D::angle_deg_triangle_law(&_centre, &_comp_tetra_atoms[&_centre, _comp_vertex_order[0]], &_comp_tetra_atoms[_comp_vertex_order[1]]);
		bond_angle_comparable[0]=180.0+ang3;
		double ang4=Point3D::angle_deg_triangle_law(&_centre, &_comp_tetra_atoms[_comp_vertex_order[1]], &_comp_tetra_atoms[_comp_vertex_order[0]]);
		bond_angle_comparable[1]=-ang4;

		if (_comp_tetra_atoms[_comp_vertex_order[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&_comp_tetra_atoms[_comp_vertex_order[1]], &_centre, &_comp_tetra_atoms[_comp_vertex_order[2]]);
			bond_angle_comparable[2]=ang-ang4;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&_comp_tetra_atoms[_comp_vertex_order[0]], &_centre, &_comp_tetra_atoms[_comp_vertex_order[2]]);
			bond_angle_comparable[2]=ang-ang3;
		}

		if ((_comp_tetra_atoms[_comp_vertex_order[3]].X()>=0.0) && (_comp_tetra_atoms[_comp_vertex_order[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_comp_tetra_atoms[_comp_vertex_order[3]], &_comp_tetra_atoms[_comp_vertex_order[2]]);
			bond_angle_comparable[3]=90-ang/2.0;
		}else if ((_comp_tetra_atoms[_comp_vertex_order[3]].X()<0.0) && (_comp_tetra_atoms[_comp_vertex_order[3]].Y()>=0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_comp_tetra_atoms[_comp_vertex_order[3]], &_comp_tetra_atoms[_comp_vertex_order[2]]);
			bond_angle_comparable[3]=90+ang/2.0;
		}else if ((_comp_tetra_atoms[_comp_vertex_order[3]].X()<0.0) && (_comp_tetra_atoms[_comp_vertex_order[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_comp_tetra_atoms[_comp_vertex_order[3]], &_comp_tetra_atoms[_comp_vertex_order[0]]);
			bond_angle_comparable[3]=180+ang/2.0;
		}else  if ((_comp_tetra_atoms[_comp_vertex_order[3]].X()>=0.0) && (_comp_tetra_atoms[_comp_vertex_order[3]].Y()<0.0)){
			double ang=Point3D::angle_deg_triangle_law(&_centre, &_comp_tetra_atoms[_comp_vertex_order[1]], &_comp_tetra_atoms[_comp_vertex_order[3]]);
			bond_angle_comparable[3]=360-ang/2.0;
		}

		for (int p=0; p<_no_of_binding_atoms; p++){
			//Write Standard Tetrahedron Vertices
			str="";

			if (_tetra_atoms[_vertex_order[p]].get_residue_name()=="HOH"){
				st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="O"){
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_tetra_atoms[_vertex_order[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_tetra_atoms[_vertex_order[p]].get_atom_name()=="OXT"){
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_tetra_atoms[_vertex_order[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_tetra_atoms[_vertex_order[p]].get_chain_name()+":"+_tetra_atoms[_vertex_order[p]].get_residue_name()+"_{"+to_string(_tetra_atoms[_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_tetra_atoms[_vertex_order[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st[p].c_str(), bond_angle_standard[p], bond_len_standard, str.c_str());

			//Write Comparable Tetrahedron Vertices
			str="";

			if (_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_name()=="HOH"){
				st1=_comp_tetra_atoms[_comp_vertex_order[p]].get_chain_name()+":"+_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_name()+"_{"+to_string(_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{blue}H_2O}{\\color{blue}"+st1+"}";
			}else{
				if (_comp_tetra_atoms[_comp_vertex_order[p]].get_atom_name()=="O"){
					st1=_comp_tetra_atoms[_comp_vertex_order[p]].get_chain_name()+":"+_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_name()+"_{"+to_string(_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+_comp_tetra_atoms[_comp_vertex_order[p]].get_atom_symbol()+"1}{\\color{blue}"+st1+"}";
				}else if (_comp_tetra_atoms[_comp_vertex_order[p]].get_atom_name()=="OXT"){
					st1=_comp_tetra_atoms[_comp_vertex_order[p]].get_chain_name()+":"+_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_name()+"_{"+to_string(_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+_comp_tetra_atoms[_comp_vertex_order[p]].get_atom_symbol()+"2}{\\color{blue}"+st1+"}";
				}else{
					st1=_comp_tetra_atoms[_comp_vertex_order[p]].get_chain_name()+":"+_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_name()+"_{"+to_string(_comp_tetra_atoms[_comp_vertex_order[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{blue}"+_comp_tetra_atoms[_comp_vertex_order[p]].get_atom_name()+"}{\\color{blue}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,blue]%s)", st[p].c_str(), bond_angle_comparable[p], bond_len_comparable, str.c_str());
		}

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_comp_metal_atom.get_atom_symbol()+" Binding Atoms}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Tetrahedral Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Tetrahedral Site to Any Atom\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); If Bond Length Exceeds Page Width, It Is Set Equal to the Page Width\\end{flushleft}";

		fprintf(outf, "%s", str.c_str());

		if (adjust_stand!=""){
			fprintf(outf, "\n\\vspace{2ex}Here Bond Length for Standard Site is Adjusted to the Page Width");
		}

		if (adjust_comp!=""){
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



	void tetraMetalicStructure::compare(tetraMetalicStructure & ctdrn){
		_dist_devt_from_cc=(_tetra.get_circum_radius()-ctdrn.get_tetrahedron().get_circum_radius());

		double* deviate_dist=new double[_no_of_sides];
		_tetra.deviation_distances(ctdrn.get_tetrahedron(), deviate_dist);
		for (int p=0; p<_no_of_sides; p++){
			_devt_dist_vert_vert[p]=deviate_dist[p];
		}

		delete[] deviate_dist;

		double* deviate_angle=new double[_no_of_angles_with_circum_centre];
		_tetra.deviation_angles_with_circum_centre(ctdrn.get_tetrahedron(), deviate_angle);
		for (int p=0; p<_no_of_angles_with_circum_centre; p++){
			_devt_angle_with_cc[p]=deviate_angle[p];
		}

		delete[] deviate_angle;

		double* deviate_angle1=new double[_no_of_vertex_vertex_angles];
		_tetra.deviation_angles_vertex_vertex(ctdrn.get_tetrahedron(), deviate_angle1);
		for (int p=0; p<_no_of_vertex_vertex_angles; p++){
			_devt_angle_vertex_vertex[p]=deviate_angle1[p];
			cout<<"\nVert-Vert Angle IN COMPARE: "<<(p+1)<<": "<<_devt_angle_vertex_vertex[p];
		}

		delete[] deviate_angle1;

		compute_mse();

		int p;
		for (p=0; p<_no_of_sides; p++){
			if (abs(_devt_dist_vert_vert[p])<3.0){
				continue;
			}else{
				break;
			}
		}

		if (p==_no_of_sides){
			_nearer=true;
		}else{
			_nearer=false;
		}

		if (_dist_devt_from_cc<0.0){
				_structure_size="larger";
		}else if (_dist_devt_from_cc>0.0){
			_structure_size="smaller";
		}else{
			_structure_size="same";
		}
	}

	void tetraMetalicStructure::compute_mse(){
		_mse[0]=sqrt((_dist_devt_from_cc * _dist_devt_from_cc));

		double sum=0.0;
		for (int i=0; i<_no_of_sides; i++){
			sum+=2*(_devt_dist_vert_vert[i]/2.0) * (_devt_dist_vert_vert[i]/2.0);
		}

		_mse[1]=sqrt((sum/(2*_no_of_sides)));
		cout<<"\n MSE Distance: "<<_mse[1];

		sum=0.0;
		for (int i=0; i<_no_of_angles_with_circum_centre; i++){
			sum+=2*(_devt_angle_with_cc[i]/2.0) * (_devt_angle_with_cc[i]/2.0);
		}

		_mse[2]=sqrt((sum/(2*_no_of_angles_with_circum_centre)));
		cout<<"\n MSE Angle (with circumcentre): "<<_mse[2];

		sum=0.0;
		for (int i=0; i<_no_of_vertex_vertex_angles; i++){
			sum+=2*(_devt_angle_vertex_vertex[i]/2.0) * (_devt_angle_vertex_vertex[i]/2.0);
		}

		_mse[3]=sqrt((sum/(2*_no_of_vertex_vertex_angles)));
		cout<<"\n MSE Angle (between vertex to vertex): "<<_mse[3];
	}

	void tetraMetalicStructure::gen_report_metal_binding_sites(tetraMetalicStructure & ori_stand, tetraMetalicStructure & ori_comp, tetraMetalicStructure & trans_comp, string mt, string standard_site_type, string comparable_site_type, FILE* fp){
		fprintf(fp, "HEADER     PROGRAM NAME: StructureDeviation   VERSION: 1.1");
		string str="TITLE      Deviation Report for Computation of Structural Deviations of Two "+mt+" ION Binding Sites";
		fprintf(fp, "\n%s", str.c_str());

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Standard Tetrahedron:                   ");
		fprintf(fp, "\nREMARK     1  Molecule Type of Standard Tetrahedron: %s", standard_site_type.c_str());
		fprintf(fp, "\nREMARK     1  Circumcentre of transformed tetrahedron: (%.3lf, %.3lf, %.3lf)          ", _tetra.get_circum_centre().X(), _tetra.get_circum_centre().Y(), _tetra.get_circum_centre().Z());
		fprintf(fp, "\nREMARK     1  Circumradius of transformed tetrahedron: %.3lf", _tetra.get_circum_radius());
		double diff1=_tetra.get_circum_centre().dist(&_metal_atom);
		fprintf(fp, "\nREMARK     1  Distance between circumcentre & metal ion: %.3lf       ", diff1);

		double diff2=ori_stand.get_tetrahedron().get_centre_of_mass().dist(&_metal_atom);
		fprintf(fp, "\nREMARK     1  Distance between Centre of Mass & metal ion: %.3lf       ", diff2);

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Distances of Atoms from Circumcentre of Tetrahedron:");

		for (int i=0; i<4; i++){
			for (int j=0; j<2; j++){
				_cc_vertex_distances_stand[i][j]=_tetra.get_circum_centre_to_vertex_distances()[i][j];
			}
		}

		gen_report_circumcentre_vertex_distances(fp, _tetra_atoms, _cc_vertex_distances_stand);

		fprintf(fp, "\nREMARK     1  Distances of Atoms-Atoms:                         ");

		for (int i=0; i<6; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_stand[i][j]=_tetra.get_vertex_vertex_distances()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, _tetra_atoms, _vertex_vertex_distances_stand);

		fprintf(fp, "\nREMARK     1  Angles made by two atoms with Circumcentre of Tetrahedron:");

		for (int i=0; i<6; i++){
			for (int j=0; j<3; j++){
				_cc_vertex_angles_stand[i][j]=_tetra.get_angles_vertices_with_circum_centre()[i][j];
			}
		}

		gen_report_circumcentre_vertex_angles(fp, _tetra_atoms, _cc_vertex_angles_stand);

		fprintf(fp, "\nREMARK     1  Angles made by two atoms with another atom of Tetrahedron:");

		for (int i=0; i<12; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_stand[i][j]=_tetra.get_vertex_vertex_angles()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, _tetra_atoms, _vertex_vertex_angles_stand);

		fprintf(fp, "\nREMARK     2");
		fprintf(fp, "\nREMARK     2  Comparable Tetrahedron:");
		fprintf(fp, "\nREMARK     1  Molecule Type of Comparable Tetrahedron: %s", comparable_site_type.c_str());
		tetrahedron _circc=trans_comp.get_tetrahedron();
		Point3D _met=Point3D(trans_comp.get_metal().X(), trans_comp.get_metal().Y(), trans_comp.get_metal().Z());
		fprintf(fp, "\nREMARK     2  Circumcentre of transformed tetrahedron: (%.3lf, %.3lf, %.3lf)   ", _circc.get_circum_centre().X(), _circc.get_circum_centre().Y(), _circc.get_circum_centre().Z());
		fprintf(fp, "\nREMARK     2  Circumradius of transformed tetrahedron: %.3lf", _circc.get_circum_radius());

		if (_nearer==true){
			fprintf(fp, "\nREMARK     2  Binding Atoms of This Site are Nearer to that of Standard Site");
		}else{
			fprintf(fp, "\nREMARK     2  Binding Atoms of Comparable Site are Far Apart as Compared to Standard Site");
		}

		if (_structure_size=="larger"){
			fprintf(fp, "\nREMARK     2  Size of This Site is Larger than the Standard Site");
		}else if (_structure_size=="smaller"){
			fprintf(fp, "\nREMARK     2  Size of This Site is smaller than the Standard Site");
		}else if(_structure_size=="same"){
			fprintf(fp, "\nREMARK     2  Size of This Site is same as that of the Standard Site");
		}

		double diff3=_circc.get_circum_centre().dist(&_met);
		fprintf(fp, "\nREMARK     2  Distance between Circumcentre & Metal ion: %.3lf", diff3);

		double diff4=ori_comp.get_tetrahedron().get_centre_of_mass().dist(&_metal_atom);
		fprintf(fp, "\nREMARK     2  Distance between Centre of Mass & metal ion: %.3lf       ", diff4);

		for (int i=0; i<4; i++){
			for (int j=0; j<2; j++){
				_cc_vertex_distances_comp[i][j]=trans_comp.get_tetrahedron().get_circum_centre_to_vertex_distances()[i][j];
			}
		}

		fprintf(fp, "\nREMARK     2  Distances of Atoms from  Circumcentre of Tetrahedron:");
		gen_report_circumcentre_vertex_distances(fp, trans_comp.get_tetra_atoms(), _cc_vertex_distances_comp);

		fprintf(fp, "\nREMARK     2  Distances of Atoms from  Metal Ion:");

		for (int i=0; i<6; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_comp[i][j]=trans_comp.get_tetrahedron().get_vertex_vertex_distances()[i][j];
			}
		}

		fprintf(fp, "\nREMARK     2  Distances of Atoms-Atoms:");
		gen_report_vertex_vertex_distances(fp, trans_comp.get_tetra_atoms(), _vertex_vertex_distances_comp);

		for (int i=0; i<6; i++){
			for (int j=0; j<3; j++){
				_cc_vertex_angles_comp[i][j]=trans_comp.get_tetrahedron().get_angles_vertices_with_circum_centre()[i][j];
			}
		}

		fprintf(fp, "\nREMARK     2  Angles made by two atoms with  Circumcentre of Tetrahedron:");
		gen_report_circumcentre_vertex_angles(fp, trans_comp.get_tetra_atoms(), _cc_vertex_angles_comp);

		for (int i=0; i<12; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_comp[i][j]=trans_comp.get_tetrahedron().get_vertex_vertex_angles()[i][j];
			}
		}

		fprintf(fp, "\nREMARK     2  Angles made by two atoms with  another atom of Tetrahedron:");
		gen_report_vertex_vertex_angles(fp, trans_comp.get_tetra_atoms(), _vertex_vertex_angles_comp);

		fprintf(fp, "\nREMARK     3");
		fprintf(fp, "\nREMARK     3  Deviation Results of Comparable Tetrahedron from Standard Tetrahedron:");
		fprintf(fp, "\nREMARK     3  Distance Differences between Circumcentre & metal ion %s: %.3lf  ", mt.c_str(), (diff1-diff2));

		fprintf(fp, "\nREMARK     3  Distance Differences from Circumcentre to any Atom %s: %.3lf  ", mt.c_str(), _dist_devt_from_cc);
		fprintf(fp, "\nREMARK     3  RMSE for Distance Deviation from Circumcentre to any Atom %s: %.3lf   ", mt.c_str(), _mse[0]);

		fprintf(fp, "\nREMARK     3  Distance Differences of Atoms-Atoms: ");
		for (int i=0; i<_no_of_sides; i++){
			fprintf(fp, "\nREMARK     3  SIDE %d: %.3lf    ", (i+1), _devt_dist_vert_vert[i]);
		}

		fprintf(fp, "\nREMARK     3  RMSE for Atom-Atom Distance %s: %.3lf   ", mt.c_str(), _mse[1]);

		fprintf(fp, "\nREMARK     3  Differences of Angles Between Two Atoms with Circumcentre:");
		for (int i=0; i<_no_of_angles_with_circum_centre; i++){
			fprintf(fp, "\nREMARK     3  Angle %d: %.3lf   ", (i+1), _devt_angle_with_cc[i]);
		}

		fprintf(fp, "\nREMARK     3  RMSE for Angle Differences CC %s: %.3lf   ", mt.c_str(), _mse[2]);

		fprintf(fp, "\nREMARK     3  Differences of Angles Between Two Atoms with Another:");
		for (int i=0; i<_no_of_vertex_vertex_angles; i++){
			fprintf(fp, "\nREMARK     3  Angle %d: %.3lf   ", (i+1), _devt_angle_vertex_vertex[i]);
		}

		fprintf(fp, "\nREMARK     3  RMSE for Angle Differences AA %s: %.3lf   ", mt.c_str(), _mse[3]);
	}

	void tetraMetalicStructure::gen_pdb_model_file(tetraMetalicStructure & ori_stand, tetraMetalicStructure & ori_comp, tetraMetalicStructure & trans_comp, FILE* fp){
		string st1="MODEL        1";
		string st2="ENDMDL";
		string st3;

		int i=(8-st1.length());
		for (int p=0; p<i; p++){
			st3+=" ";
		}

		fprintf(fp, "%s%s", st1.c_str(), st3.c_str());
		fprintf(fp, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", ori_stand.get_metal().get_atom_type().c_str(), ori_stand.get_metal().get_serial_no(), ori_stand.get_metal().get_atom_name().c_str(), ori_stand.get_metal().get_residue_name().c_str(), ori_stand.get_metal().get_chain_name().c_str(), ori_stand.get_metal().get_residue_no(), ori_stand.get_metal().X(), ori_stand.get_metal().Y(), ori_stand.get_metal().Z(), ori_stand.get_metal().get_occupancy(), ori_stand.get_metal().get_temp_factor(), ori_stand.get_metal().get_atom_symbol().c_str());
		for (int i=0; i<_no_of_binding_atoms; i++){
			fprintf(fp, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", ori_stand.get_tetra_atoms()[i].get_atom_type().c_str(), ori_stand.get_tetra_atoms()[i].get_serial_no(), ori_stand.get_tetra_atoms()[i].get_atom_name().c_str(), ori_stand.get_tetra_atoms()[i].get_residue_name().c_str(), ori_stand.get_tetra_atoms()[i].get_chain_name().c_str(), ori_stand.get_tetra_atoms()[i].get_residue_no(), ori_stand.get_tetra_atoms()[i].X(), ori_stand.get_tetra_atoms()[i].Y(), ori_stand.get_tetra_atoms()[i].Z(), ori_stand.get_tetra_atoms()[i].get_occupancy(), ori_stand.get_tetra_atoms()[i].get_temp_factor(), ori_stand.get_tetra_atoms()[i].get_atom_symbol().c_str());
		}

		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", ori_stand.get_tetra_atoms()[0].get_serial_no(), ori_stand.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", ori_stand.get_tetra_atoms()[1].get_serial_no(), ori_stand.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", ori_stand.get_tetra_atoms()[2].get_serial_no(), ori_stand.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", ori_stand.get_tetra_atoms()[3].get_serial_no(), ori_stand.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld%5ld%5ld%5ld                                                 ", ori_stand.get_metal().get_serial_no(), ori_stand.get_tetra_atoms()[0].get_serial_no(), ori_stand.get_tetra_atoms()[1].get_serial_no(), ori_stand.get_tetra_atoms()[2].get_serial_no(), ori_stand.get_tetra_atoms()[3].get_serial_no());
		fprintf(fp, "\n%s", st2.c_str());

		st1="MODEL        2";
		fprintf(fp, "\n%s%s", st1.c_str(), st3.c_str());
		fprintf(fp, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", ori_comp.get_metal().get_atom_type().c_str(), ori_comp.get_metal().get_serial_no(), ori_comp.get_metal().get_atom_name().c_str(), ori_comp.get_metal().get_residue_name().c_str(), ori_comp.get_metal().get_chain_name().c_str(), ori_comp.get_metal().get_residue_no(), ori_comp.get_metal().X(), ori_comp.get_metal().Y(), ori_comp.get_metal().Z(), ori_comp.get_metal().get_occupancy(), ori_comp.get_metal().get_temp_factor(), ori_comp.get_metal().get_atom_symbol().c_str());
		for (int i=0; i<_no_of_binding_atoms; i++){
			fprintf(fp, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", ori_comp.get_tetra_atoms()[i].get_atom_type().c_str(), ori_comp.get_tetra_atoms()[i].get_serial_no(), ori_comp.get_tetra_atoms()[i].get_atom_name().c_str(), ori_comp.get_tetra_atoms()[i].get_residue_name().c_str(), ori_comp.get_tetra_atoms()[i].get_chain_name().c_str(), ori_comp.get_tetra_atoms()[i].get_residue_no(), ori_comp.get_tetra_atoms()[i].X(), ori_comp.get_tetra_atoms()[i].Y(), ori_comp.get_tetra_atoms()[i].Z(), ori_comp.get_tetra_atoms()[i].get_occupancy(), ori_comp.get_tetra_atoms()[i].get_temp_factor(), ori_comp.get_tetra_atoms()[i].get_atom_symbol().c_str());
		}

		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", ori_comp.get_tetra_atoms()[0].get_serial_no(), ori_comp.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", ori_comp.get_tetra_atoms()[1].get_serial_no(), ori_comp.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", ori_comp.get_tetra_atoms()[2].get_serial_no(), ori_comp.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", ori_comp.get_tetra_atoms()[3].get_serial_no(), ori_comp.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld%5ld%5ld%5ld                                                 ", ori_comp.get_metal().get_serial_no(), ori_comp.get_tetra_atoms()[0].get_serial_no(), ori_comp.get_tetra_atoms()[1].get_serial_no(), ori_comp.get_tetra_atoms()[2].get_serial_no(), ori_comp.get_tetra_atoms()[3].get_serial_no());
		fprintf(fp, "\n%s", st2.c_str());

		st1="MODEL        3";
		fprintf(fp, "\n%s%s", st1.c_str(), st3.c_str());
		fprintf(fp, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _metal_atom.get_atom_type().c_str(), _metal_atom.get_serial_no(), _metal_atom.get_atom_name().c_str(), _metal_atom.get_residue_name().c_str(), _metal_atom.get_chain_name().c_str(), _metal_atom.get_residue_no(), _metal_atom.X(), _metal_atom.Y(), _metal_atom.Z(), _metal_atom.get_occupancy(), _metal_atom.get_temp_factor(), _metal_atom.get_atom_symbol().c_str());
		for (int i=0; i<_no_of_binding_atoms; i++){
			fprintf(fp, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _tetra_atoms[i].get_atom_type().c_str(), _tetra_atoms[i].get_serial_no(), _tetra_atoms[i].get_atom_name().c_str(), _tetra_atoms[i].get_residue_name().c_str(), _tetra_atoms[i].get_chain_name().c_str(), _tetra_atoms[i].get_residue_no(), _tetra_atoms[i].X(), _tetra_atoms[i].Y(), _tetra_atoms[i].Z(), _tetra_atoms[i].get_occupancy(), _tetra_atoms[i].get_temp_factor(), _tetra_atoms[i].get_atom_symbol().c_str());
		}

		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", _tetra_atoms[0].get_serial_no(), _metal_atom.get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", _tetra_atoms[1].get_serial_no(), _metal_atom.get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", _tetra_atoms[2].get_serial_no(), _metal_atom.get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", _tetra_atoms[3].get_serial_no(), _metal_atom.get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld%5ld%5ld%5ld                                                 ", _metal_atom.get_serial_no(), _tetra_atoms[0].get_serial_no(), _tetra_atoms[1].get_serial_no(), _tetra_atoms[2].get_serial_no(), _tetra_atoms[3].get_serial_no());
		fprintf(fp, "\n%s", st2.c_str());

		st1="MODEL        4";
		fprintf(fp, "\n%s%s", st1.c_str(), st3.c_str());
		fprintf(fp, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", trans_comp.get_metal().get_atom_type().c_str(), trans_comp.get_metal().get_serial_no(), trans_comp.get_metal().get_atom_name().c_str(), trans_comp.get_metal().get_residue_name().c_str(), trans_comp.get_metal().get_chain_name().c_str(), trans_comp.get_metal().get_residue_no(), trans_comp.get_metal().X(), trans_comp.get_metal().Y(), trans_comp.get_metal().Z(), trans_comp.get_metal().get_occupancy(), trans_comp.get_metal().get_temp_factor(), trans_comp.get_metal().get_atom_symbol().c_str());
		for (int i=0; i<_no_of_binding_atoms; i++){
			fprintf(fp, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", trans_comp.get_tetra_atoms()[i].get_atom_type().c_str(), trans_comp.get_tetra_atoms()[i].get_serial_no(), trans_comp.get_tetra_atoms()[i].get_atom_name().c_str(), trans_comp.get_tetra_atoms()[i].get_residue_name().c_str(), trans_comp.get_tetra_atoms()[i].get_chain_name().c_str(), trans_comp.get_tetra_atoms()[i].get_residue_no(), trans_comp.get_tetra_atoms()[i].X(), trans_comp.get_tetra_atoms()[i].Y(), trans_comp.get_tetra_atoms()[i].Z(), trans_comp.get_tetra_atoms()[i].get_occupancy(), trans_comp.get_tetra_atoms()[i].get_temp_factor(), trans_comp.get_tetra_atoms()[i].get_atom_symbol().c_str());
		}

		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", trans_comp.get_tetra_atoms()[0].get_serial_no(), trans_comp.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", trans_comp.get_tetra_atoms()[1].get_serial_no(), trans_comp.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", trans_comp.get_tetra_atoms()[2].get_serial_no(), trans_comp.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld                                                                ", trans_comp.get_tetra_atoms()[3].get_serial_no(), trans_comp.get_metal().get_serial_no());
		fprintf(fp, "\nCONECT%5ld%5ld%5ld%5ld%5ld                                                 ", trans_comp.get_metal().get_serial_no(), trans_comp.get_tetra_atoms()[0].get_serial_no(), trans_comp.get_tetra_atoms()[1].get_serial_no(), trans_comp.get_tetra_atoms()[2].get_serial_no(), trans_comp.get_tetra_atoms()[3].get_serial_no());
		fprintf(fp, "\n%s", st2.c_str());

		fprintf(fp, "\nEND");
	}

	void tetraMetalicStructure::gen_report_circumcentre_vertex_distances(FILE* fp, atom atm[], double cc_vert_distances[][2]){
		for (int i=0; i<_no_of_binding_atoms; i++){
			fprintf(fp, "\nDISTCR %ld %8.3lf", atm[(int) cc_vert_distances[i][0]].get_serial_no(), cc_vert_distances[i][1]);
		}
	}

	void tetraMetalicStructure::gen_report_vertex_vertex_distances(FILE* fp, atom atm[], double vert_vert_distances[][3]){
		for (int i=0; i<_no_of_sides; i++){
			fprintf(fp, "\nDISTRR %ld %ld %8.3lf", atm[(int) vert_vert_distances[i][0]].get_serial_no(), atm[(int) vert_vert_distances[i][1]].get_serial_no(), vert_vert_distances[i][2]);
		}
	}

	void tetraMetalicStructure::gen_report_circumcentre_vertex_angles(FILE* fp, atom atm[], double circumcentre_vert_angles[][3]){
		for (int i=0; i<_no_of_angles_with_circum_centre; i++){
			fprintf(fp, "\nANGLEC %ld %ld %8.3lf", atm[(int) circumcentre_vert_angles[i][0]].get_serial_no(), atm[(int) circumcentre_vert_angles[i][1]].get_serial_no(), circumcentre_vert_angles[i][2]);
		}
	}

	void tetraMetalicStructure::gen_report_vertex_vertex_angles(FILE* fp, atom atm[], double vert_vert_angles[][4]){
		for (int i=0; i<_no_of_vertex_vertex_angles; i++){
			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", atm[(int) vert_vert_angles[i][0]].get_serial_no(), atm[(int) vert_vert_angles[i][1]].get_serial_no(), atm[(int) vert_vert_angles[i][2]].get_serial_no(), vert_vert_angles[i][3]);
		}
	}
