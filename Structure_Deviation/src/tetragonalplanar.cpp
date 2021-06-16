#include "tetragonalplanar.h"

void tetraplMetalicStructure::print(atom mt, atom* non_mt){
	mt.fprint(stdout);

	for (int i=0; i<_no_of_binding_atoms; i++){
		non_mt[i].fprint(stdout);
	}
}

void tetraplMetalicStructure::print(Point3D* vert){
	for (int i=0; i<_no_of_binding_atoms; i++){
		cout<<"\nVertex "<<i<<": "<<vert[i].X()<<"  "<<vert[i].Y()<<"  "<<vert[i].Z();
	}
}

tetraplMetalicStructure::tetraplMetalicStructure(){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetragonal Planar No. of Vertices should be 4!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=4;
	assert(_no_of_sides==4 && "For Tetragonal Planar No. of Sides should be 4!!");

	_no_of_vertex_vertex_angles=4;
	assert(_no_of_vertex_vertex_angles==4 && "For Tetragonal Planar No. of Angles of Two Vertices with Another Vertex should be 4!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Tetragonal Planar No. of Features to Compare should be 2!!");

	int no_of_adjacent_pts_to_any_pt=2;
	int no_of_non_adjacent_pts_to_any_pt=_no_of_binding_atoms-no_of_adjacent_pts_to_any_pt;
	_no_of_triangles=no_of_adjacent_pts_to_any_pt+no_of_non_adjacent_pts_to_any_pt-1;
}

tetraplMetalicStructure::tetraplMetalicStructure(atom* new_atm, atom metal_atom){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetragonal Planar No. of Vertices should be 4!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=4;
	assert(_no_of_sides==4 && "For Tetragonal Planar No. of Sides should be 4!!");

	_no_of_vertex_vertex_angles=4;
	assert(_no_of_vertex_vertex_angles==4 && "For Tetragonal Planar No. of Angles of Two Vertices with Another Vertex should be 4!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Tetragonal Planar No. of Features to Compare should be 2!!");

	int no_of_adjacent_pts_to_any_pt=2;
	int no_of_non_adjacent_pts_to_any_pt=_no_of_binding_atoms-no_of_adjacent_pts_to_any_pt;
	_no_of_triangles=no_of_adjacent_pts_to_any_pt+no_of_non_adjacent_pts_to_any_pt-1;

	_metal_atom=metal_atom;

	for (int i=0; i<_no_of_binding_atoms; i++){
		_tetrapl_atoms[i]=new_atm[i];
	}

	for (int i=0; i<_no_of_binding_atoms; i++){
		_tetrapl_vertices[i]=Point3D(_tetrapl_atoms[i].X(), _tetrapl_atoms[i].Y(), _tetrapl_atoms[i].Z());
	}

	_tetrapl=tetragonalplanar(_tetrapl_vertices);
	print(_metal_atom, _tetrapl_atoms);

	for (int i=0; i<_no_of_binding_atoms; i++){
		_vertex_order[i]=_tetrapl.get_vertex_order()[i];
	}
}

tetraplMetalicStructure::tetraplMetalicStructure(vector<direct_ligand> new_atm, atom metal_atom){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetragonal Planar No. of Vertices should be 4!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=4;
	assert(_no_of_sides==4 && "For Tetragonal Planar No. of Sides should be 4!!");

	_no_of_vertex_vertex_angles=4;
	assert(_no_of_vertex_vertex_angles==4 && "For Tetragonal Planar No. of Angles of Two Vertices with Another Vertex should be 4!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Tetragonal Planar No. of Features to Compare should be 2!!");

	int no_of_adjacent_pts_to_any_pt=2;
	int no_of_non_adjacent_pts_to_any_pt=_no_of_binding_atoms-no_of_adjacent_pts_to_any_pt;
	_no_of_triangles=no_of_adjacent_pts_to_any_pt+no_of_non_adjacent_pts_to_any_pt-1;

	_metal_atom=metal_atom;

	for (int i=0; i<new_atm.size(); i++){
		_tetrapl_atoms[i]=new_atm.at(i).get_atom_direct_contact();
	}

	for (int i=0; i<_no_of_binding_atoms; i++){
		_tetrapl_vertices[i]=Point3D(_tetrapl_atoms[i].X(), _tetrapl_atoms[i].Y(), _tetrapl_atoms[i].Z());
	}

	_tetrapl=tetragonalplanar(_tetrapl_vertices);
	print(_metal_atom, _tetrapl_atoms);

	for (int i=0; i<_no_of_binding_atoms; i++){
		_vertex_order[i]=_tetrapl.get_vertex_order()[i];
	}
}

tetraplMetalicStructure::tetraplMetalicStructure(triangle lt, atom atom_tri1[], triangle rt, atom atom_tri2[], atom met){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetragonal Planar No. of Vertices should be 4!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=4;
	assert(_no_of_sides==4 && "For Tetragonal Planar No. of Sides should be 4!!");

	_no_of_vertex_vertex_angles=4;
	assert(_no_of_vertex_vertex_angles==4 && "For Tetragonal Planar No. of Angles of Two Vertices with Another Vertex should be 4!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Tetragonal Planar No. of Features to Compare should be 2!!");

	int no_of_adjacent_pts_to_any_pt=2;
	int no_of_non_adjacent_pts_to_any_pt=_no_of_binding_atoms-no_of_adjacent_pts_to_any_pt;
	_no_of_triangles=no_of_adjacent_pts_to_any_pt+no_of_non_adjacent_pts_to_any_pt-1;

	_metal_atom=met;
	_left_tri=lt;
	_right_tri=rt;

	for (int i=0; i<3; i++){
		_cc_shifted_vertices_left[i]=_left_tri.get_vertices()[i];
	}

	for (int i=0; i<3; i++){
		_cc_shifted_vertices_right[i]=_right_tri.get_vertices()[i];
	}

	for (int i=0; i<3; i++){
		_cc_shifted_left_triangle_atoms[i]=atom_tri1[i];
	}

	for (int i=0; i<3; i++){
		_cc_shifted_right_triangle_atoms[i]=atom_tri2[i];
	}

	for (int i=0; i<3; i++){
		_left_vertex_order[i]=_left_tri.get_vertex_order()[i];
	}

	for (int i=0; i<3; i++){
		_right_vertex_order[i]=_right_tri.get_vertex_order()[i];
	}
}

tetraplMetalicStructure::tetraplMetalicStructure(triplMetalicStructure tps_left, triplMetalicStructure tps_right, atom metal){
	_no_of_binding_atoms=4;
	assert(_no_of_binding_atoms==4 && "For Tetragonal Planar No. of Vertices should be 4!!");

	_no_of_coords=3;
	assert(_no_of_coords==3 && "No. of Coordinates must be 3!!");

	_no_of_sides=4;
	assert(_no_of_sides==4 && "For Tetragonal Planar No. of Sides should be 4!!");

	_no_of_vertex_vertex_angles=4;
	assert(_no_of_vertex_vertex_angles==4 && "For Tetragonal Planar No. of Angles of Two Vertices with Another Vertex should be 4!!");

	_no_of_features_to_compare=2;
	assert(_no_of_features_to_compare==2 && "For Tetragonal Planar No. of Features to Compare should be 2!!");

	int no_of_adjacent_pts_to_any_pt=2;
	int no_of_non_adjacent_pts_to_any_pt=_no_of_binding_atoms-no_of_adjacent_pts_to_any_pt;
	_no_of_triangles=no_of_adjacent_pts_to_any_pt+no_of_non_adjacent_pts_to_any_pt-1;

	_metal_atom=metal;
	_transformed_tri_left=tps_left;
	_left_tri=tps_left.get_triangle();

	for (int i=0; i<3; i++){
		_cc_shifted_vertices_left[i]=tps_left.get_cc_shifted_tripl_vertices()[i];
	}

	for (int i=0; i<3; i++){
		_cc_shifted_left_triangle_atoms[i]=tps_left.get_tripl_atoms()[i];
	}

	for (int i=0; i<3; i++){
		_left_vertex_order[i]=tps_left.get_triangle().get_vertex_order()[i];
	}

	_transformed_tri_right=tps_right;
	_right_tri=tps_right.get_triangle();

	for (int i=0; i<3; i++){
		_cc_shifted_vertices_right[i]=tps_right.get_cc_shifted_tripl_vertices()[i];
	}

	for (int i=0; i<3; i++){
		_cc_shifted_right_triangle_atoms[i]=tps_right.get_tripl_atoms()[i];
		cout<<"\nPrint 1: "<<_cc_shifted_right_triangle_atoms[i].get_atom_name();
		cout<<"\nPrint 2: "<<_cc_shifted_right_triangle_atoms[i].get_atom_symbol();
		cout<<"\nPrint 3: "<<_cc_shifted_right_triangle_atoms[i].get_residue_name();
		cout<<"\nPrint 4: "<<_cc_shifted_right_triangle_atoms[i].get_residue_no();
	}

	for (int i=0; i<3; i++){
		_right_vertex_order[i]=tps_right.get_triangle().get_vertex_order()[i];
	}

	cout<<"\nVERTEX ORDER After Orientation Left: "<<_left_vertex_order[0]<<"  "<<_left_vertex_order[1]<<"  "<<_left_vertex_order[2];
	cout<<"\nVERTEX ORDER After Orientation Right: "<<_right_vertex_order[0]<<"  "<<_right_vertex_order[1]<<"  "<<_right_vertex_order[2];
}

tetraplMetalicStructure tetraplMetalicStructure::apply_transformations(){
	cout<<"\n Original Vertices";
	print(_tetrapl_vertices);

	for (int i=0; i<_no_of_binding_atoms; i++){
		_vertex_order[i]=_tetrapl.get_vertex_order()[i];
	}

	find_largest_side_and_update_vertex_order();
	tetraplMetalicStructure tms=make_transformation();

	return tms;
}

atom tetraplMetalicStructure::get_metal(){
	return _metal_atom;
}

Point3D* tetraplMetalicStructure::get_tetrapl_vertices(){
		return _tetrapl_vertices;
	}

atom* tetraplMetalicStructure::get_tetrapl_atoms(){
	return _tetrapl_atoms;
}

tetragonalplanar tetraplMetalicStructure::get_tetragonalplanar(){
	return _tetrapl;
}

int* tetraplMetalicStructure::get_vertex_order(){
	return _vertex_order;
}

triangle tetraplMetalicStructure::get_left_triangle(){
	return _left_tri;
}

triangle tetraplMetalicStructure::get_right_triangle(){
	return _right_tri;
}

atom* tetraplMetalicStructure::get_cc_shifted_left_triangle_atoms(){
	return _cc_shifted_left_triangle_atoms;
}

atom* tetraplMetalicStructure::get_cc_shifted_right_triangle_atoms(){
	return _cc_shifted_right_triangle_atoms;
}

int tetraplMetalicStructure::get_no_of_binding_atoms(){
	return _no_of_binding_atoms;
}

double* tetraplMetalicStructure::get_side_deviations(){
	return _devt_dist_vert_vert;
}

double* tetraplMetalicStructure::get_angle_deviations_vertex_vertex(){
	return _devt_angle_vertex_vertex;
}

double* tetraplMetalicStructure::get_mse(){
	return _mse;
}

string tetraplMetalicStructure::get_structure_size(){
	return _structure_size;
}

void tetraplMetalicStructure::comparison(tetraplMetalicStructure & comp, string out_comparable_name){
	double dist_stand_tetra_planar[4][3];
	for (int i=0; i<4; i++){
		for (int j=0; j<3; j++){
			dist_stand_tetra_planar[i][j]=_tetrapl.get_vertex_vertex_distances()[i][j];
		}
	}

	Point3D* comp_vertices=comp.get_tetrapl_vertices();
	measure_direction_cosines();

	int k=0;
	double proj[4];

	double prx=	_direction_cosines[0][0]*(comp.get_tetrapl_vertices()[1].X()-comp.get_tetrapl_vertices()[0].X());
	double pry= _direction_cosines[0][1]*(comp.get_tetrapl_vertices()[1].Y()-comp.get_tetrapl_vertices()[0].Y());
	double prz=	_direction_cosines[0][2]*(comp.get_tetrapl_vertices()[1].Z()-comp.get_tetrapl_vertices()[0].Z());

	proj[0]=(prx+pry+prz);

	prx= _direction_cosines[1][0]*(comp.get_tetrapl_vertices()[2].X()-comp.get_tetrapl_vertices()[1].X());
	pry= _direction_cosines[1][1]*(comp.get_tetrapl_vertices()[2].Y()-comp.get_tetrapl_vertices()[1].Y());
	prz= _direction_cosines[1][2]*(comp.get_tetrapl_vertices()[2].Z()-comp.get_tetrapl_vertices()[1].Z());

	proj[1]=(prx+pry+prz);

	prx= _direction_cosines[2][0]*(comp.get_tetrapl_vertices()[3].X()-comp.get_tetrapl_vertices()[2].X());
	pry= _direction_cosines[2][1]*(comp.get_tetrapl_vertices()[3].Y()-comp.get_tetrapl_vertices()[2].Y());
	prz= _direction_cosines[2][2]*(comp.get_tetrapl_vertices()[3].Z()-comp.get_tetrapl_vertices()[2].Z());

	proj[2]=(prx+pry+prz);

	prx= _direction_cosines[3][0]*(comp.get_tetrapl_vertices()[0].X()-comp.get_tetrapl_vertices()[3].X());
	pry= _direction_cosines[3][1]*(comp.get_tetrapl_vertices()[0].Y()-comp.get_tetrapl_vertices()[3].Y());
	prz= _direction_cosines[3][2]*(comp.get_tetrapl_vertices()[0].Z()-comp.get_tetrapl_vertices()[3].Z());

	proj[3]=(prx+pry+prz);

	double dist_dev_vert_vert[4];
	FILE* outf=fopen(out_comparable_name.c_str(), "w");

	cout<<"\n\nComparison Dist:\n";
	for (int i=0; i<_no_of_sides; i++){
		double diff=dist_stand_tetra_planar[i][2]-proj[i];
		dist_dev_vert_vert[i]=diff;
		cout<<"Side "<<(i+1)<<": "<<dist_dev_vert_vert[i]<<endl;
	}

	double sum=0.0;
	for (int p=0; p<4; p++){
		sum+=(dist_dev_vert_vert[p] * dist_dev_vert_vert[p]);
	}

	double mse_dist=(sum/4.0);
	cout<<"\n MSE Distance: "<<mse_dist;
	fprintf(outf, "\nOriginal MSE Dist: %8.3lf", mse_dist);

	cout<<"\n\nComparison Angle:\n";
	double angle_stand_tetra_planar[4][4];
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			angle_stand_tetra_planar[i][j]=_tetrapl.get_vertex_vertex_angles()[i][j];
		}
	}

	double angle_comp_tetra_planar[4][4];
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			angle_comp_tetra_planar[i][j]=comp.get_tetragonalplanar().get_vertex_vertex_angles()[i][j];
		}
	}

	double angle_dev_vert_vert[4];

	for (int p=0; p<_no_of_vertex_vertex_angles; p++){
		angle_dev_vert_vert[p]=angle_stand_tetra_planar[p][3]-angle_comp_tetra_planar[p][3];
		cout<<"\nAngle "<<(p+1)<<": "<<angle_dev_vert_vert[p];
	}

	sum=0.0;
	for (int p=0; p<4; p++){
		sum+=(angle_dev_vert_vert[p] * angle_dev_vert_vert[p]);
	}

	double mse_angle=(sum/4.0);
	cout<<"\n MSE Angle: "<<mse_angle;
	fprintf(outf, "\nMSE Angle: %8.3lf", mse_angle);

	fclose(outf);
}

void tetraplMetalicStructure:: find_largest_side_and_update_vertex_order(){
	double dist[4][3];
	for (int i=0; i<4; i++){
		for (int j=0; j<3; j++){
			dist[i][j]=_tetrapl.get_vertex_vertex_distances()[i][j];
		}
	}

	int* _initial_vertex_order=new int[4];

	for (int i=0; i<4; i++){
		_initial_vertex_order[i]=_tetrapl.get_vertex_order()[i];
	}

	if ((dist[0][2]>=dist[1][2]) && (dist[0][2]>=dist[2][2]) && (dist[0][2]>=dist[3][2])){
		_vertex_order[0]=_initial_vertex_order[0];
		_vertex_order[1]=_initial_vertex_order[1];

		if (dist[1][2]>=dist[3][2]){
			_vertex_order[2]=_initial_vertex_order[2];
			_vertex_order[3]=_initial_vertex_order[3];
		}else{
			_vertex_order[0]=_initial_vertex_order[1];
			_vertex_order[1]=_initial_vertex_order[0];
			_vertex_order[2]=_initial_vertex_order[3];
			_vertex_order[3]=_initial_vertex_order[2];
		}
	}else if ((dist[1][2]>=dist[0][2]) && (dist[1][2]>=dist[2][2]) && (dist[1][2]>=dist[3][2])){
		_vertex_order[0]=_initial_vertex_order[1];
		_vertex_order[1]=_initial_vertex_order[2];

		if (dist[2][2]>=dist[0][2]){
			_vertex_order[2]=_initial_vertex_order[3];
			_vertex_order[3]=_initial_vertex_order[0];
		}else{
			_vertex_order[0]=_initial_vertex_order[2];
			_vertex_order[1]=_initial_vertex_order[1];
			_vertex_order[2]=_initial_vertex_order[0];
			_vertex_order[3]=_initial_vertex_order[3];
		}
	}else if ((dist[2][2]>=dist[0][2]) && (dist[2][2]>=dist[1][2]) && (dist[2][2]>=dist[3][2])){
		_vertex_order[0]=_initial_vertex_order[2];
		_vertex_order[1]=_initial_vertex_order[3];

		if (dist[3][2]>=dist[1][2]){
			_vertex_order[2]=_initial_vertex_order[0];
			_vertex_order[3]=_initial_vertex_order[1];
		}else{
			_vertex_order[0]=_initial_vertex_order[3];
			_vertex_order[1]=_initial_vertex_order[2];
			_vertex_order[2]=_initial_vertex_order[1];
			_vertex_order[3]=_initial_vertex_order[0];
		}
	}else if ((dist[3][2]>=dist[0][2]) && (dist[3][2]>=dist[1][2]) && (dist[3][2]>=dist[2][2])){
		_vertex_order[0]=_initial_vertex_order[3];
		_vertex_order[1]=_initial_vertex_order[0];

		if (dist[0][2]>=dist[2][2]){
			_vertex_order[2]=_initial_vertex_order[1];
			_vertex_order[3]=_initial_vertex_order[2];
		}else{
			_vertex_order[0]=_initial_vertex_order[0];
			_vertex_order[1]=_initial_vertex_order[3];
			_vertex_order[2]=_initial_vertex_order[2];
			_vertex_order[3]=_initial_vertex_order[1];
		}
	}
}

	tetraplMetalicStructure tetraplMetalicStructure::make_transformation(){
		Point3D* left_triangle_vertices=new Point3D[3];
		atom left_triangle_atoms[3];
		Point3D* right_triangle_vertices=new Point3D[3];
		atom right_triangle_atoms[3];

		left_triangle_vertices[0]=_tetrapl_vertices[_vertex_order[0]];
		left_triangle_vertices[1]=_tetrapl_vertices[_vertex_order[1]];
		left_triangle_vertices[2]=_tetrapl_vertices[_vertex_order[3]];
		left_triangle_atoms[0]=_tetrapl_atoms[_vertex_order[0]];
		left_triangle_atoms[1]=_tetrapl_atoms[_vertex_order[1]];
		left_triangle_atoms[2]=_tetrapl_atoms[_vertex_order[3]];

		triplMetalicStructure leftTstruct=triplMetalicStructure(left_triangle_atoms, left_triangle_vertices, _metal_atom, 0, 1, 2);
		triplMetalicStructure  left_transformed=leftTstruct.transform_triangle();

		right_triangle_vertices[0]=_tetrapl_vertices[_vertex_order[3]];
		right_triangle_vertices[1]=_tetrapl_vertices[_vertex_order[1]];
		right_triangle_vertices[2]=_tetrapl_vertices[_vertex_order[2]];
		right_triangle_atoms[0]=_tetrapl_atoms[_vertex_order[3]];
		right_triangle_atoms[1]=_tetrapl_atoms[_vertex_order[1]];
		right_triangle_atoms[2]=_tetrapl_atoms[_vertex_order[2]];

		triplMetalicStructure rightTstruct=triplMetalicStructure(right_triangle_atoms, right_triangle_vertices, _metal_atom, 0, 1, 2);
		triplMetalicStructure  right_transformed=rightTstruct.transform_triangle();
		tetraplMetalicStructure tps=tetraplMetalicStructure(left_transformed, right_transformed, _metal_atom);

		return tps;
	}

	void tetraplMetalicStructure::measure_direction_cosines(){
		double x=_tetrapl_vertices[1].X() - _tetrapl_vertices[0].X();
		double y=_tetrapl_vertices[1].Y() - _tetrapl_vertices[0].Y();
		double z=_tetrapl_vertices[1].Z() - _tetrapl_vertices[0].Z();
		double r=sqrt(x*x+y*y+z*z);

		double l=x/r;
		double m=y/r;
		double n=z/r;

		_direction_cosines[0][0]=l;
		_direction_cosines[0][1]=m;
		_direction_cosines[0][2]=n;

		x=_tetrapl_vertices[2].X() - _tetrapl_vertices[1].X();
		y=_tetrapl_vertices[2].Y() - _tetrapl_vertices[1].Y();
		z=_tetrapl_vertices[2].Z() - _tetrapl_vertices[1].Z();
		r=sqrt(x*x+y*y+z*z);

		l=x/r;
		m=y/r;
		n=z/r;

		_direction_cosines[1][0]=l;
		_direction_cosines[1][1]=m;
		_direction_cosines[1][2]=n;

		x=_tetrapl_vertices[3].X() - _tetrapl_vertices[2].X();
		y=_tetrapl_vertices[3].Y() - _tetrapl_vertices[2].Y();
		z=_tetrapl_vertices[3].Z() - _tetrapl_vertices[2].Z();
		r=sqrt(x*x+y*y+z*z);

		l=x/r;
		m=y/r;
		n=z/r;

		_direction_cosines[2][0]=l;
		_direction_cosines[2][1]=m;
		_direction_cosines[2][2]=n;

		x=_tetrapl_vertices[0].X() - _tetrapl_vertices[3].X();
		y=_tetrapl_vertices[0].Y() - _tetrapl_vertices[3].Y();
		z=_tetrapl_vertices[0].Z() - _tetrapl_vertices[3].Z();
		r=sqrt(x*x+y*y+z*z);

		l=x/r;
		m=y/r;
		n=z/r;

		_direction_cosines[3][0]=l;
		_direction_cosines[3][1]=m;
		_direction_cosines[3][2]=n;
	}

	void tetraplMetalicStructure::write_atoms_pdb_format(string fname){
		for (int i=0; i<_no_of_binding_atoms; i++){
			_vertex_order[i]=_tetrapl.get_vertex_order()[i];
		}

		FILE* outf=fopen(fname.c_str(), "w");
		fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _metal_atom.get_atom_type().c_str(), _metal_atom.get_serial_no(), _metal_atom.get_atom_name().c_str(), _metal_atom.get_residue_name().c_str(), _metal_atom.get_chain_name().c_str(), _metal_atom.get_residue_no(), _metal_atom.X(), _metal_atom.Y(), _metal_atom.Z(), _metal_atom.get_occupancy(), _metal_atom.get_temp_factor(), _metal_atom.get_atom_symbol().c_str());
		for (int i=0; i<_no_of_binding_atoms; i++){
			fprintf(outf, "\n%-6s%5ld %4s %3s %s%4ld    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s  ", _tetrapl_atoms[_vertex_order[i]].get_atom_type().c_str(), _tetrapl_atoms[_vertex_order[i]].get_serial_no(), _tetrapl_atoms[_vertex_order[i]].get_atom_name().c_str(), _tetrapl_atoms[_vertex_order[i]].get_residue_name().c_str(), _tetrapl_atoms[_vertex_order[i]].get_chain_name().c_str(), _tetrapl_atoms[_vertex_order[i]].get_residue_no(), _tetrapl_atoms[_vertex_order[i]].X(), _tetrapl_atoms[_vertex_order[i]].Y(), _tetrapl_atoms[_vertex_order[i]].Z(), _tetrapl_atoms[_vertex_order[i]].get_occupancy(), _tetrapl_atoms[_vertex_order[i]].get_temp_factor(), _tetrapl_atoms[_vertex_order[i]].get_atom_symbol().c_str());
		}

		for (int i=0; i<_no_of_binding_atoms-1; i++){
			fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapl_atoms[_vertex_order[i]].get_serial_no(), _tetrapl_atoms[_vertex_order[i+1]].get_serial_no());
			fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapl_atoms[_vertex_order[i+1]].get_serial_no(), _tetrapl_atoms[_vertex_order[i]].get_serial_no());
		}

		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapl_atoms[_vertex_order[_no_of_binding_atoms-1]].get_serial_no(), _tetrapl_atoms[_vertex_order[0]].get_serial_no());
		fprintf(outf, "\nCONECT%5ld%5ld                                                                ", _tetrapl_atoms[_vertex_order[0]].get_serial_no(), _tetrapl_atoms[_vertex_order[_no_of_binding_atoms-1]].get_serial_no());
		fclose(outf);
	}

	void tetraplMetalicStructure::write_atoms_combined_tex_format_left(tetraplMetalicStructure & comp, string fname){
		string out_file_name=fname+".tex";
		FILE* outf=fopen(out_file_name.c_str(), "w");
		fprintf(outf, "\\documentclass[a4paper]{article}\n");
		fprintf(outf, "\\usepackage[utf8]{inputenc}\n");
		fprintf(outf, "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}");
		fprintf(outf, "\\usepackage{chemfig}\n");
		fprintf(outf, "\\begin{document}\n");
		fprintf(outf, "\\setbondoffset{3pt}\n");
		string st_left[3], st1, str;
		double bond_angle_standard_left[3];
		double bond_angle_comparable_left[3];
		atom* _comp_tetrapl_atoms_tri_left=comp.get_cc_shifted_left_triangle_atoms();
		//atom* _comp_tetrapl_atoms_tri_right=comp.get_cc_shifted_right_triangle_atoms();

		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		triangle left_comp=comp.get_left_triangle();

		Point3D* left_vertices_standard=_left_tri.get_vertices();
		Point3D* left_vertices_comparable=left_comp.get_vertices();

		int* stand_vertex_order_left=_left_tri.get_vertex_order();
		int* comp_vertex_order_left=left_comp.get_vertex_order();

		double bond_len_standard_left=_left_tri.get_circum_radius();
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

		st_left[0]="<:";
		st_left[1]="<:";
		st_left[2]="<:";

		//Compute bond angles for Standard Left Triangle
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

		//Compute bond angles for Comparable Left Triangle
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

		for (int p=0; p<3; p++){
			//Write Standard Triangle Left Vertices
			str="";

			if (_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_name()=="HOH"){
				st1=_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_atom_name()=="O"){
					st1=_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_atom_name()=="OXT"){
					st1=_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_chain_name()+":"+_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_name()+"_{"+to_string(_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_cc_shifted_left_triangle_atoms[stand_vertex_order_left[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_left[p].c_str(), bond_angle_standard_left[p], bond_len_standard_left, str.c_str());

			//Write Comparable Triangle left Vertices
			str="";
			atom* comp_left_tetra_atoms=comp.get_cc_shifted_left_triangle_atoms();
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

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms Left}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Trigonal Bipyramidal Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Trigonal Bipyramidal Site to Any Atom\\par\n";
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

	void tetraplMetalicStructure::write_atoms_combined_tex_format_right(tetraplMetalicStructure & comp, string fname){
		string out_file_name=fname+".tex";
		FILE* outf=fopen(out_file_name.c_str(), "w");
		fprintf(outf, "\\documentclass[a4paper]{article}\n");
		fprintf(outf, "\\usepackage[utf8]{inputenc}\n");
		fprintf(outf, "\\usepackage[top=0.5in, bottom=0.5in, left=0.5in, right=0.5in]{geometry}");
		fprintf(outf, "\\usepackage{chemfig}\n");
		fprintf(outf, "\\begin{document}\n");
		fprintf(outf, "\\setbondoffset{3pt}\n");
		string st_right[3], st1, str;
		double bond_angle_standard_right[3];
		double bond_angle_comparable_right[3];
		atom* _comp_tetrapl_atoms_tri_right=comp.get_cc_shifted_right_triangle_atoms();

		Point3D _centre=Point3D(0.0, 0.0, 0.0);

		triangle right_comp=comp.get_right_triangle();

		Point3D* right_vertices_standard=_right_tri.get_vertices();
		Point3D* right_vertices_comparable=right_comp.get_vertices();

		int* stand_vertex_order_right=_right_tri.get_vertex_order();
		int* comp_vertex_order_right=right_comp.get_vertex_order();

		double bond_len_standard_right=_right_tri.get_circum_radius();
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

		st_right[0]="<:";
		st_right[1]="<:";
		st_right[2]="<:";

		//Compute bond angles for Standard Right Triangle
		double ang1=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_standard[stand_vertex_order_right[0]], &right_vertices_standard[stand_vertex_order_right[1]]);
		bond_angle_standard_right[0]=180.0+ang1;
		double ang2=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_standard[stand_vertex_order_right[1]], &right_vertices_standard[stand_vertex_order_right[0]]);
		bond_angle_standard_right[1]=-ang2;

		if (right_vertices_standard[stand_vertex_order_right[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&right_vertices_standard[stand_vertex_order_right[1]], &_centre, &right_vertices_standard[stand_vertex_order_right[2]]);
			bond_angle_standard_right[2]=ang-ang2;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&right_vertices_standard[stand_vertex_order_right[0]], &_centre, &right_vertices_standard[stand_vertex_order_right[2]]);
			bond_angle_standard_right[2]=ang-ang1;
		}

		//Compute bond angles for Comparable Right Triangle
		double ang3=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_comparable[comp_vertex_order_right[0]], &right_vertices_comparable[comp_vertex_order_right[1]]);
		bond_angle_comparable_right[0]=180.0+ang3;
		double ang4=Point3D::angle_deg_triangle_law(&_centre, &right_vertices_comparable[comp_vertex_order_right[1]], &right_vertices_comparable[comp_vertex_order_right[0]]);
		bond_angle_comparable_right[1]=-ang4;

		if (right_vertices_comparable[comp_vertex_order_right[2]].X()>=0.0){
			double ang=Point3D::angle_deg_triangle_law(&right_vertices_comparable[comp_vertex_order_right[1]], &_centre, &right_vertices_comparable[comp_vertex_order_right[2]]);
			bond_angle_comparable_right[2]=ang-ang4;
		}else{
			double ang=Point3D::angle_deg_triangle_law(&right_vertices_comparable[comp_vertex_order_right[0]], &_centre, &right_vertices_comparable[comp_vertex_order_right[2]]);
			bond_angle_comparable_right[2]=ang-ang3;
		}

		for (int p=0; p<3; p++){
			//Write Standard Triangle Right Vertices
			str="";

			if (_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_name()=="HOH"){
				st1=_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
				str="\\chembelow[4 pt]{\\color{red}H_2O}{\\color{red}"+st1+"}";
			}else{
				if (_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_atom_name()=="O"){
					st1=_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_atom_symbol()+"1}{\\color{red}"+st1+"}";
				}else if (_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_atom_name()=="OXT"){
					st1=_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_atom_symbol()+"2}{\\color{red}"+st1+"}";
				}else{
					st1=_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_chain_name()+":"+_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_name()+"_{"+to_string(_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_residue_no())+"}";
					str="\\chembelow[4 pt]{\\color{red}"+_cc_shifted_right_triangle_atoms[stand_vertex_order_right[p]].get_atom_name()+"}{\\color{red}"+st1+"}";
				}
			}

			fprintf(outf, "(%s[::%lf,%lf,,,red]%s)", st_right[p].c_str(), bond_angle_standard_right[p], bond_len_standard_right, str.c_str());

			//Write Comparable Triangle Right Vertices
			str="";
			atom* comp_right_tetra_atoms=comp.get_cc_shifted_right_triangle_atoms();
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

		str="}}{Combined "+_metal_atom.get_atom_symbol()+"-"+_metal_atom.get_atom_symbol()+" Binding Atoms Right}";

		fprintf(outf, "%s", str.c_str());

		str="\\vspace{2ex}\n\\center{Legend}\n\\begin{flushleft}\\tikz[baseline=0.6ex]\\fill[red] (0,0) rectangle (0.5cm,3ex); Standard Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\fill[blue] (0,0) rectangle (0.5cm,3ex); Comparable Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Angle = Angle Between Two Atoms with  Circumcentre of the Trigonal Bipyramidal Site\\par\n";
		str+="\\tikz[baseline=0.6ex]\\draw (0,0) rectangle (0.5cm,3ex); Bond Length = Distance from Circumcentre of the Trigonal Bipyramidal Site to Any Atom\\par\n";
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

void tetraplMetalicStructure::compare_with_standard(tetraplMetalicStructure & comp){
	double *dev_dist1, *dev_dist2;
	dev_dist1=new double[3];
	dev_dist2=new double[3];

	_left_tri.deviation_distances(comp.get_left_triangle(), _left_tri.measure_direction_cosines_of_sides(), dev_dist1);
	_right_tri.deviation_distances(comp.get_right_triangle(), _right_tri.measure_direction_cosines_of_sides(), dev_dist2);

	_devt_dist_vert_vert[0]=dev_dist1[0];
	_devt_dist_vert_vert[1]=dev_dist2[1];
	_devt_dist_vert_vert[2]=dev_dist2[2];
	_devt_dist_vert_vert[3]=dev_dist1[2];

	delete[] dev_dist1;
	delete[] dev_dist2;

	double *dev_angle1, *dev_angle2;
	dev_angle1=new double[3];
	dev_angle2=new double[3];

	_left_tri.deviation_angles_vertex_vertex(comp.get_left_triangle(), dev_angle1);
	_right_tri.deviation_angles_vertex_vertex(comp.get_right_triangle(), dev_angle2);

	_devt_angle_vertex_vertex[0]=dev_angle1[0];
	_devt_angle_vertex_vertex[1]=dev_angle1[1]+dev_angle2[1];
	_devt_angle_vertex_vertex[2]=dev_angle2[2];
	_devt_angle_vertex_vertex[3]=dev_angle1[2]+dev_angle2[0];

	delete[] dev_angle1;
	delete[] dev_angle2;

	compute_mse();

	_dist_devt_from_cc_left=_left_tri.get_circum_radius()-comp.get_left_triangle().get_circum_radius();
	_dist_devt_from_cc_right=_right_tri.get_circum_radius()-comp.get_right_triangle().get_circum_radius();

	if (_dist_devt_from_cc_left<0.0 && _dist_devt_from_cc_right<0.0){
			_structure_size="smaller";
	}else if (_dist_devt_from_cc_left>0.0 || _dist_devt_from_cc_right>0.0){
		_structure_size="larger";
	}else if (_dist_devt_from_cc_left==0.0 && _dist_devt_from_cc_right==0.0){
		_structure_size="same";
	}else{
		_structure_size="irregular";
	}
}

void tetraplMetalicStructure::compute_mse(){
	double sum=0.0;
	for (int i=0; i<_no_of_sides; i++){
		sum+=2*(_devt_dist_vert_vert[i]/2.0) * (_devt_dist_vert_vert[i]/2.0);
		cout<<"\nDevt Dist: "<<_devt_dist_vert_vert[i];
	}

	_mse[0]=sqrt((sum/(2*_no_of_sides)));
	cout<<"\n MSE Distance: "<<_mse[0];

	sum=0.0;
	for (int i=0; i<_no_of_vertex_vertex_angles; i++){
		sum+=2*(_devt_angle_vertex_vertex[i]/2.0) * (_devt_angle_vertex_vertex[i]/2.0);
		cout<<"\nDevt Angle: "<<_devt_angle_vertex_vertex[i];
	}

	_mse[1]=sqrt((sum/(2*_no_of_vertex_vertex_angles)));
	cout<<"\n MSE Angle (between vertex to vertex): "<<_mse[1];
}

	void tetraplMetalicStructure::gen_report_metal_binding_sites(tetraplMetalicStructure & ori_stand, tetraplMetalicStructure & ori_comp, tetraplMetalicStructure & trans_comp, string mt, string standard_site_type, string comparable_site_type, FILE* fp){
		fprintf(fp, "HEADER     PROGRAM NAME: StructureDeviation VERSION: 1.1");
		string str="TITLE      Deviation Report for Computation of Structural Deviations of Two Tetragonal Planar "+mt+" ION Binding Sites";
		fprintf(fp, "\n%s", str.c_str());

		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Standard Tetragonal Planar:                   ");
		fprintf(fp, "\nREMARK     1  Molecule Type of Standard Tetragonal Planar: %s", standard_site_type.c_str());
		fprintf(fp, "\nREMARK     1");

		fprintf(fp, "\nREMARK     1  Distances of Atoms-Atoms:                         ");
		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_stand[i][j]=ori_stand.get_tetragonalplanar().get_vertex_vertex_distances()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_stand.get_tetrapl_atoms(), _vertex_vertex_distances_stand);

		fprintf(fp, "\nREMARK     1  Angles made by two atoms with another atom of Tetragonal Planar:");
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_stand[i][j]=ori_stand.get_tetragonalplanar().get_vertex_vertex_angles()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_stand.get_tetrapl_atoms(), _vertex_vertex_angles_stand);

		fprintf(fp, "\nREMARK     2");
		fprintf(fp, "\nREMARK     2  Comparable Triangle:");
		fprintf(fp, "\nREMARK     1  Molecule Type of Comparable Tetragonal Planar: %s", comparable_site_type.c_str());

		if (_structure_size=="larger"){
			fprintf(fp, "\nREMARK     2  Size of This Site is Larger than the Standard Site");
		}else if (_structure_size=="smaller"){
			fprintf(fp, "\nREMARK     2  Size of This Site is smaller than the Standard Site");
		}else if (_structure_size=="same"){
			fprintf(fp, "\nREMARK     2  Size of This Site is same as that of the Standard Site");
		}else if(_structure_size=="irregular"){
			fprintf(fp, "\nREMARK     2  Size of This Site is irregular in comparison to the Standard Site");
		}

		fprintf(fp, "\nREMARK     2  Inter-atomic Distances and Angles:");
		fprintf(fp, "\nREMARK     2  Distances of Atoms-Atoms:");
		for (int i=0; i<4; i++){
			for (int j=0; j<3; j++){
				_vertex_vertex_distances_comp[i][j]=ori_comp.get_tetragonalplanar().get_vertex_vertex_distances()[i][j];
			}
		}

		gen_report_vertex_vertex_distances(fp, ori_comp.get_tetrapl_atoms(), _vertex_vertex_distances_comp);

		fprintf(fp, "\nREMARK     2  Angles made by two atoms with  another atom of Tetragonal Planar:");
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				_vertex_vertex_angles_comp[i][j]=ori_comp.get_tetragonalplanar().get_vertex_vertex_angles()[i][j];
			}
		}

		gen_report_vertex_vertex_angles(fp, ori_comp.get_tetrapl_atoms(), _vertex_vertex_angles_comp);

		fprintf(fp, "\nREMARK     3");
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

	void tetraplMetalicStructure::gen_report_vertex_vertex_distances(FILE* fp, atom atm[], double vert_vert_distances[][3]){
		for (int i=0; i<_no_of_sides; i++){
			fprintf(fp, "\nDISTRR %ld %ld %8.3lf", atm[(int) vert_vert_distances[i][0]].get_serial_no(), atm[(int) vert_vert_distances[i][1]].get_serial_no(), vert_vert_distances[i][2]);
		}
	}

	void tetraplMetalicStructure::gen_report_vertex_vertex_angles(FILE* fp, atom atm[], double vert_vert_angles[][4]){
		for (int i=0; i<_no_of_vertex_vertex_angles; i++){
			fprintf(fp, "\nANGLEV %ld %ld %ld %8.3lf", atm[(int) vert_vert_angles[i][0]].get_serial_no(), atm[(int) vert_vert_angles[i][1]].get_serial_no(), atm[(int) vert_vert_angles[i][2]].get_serial_no(), vert_vert_angles[i][3]);
		}
	}
