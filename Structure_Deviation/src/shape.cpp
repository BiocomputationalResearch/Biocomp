#include "coordinationsite.h"

void metal_binding_site::determine_shape(){	//determines shape of binding site
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());

	if (_inner_coord_no==1){
		//cout<<"Start Point";
		_shape_of_inner_sphere="point";
		_no_of_base_atoms=1;
	}else if (_inner_coord_no==2){
		atom atm1=_inner_coord_sphere.at(0).get_atom_direct_contact();
		Point3D p1=Point3D(atm1.X(), atm1.Y(), atm1.Z());
		atom atm2=_inner_coord_sphere.at(1).get_atom_direct_contact();
		Point3D p2=Point3D(atm2.X(), atm2.Y(), atm2.Z());
		double angle=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);

		if (angle==180.0){
			_shape_of_inner_sphere="linear";
		}else if ((angle>0) && (angle<180.0)){
			_shape_of_inner_sphere="bent";
		}

		_no_of_base_atoms=2;
	}else if (_inner_coord_no==3){
		bool is_t_shaped=false;
		atom atmp=_inner_coord_sphere.at(0).get_atom_direct_contact();
		atom atmq=_inner_coord_sphere.at(1).get_atom_direct_contact();
		atom atmr=_inner_coord_sphere.at(2).get_atom_direct_contact();
		Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
		Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
		Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
		double angle1=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);
		double angle2=Point3D::angle_deg_triangle_law(&p2, &pm, &p3);
		double angle3=Point3D::angle_deg_triangle_law(&p3, &pm, &p1);

		if (angle1>175.0 && angle1<185.0){
			if ((angle2==(angle1/2.0)) && (angle2==angle3)){
				is_t_shaped=true;
			}
		}else if (angle2>175.0 && angle2<185.0){
			if ((angle3==(angle2/2.0)) && (angle3==angle1)){
				is_t_shaped=true;
			}
		}else if (angle3>170.0 && angle3<190.0){
			if ((angle1==(angle3/2.0)) && (angle1==angle2)){
				is_t_shaped=true;
			}
		}

		if (is_t_shaped==true){
			_shape_of_inner_sphere="t_shaped";
		}else{
			vector<Point3D> points{p1, p2, p3, pm};
			bool is_coplanar=check_for_coplanarity(points);
			if (is_coplanar==true){
				_shape_of_inner_sphere="trigonal_planar";
			}else{
				_shape_of_inner_sphere="trigonal_pyramidal";
			}
		}

		_no_of_base_atoms=3;
	}else{
		vector<int> _planar_atom_nos;
		int no_of_non_planar;
		_planar_atom_nos=count_planar_atoms();

		if (_planar_atom_nos.size()!=0){
			no_of_non_planar=_inner_coord_no-_planar_atom_nos.size();
		}else{
			no_of_non_planar=_inner_coord_no;
		}

		int* pt=new int[no_of_non_planar];

		if (_planar_atom_nos.size()!=0){
			_no_of_base_atoms=_planar_atom_nos.size();
			find_non_planar_atom_indices(_planar_atom_nos, pt);
			compute_mean_angle_planar_atoms_with_metal(_planar_atom_nos);

			if (no_of_non_planar>0){
				compute_mean_angle_non_planar_atoms_with_metal(pt, no_of_non_planar);
				compute_mean_angle_planar_and_non_planar_atoms_with_metal(pt, no_of_non_planar, _planar_atom_nos);
			}
		}else{
			_no_of_base_atoms=3;

			for (int x=0; x<no_of_non_planar; x++){
				pt[x]=x;
			}

			compute_mean_angle_non_planar_atoms_with_metal(pt, no_of_non_planar);
		}

		if (_planar_atom_nos.size()==_inner_coord_no){
			vector<Point3D> points;
			for (int x=0; x<_planar_atom_nos.size(); x++){
				atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(x)).get_atom_direct_contact();
				Point3D ptt=Point3D(atmp.X(), atmp.Y(), atmp.Z());
				points.push_back(ptt);
			}

			points.push_back(pm);

			bool is_coplanar=check_for_coplanarity(points);

			if (_inner_coord_no==4){
				bool square_planar=is_square(points.at(0), points.at(1), points.at(2), points.at(3));

				if (is_coplanar==true){
					if (square_planar==true){
						_shape_of_inner_sphere="square_planar_metal_on_plane";
					}else{
						_shape_of_inner_sphere="tetragonal_planar_metal_on_plane";
					}
				}else{
					if (square_planar==true){
						_shape_of_inner_sphere="square_planar_metal_outside";
					}else{
						_shape_of_inner_sphere="tetragonal_planar_metal_outside";
					}
				}
			}/*else if (_inner_coord_no==5){
				if (is_coplanar==true){
					_shape_of_inner_sphere="pentagonal_planar_metal_on_plane";
				}else{
					_shape_of_inner_sphere="pentagonal_planar_metal_outside";
				}
			}else if (_inner_coord_no==6){
				if (is_coplanar==true){
					_shape_of_inner_sphere="hexagonal_planar_metal_on_plane";
				}else{
					_shape_of_inner_sphere="hexagonal_planar_metal_outside";
				}
			}else if (_inner_coord_no==7){
				if (is_coplanar==true){
					_shape_of_inner_sphere="heptagonal_planar_metal_on_plane";
				}else{
					_shape_of_inner_sphere="heptagonal_planar_metal_outside";
				}
			}else if (_inner_coord_no==8){
				if (is_coplanar==true){
					_shape_of_inner_sphere="octagonal_planar_metal_on_plane";
				}else{
					_shape_of_inner_sphere="octagonal_planar_metal_outside";
				}
			}else if (_inner_coord_no==9){
				if (is_coplanar==true){
					_shape_of_inner_sphere="nonagonal_planar_metal_on_plane";
				}else{
					_shape_of_inner_sphere="nonagonal_planar_metal_outside";
				}
			}*/else{
				/*if (is_coplanar==true){
					_shape_of_inner_sphere="more_than_9_atoms_planar_metal_on_plane";
				}else{
					_shape_of_inner_sphere="more_than_9_atoms_metal_outside";
				}*/

				if (is_coplanar==true){
					_shape_of_inner_sphere="more_than_4_atoms_planar_metal_on_plane";
				}else{
					_shape_of_inner_sphere="more_than_4_atoms_metal_outside";
				}
			}
		}else if ((_inner_coord_no-_planar_atom_nos.size())==1){
			if (_inner_coord_no==5){
				_shape_of_inner_sphere="tetragonal_pyramidal";
			}/*else if (_inner_coord_no==6){
				_shape_of_inner_sphere="pentagonal_pyramidal";
			}else if (_inner_coord_no==7){
				_shape_of_inner_sphere="hexagonal_pyramidal";
			}else if (_inner_coord_no==8){
				_shape_of_inner_sphere="heptagonal_pyramidal";
			}else if (_inner_coord_no==9){
				_shape_of_inner_sphere="octagonal_pyramidal";
			}*/else{
				//_shape_of_inner_sphere="more_than_9_atoms_pyramidal_base_to_pt_diff_1";
				_shape_of_inner_sphere="more_than_5_atoms_pyramidal_base_to_pt_diff_1";
			}
		}else if (_inner_coord_no==4){
			bool is_perfect_seesaw=test_for_perfect_seesaw();
			bool is_seesaw=test_for_seesaw();
			if (is_perfect_seesaw==true){
				_shape_of_inner_sphere="perfect_seesaw";
			}else if (is_seesaw==true){
				_shape_of_inner_sphere="seesaw";
			}else{
				atom atmp=_inner_coord_sphere.at(0).get_atom_direct_contact();
				atom atmq=_inner_coord_sphere.at(1).get_atom_direct_contact();
				atom atmr=_inner_coord_sphere.at(2).get_atom_direct_contact();
				atom atms=_inner_coord_sphere.at(3).get_atom_direct_contact();
				Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
				Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
				Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
				Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
				double angle1=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);
				double angle2=Point3D::angle_deg_triangle_law(&p1, &pm, &p3);
				double angle3=Point3D::angle_deg_triangle_law(&p1, &pm, &p4);
				double angle4=Point3D::angle_deg_triangle_law(&p2, &pm, &p3);
				double angle5=Point3D::angle_deg_triangle_law(&p2, &pm, &p4);
				double angle6=Point3D::angle_deg_triangle_law(&p3, &pm, &p4);

				if ((angle1==109.4712) && (angle2==109.4712) && (angle3==109.4712) && (angle4==109.4712) && (angle5==109.4712) && (angle6==109.4712)){
					_shape_of_inner_sphere="regular_tetrahedral";
				}else{
					_shape_of_inner_sphere="tetrahedral";
				}
			}
		}else if (_inner_coord_no==5){
			if (_planar_atom_nos.size()==0){
				bool found_shape=false;

				for (int x=0; x<_inner_coord_no-2; x++){
					atom atmp=_inner_coord_sphere.at(x).get_atom_direct_contact();
					Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());

					for (int y=(x+1); y<_inner_coord_no-1; y++){
						atom atmq=_inner_coord_sphere.at(y).get_atom_direct_contact();
						Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());

						for (int z=(y+1); z<_inner_coord_no; z++){
							atom atmr=_inner_coord_sphere.at(z).get_atom_direct_contact();
							Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
							vector<Point3D> points{p1, p2, p3, pm};
							bool is_coplanar=check_for_coplanarity(points);
							if (is_coplanar==true){
								found_shape=true;
								_planar_atoms.push_back(_inner_coord_sphere.at(x).get_atom_direct_contact().get_serial_no());
								_planar_atoms.push_back(_inner_coord_sphere.at(y).get_atom_direct_contact().get_serial_no());
								_planar_atoms.push_back(_inner_coord_sphere.at(z).get_atom_direct_contact().get_serial_no());

								_planar_atom_nos.push_back(x);
								_planar_atom_nos.push_back(y);
								_planar_atom_nos.push_back(z);
								no_of_non_planar=_inner_coord_no-_planar_atom_nos.size();
								find_non_planar_atom_indices(_planar_atom_nos, pt);
								compute_mean_angle_planar_atoms_with_metal(_planar_atom_nos);
								compute_mean_angle_non_planar_atoms_with_metal(pt, no_of_non_planar);
								compute_mean_angle_planar_and_non_planar_atoms_with_metal(pt, no_of_non_planar, _planar_atom_nos);

								break;
							}
						}

						if (found_shape==true){
							break;
						}
					}

					if (found_shape==true){
						break;
					}
				}

				if (found_shape==true){
					bool same_side=is_same_side(pt, _planar_atom_nos);	//is non planar atoms are on same side of the plane

					if (same_side==true){
						_shape_of_inner_sphere="hemispheric_coord_5";
					}else{
						_shape_of_inner_sphere="trigonal_bipyramidal";
					}
				}else{
					_shape_of_inner_sphere="holospheric_coord_5";
				}
			}else{
				_shape_of_inner_sphere="holospheric_coord_5";
			}
		}else if (_inner_coord_no==6){
			if (_planar_atom_nos.size()!=0){
				if (_planar_atom_nos.size()==4){
					bool same_side=is_same_side(pt, _planar_atom_nos);	//is non planar atoms are on same side of the plane

					if (same_side==true){
						bool tri_pris=check_for_trigonal_prismatic(pt, _planar_atom_nos);

						if (tri_pris==true){
							_shape_of_inner_sphere="trigonal_prismatic";
						}else{
							_shape_of_inner_sphere="hemispheric_planar_4_coord_6";
						}

					}else{
						vector<Point3D> points;
						for (int x=0; x<_planar_atom_nos.size(); x++){
							atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(x)).get_atom_direct_contact();
							Point3D ptt=Point3D(atmp.X(), atmp.Y(), atmp.Z());
							points.push_back(ptt);
						}

						points.push_back(pm);

						bool is_coplanar=check_for_coplanarity(points);

						if (is_coplanar==true){
							bool octa=check_for_octahedral(pt, _planar_atom_nos);
							if (octa==true){
								_shape_of_inner_sphere="octahedral";
							}else{
								_shape_of_inner_sphere="tetragonal_bipyramidal";
							}
						}else{
							_shape_of_inner_sphere="holospheric_planar_4_coord_6";
						}
					}
				}
			}else{
				_shape_of_inner_sphere="holospheric_coord_6";
			}
		}/*else if (_inner_coord_no==7){
			if (_planar_atom_nos.size()!=0){
				bool same_side=is_same_side(pt, _planar_atom_nos);

				if (_planar_atom_nos.size()==4){
					if (same_side==true){
						_shape_of_inner_sphere="hemispheric_planar_4_coord_7";
					}else{
						bool capped_tri_pris=false;
						int* temp=new int[2];
						for (int x=0; x<no_of_non_planar-1; x++){
							temp[0]=pt[x];
							for (int y=(x+1); y<no_of_non_planar; y++){
								temp[1]=pt[y];
								bool tri_pris=check_for_trigonal_prismatic(temp, _planar_atom_nos);
								if (tri_pris==true){
									for (int z=0; z<no_of_non_planar; z++){
										if (z!=x && z!=y){
											capped_tri_pris=check_for_capped_trigonal_prismatic(z, _planar_atom_nos);
											if (capped_tri_pris==true){
												_shape_of_inner_sphere="capped_trigonal_prismatic";
												break;
											}
										}
									}
								}

								if (capped_tri_pris==true){
									break;
								}
							}

							if (capped_tri_pris==true){
								break;
							}
						}

						bool capped_octa=false;

						if (capped_tri_pris==false){
							vector<Point3D> points;
							for (int x=0; x<_planar_atom_nos.size(); x++){
								atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(x)).get_atom_direct_contact();
								Point3D ptt=Point3D(atmp.X(), atmp.Y(), atmp.Z());
								points.push_back(ptt);
							}

							points.push_back(pm);

							bool is_coplanar=check_for_coplanarity(points);

							if (is_coplanar==true){
								for (int x=0; x<no_of_non_planar-1; x++){
									temp[0]=pt[x];
									for (int y=(x+1); y<no_of_non_planar; y++){
										temp[1]=pt[y];
										bool octa=check_for_octahedral(temp, _planar_atom_nos);
										if (octa==true){
											for (int z=0; z<no_of_non_planar; z++){
												if (z!=x && z!=y){
													capped_octa=check_for_capped_octahedral(z, temp, _planar_atom_nos);
													if (capped_octa==true){
														_shape_of_inner_sphere="capped_octahedral";
														break;
													}
												}
											}
										}

										if (capped_octa==true){
											break;
										}
									}

									if (capped_octa==true){
										break;
									}
								}

								if (capped_octa==false){
									_shape_of_inner_sphere="holospheric_planar_4_coord_7";
								}
							}else{
								_shape_of_inner_sphere="holospheric_planar_4_coord_7";
							}
						}

						delete[] temp;
					}
				}else if (_planar_atom_nos.size()==5){
					if (same_side==true){
						_shape_of_inner_sphere="hemispheric_planar_5_coord_7";
					}else{
						vector<Point3D> points;
						for (int x=0; x<_planar_atom_nos.size(); x++){
							atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(x)).get_atom_direct_contact();
							Point3D ptt=Point3D(atmp.X(), atmp.Y(), atmp.Z());
							points.push_back(ptt);
						}

						points.push_back(pm);

						bool is_coplanar=check_for_coplanarity(points);

						if (is_coplanar==true){
							_shape_of_inner_sphere="pentagonal_bipyramidal";
						}else{
							_shape_of_inner_sphere="holospheric_planar_5_coord_7";
						}
					}
				}
			}else{
				_shape_of_inner_sphere="holospheric_coord_7";
			}
		}else if (_inner_coord_no==8){
			if (_planar_atom_nos.size()!=0){
				if (_planar_atom_nos.size()==4){
					bool bicapped_pris_found=false;
					bool capped_tri_pris=false;
					int* temp=new int[2];
					vector<int> non_pl;

					for (int x=0; x<no_of_non_planar-1; x++){
						temp[0]=pt[x];
						for (int y=(x+1); y<no_of_non_planar; y++){
							temp[1]=pt[y];
							bool capped_tri_pris=check_for_trigonal_prismatic(temp, _planar_atom_nos);
							if (capped_tri_pris==true){
								for (int z=0; z<no_of_non_planar; z++){
									if (z!=x && z!=y){
										non_pl.push_back(z);
									}
								}

								bool tetra_face_bicap=check_position_for_tetragonal_face_bicapped(non_pl, _planar_atom_nos, temp);
								if (tetra_face_bicap==true){
									_shape_of_inner_sphere="tetragonal_face_bicapped_trigonal_prismatic";
									bicapped_pris_found=true;
									break;
								}else{
									bool tri_face_bicap=check_position_for_triangular_face_bicapped(non_pl, _planar_atom_nos, temp);
									if (tri_face_bicap==true){
										_shape_of_inner_sphere="trigonal_face_bicapped_trigonal_prismatic";
										bicapped_pris_found=true;
										break;
									}
								}
							}
						}

						if (bicapped_pris_found==true){
							break;
						}
					}

					if (bicapped_pris_found==false){
						bool same_side=is_same_side(pt, _planar_atom_nos);

						if (same_side==true){
							bool cube=is_cube(pt, _planar_atom_nos);
							if (cube==true){
								_shape_of_inner_sphere="cubic";
							}else{
								bool sq_prism=is_square_prism(pt, _planar_atom_nos);
								if (sq_prism==true){
									_shape_of_inner_sphere="square_prismatic";
								}else{
									bool sq_ant_pris=check_position_for_square_antiprism(pt, _planar_atom_nos);
									if (sq_ant_pris==true){
										_shape_of_inner_sphere="square_anti_prismatic";
									}else{
										bool dod=check_for_dodecahedron(pt, _planar_atom_nos);
										if (dod==true){
											_shape_of_inner_sphere="dodecahedral";
										}else{
											_shape_of_inner_sphere="hemispheric_planar_4_coord_8";
										}
									}
								}
							}
						}else{
							_shape_of_inner_sphere="holospheric_planar_4_coord_8";
						}
					}

					delete[] temp;
				}else if (_planar_atom_nos.size()==5){
					bool same_side=is_same_side(pt, _planar_atom_nos);

					if (same_side==true){
						_shape_of_inner_sphere="hemispheric_planar_5_coord_8";
					}else{
						_shape_of_inner_sphere="holospheric_planar_5_coord_8";
					}
				}else if (_planar_atom_nos.size()==6){
					bool same_side=is_same_side(pt, _planar_atom_nos);

					if (same_side==true){
						_shape_of_inner_sphere="hemispheric_planar_6_coord_8";
					}else{
						vector<Point3D> points;
						for (int x=0; x<_planar_atom_nos.size(); x++){
							atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(x)).get_atom_direct_contact();
							Point3D ptt=Point3D(atmp.X(), atmp.Y(), atmp.Z());
							points.push_back(ptt);
						}

						points.push_back(pm);

						bool is_coplanar=check_for_coplanarity(points);

						if (is_coplanar==true){
							_shape_of_inner_sphere="hexagonal_bipyramidal";
						}else{
							_shape_of_inner_sphere="holospheric_planar_6_coord_8";
						}
					}
				}
			}else{
				_shape_of_inner_sphere="holospheric_coord_8";
			}
		}else if (_inner_coord_no==9){
			if (_planar_atom_nos.size()!=0){
				if (_planar_atom_nos.size()==4){
					bool tricapped_pris_found=false;
					bool capped_tri_pris=false;
					int* temp=new int[2];
					vector<int> non_pl;

					for (int x=0; x<no_of_non_planar-1; x++){
						temp[0]=pt[x];
						for (int y=(x+1); y<no_of_non_planar; y++){
							temp[1]=pt[y];
							bool capped_tri_pris=check_for_trigonal_prismatic(temp, _planar_atom_nos);
							if (capped_tri_pris==true){
								for (int z=0; z<no_of_non_planar; z++){
									if (z!=x && z!=y){
										non_pl.push_back(z);
									}
								}

								bool tetra_face_bicap=check_position_for_tricapped(non_pl, _planar_atom_nos, temp);
								if (tetra_face_bicap==true){
									_shape_of_inner_sphere="tetragonal_face_tricapped_trigonal_prismatic";
									tricapped_pris_found=true;
									break;
								}
							}
						}

						if (tricapped_pris_found==true){
							break;
						}
					}

					if (tricapped_pris_found==false){
						bool same_side=is_same_side(pt, _planar_atom_nos);

						if (same_side==true){
							bool capped_sq_ant_pris=check_for_capped_square_anti_prism(pt, _planar_atom_nos, no_of_non_planar);

							if (capped_sq_ant_pris==true){
								_shape_of_inner_sphere="capped_square_antiprism";
							}else{
								_shape_of_inner_sphere="hemispheric_planar_4_coord_9";
							}
						}else{
							_shape_of_inner_sphere="holospheric_planar_4_coord_9";
						}
					}

					delete[] temp;
				}else if (_planar_atom_nos.size()==5){
					bool same_side=is_same_side(pt, _planar_atom_nos);

					if (same_side==true){
						_shape_of_inner_sphere="hemispheric_planar_5_coord_9";
					}else{
						_shape_of_inner_sphere="holospheric_planar_5_coord_9";
					}
				}else if (_planar_atom_nos.size()==6){
					bool same_side=is_same_side(pt, _planar_atom_nos);

					if (same_side==true){
						_shape_of_inner_sphere="hemispheric_planar_6_coord_9";
					}else{
						_shape_of_inner_sphere="holospheric_planar_6_coord_9";
					}
				}else if (_planar_atom_nos.size()==7){
					bool same_side=is_same_side(pt, _planar_atom_nos);

					if (same_side==true){
						_shape_of_inner_sphere="hemispheric_planar_7_coord_9";
					}else{
						_shape_of_inner_sphere="holospheric_planar_7_coord_9";
					}
				}
			}else{
				_shape_of_inner_sphere="holospheric_coord_9";
			}
		}*/else if (_inner_coord_no>6){
			//_shape_of_inner_sphere="holospheric_coord_more_than_9";
			_shape_of_inner_sphere="shape_coord_more_than_6";
		}

		delete[] pt;
	}
}

vector<int> metal_binding_site::count_planar_atoms(){
	double det;
	int size=_inner_coord_sphere.size();
	vector<int> _planar_atom_nos;
	vector<Point3D> all_points;

	for (int x=0; x<size; x++){
		atom atmp=_inner_coord_sphere.at(x).get_atom_direct_contact();
		Point3D pt=Point3D(atmp.X(), atmp.Y(), atmp.Z());
		all_points.push_back(pt);
	}

	vector<int> range;
	for (int x=0; x<size; x++){
		range.push_back(x);
	}

	vector<vector<int> > all_comb;

	store_combination(range, size, 4, all_comb);

	vector<vector<int> > probable_planar;
	vector<double> planar_det;

	cout<<"\n\nCompute Determinant: ";
	for (auto comb: all_comb){
		vector<Point3D> sub_points;

		for (auto c: comb){
			sub_points.push_back(all_points.at(c));
		}

		bool is_coplanar=check_for_coplanarity(sub_points, planar_det);

		if (is_coplanar==true){
			cout<<"\nPlanar Atoms:";
			for (auto c: comb){
				cout<<c<<", ";
			}

			probable_planar.push_back(comb);
		}
	}

	if (probable_planar.size()!=0){
		for (int x=0; x<planar_det.size()-1; x++){
			for (int y=(x+1); y<planar_det.size(); y++){
				if (planar_det.at(x)>planar_det.at(y)){
					double temp=planar_det.at(x);
					planar_det.at(x)=planar_det.at(y);
					planar_det.at(y)=temp;

					vector<int> swap1=probable_planar.at(x);
					probable_planar.at(x).clear();
					vector<int> swap2=probable_planar.at(y);
					probable_planar.at(y).clear();
					for (auto pp: swap2){
						probable_planar.at(x).push_back(pp);
					}

					for (auto pp: swap1){
						probable_planar.at(y).push_back(pp);
					}
				}
			}
		}

		vector<int> planar_index(probable_planar.at(0).begin(), probable_planar.at(0).end());
		vector<int> non_planar_index;
		for (int x=0; x<size; x++){
			int y;
			for (y=0; y<planar_index.size(); y++){
				if (x==planar_index.at(y)){
					break;
				}
			}

			if (y==planar_index.size()){
				non_planar_index.push_back(x);
			}
		}

		cout<<"\nPre planar: ";
		for (auto pi: planar_index){
			cout<<pi<<", ";
		}

		cout<<"\nPre non planar: ";
		for (auto pi: non_planar_index){
			cout<<pi<<", ";
		}

		for (int x=0; x<non_planar_index.size(); x++){
			vector<Point3D> sub_points;

			for (auto pi: planar_index){
				sub_points.push_back(all_points.at(pi));
			}

			sub_points.push_back(all_points.at(non_planar_index.at(x)));
			bool is_coplanar1=check_for_coplanarity(sub_points);

			if (is_coplanar1==true){
				planar_index.push_back(non_planar_index.at(x));
			}
		}

		for (auto pi: planar_index){
			_planar_atom_nos.push_back(pi);
		}

		cout<<"\nPost planar: ";
		for (auto pi: planar_index){
			cout<<pi<<", ";
		}
	}

	if (_planar_atom_nos.size()!=0){
		for (auto pl: _planar_atom_nos){
			_planar_atoms.push_back(_inner_coord_sphere.at(pl).get_atom_direct_contact().get_serial_no());
		}
	}

	//cout<<"\n\nPlanar_Atom_nos: "<<_planar_atom_nos.size();
	return _planar_atom_nos;
}

void metal_binding_site::store_combination(vector<int> range, int size, int x, vector<vector<int> > & all_comb){
	int temp[x];
	all_combinations(range, temp, 0, size-1, 0, x, all_comb);
}

void metal_binding_site::all_combinations(vector<int> range, int temp[], int start, int end, int index, int x, vector<vector<int> > & all_comb){
	if (index==x){
		vector<int> st;
		for (int p=0; p<x; p++){
			st.push_back(temp[p]);
		}

		all_comb.push_back(st);
		return;
	}

	for (int i=start; i<=end && end-i+1>=x-index; i++){
		temp[index]=range.at(i);
		all_combinations(range, temp, i+1, end, index+1, x, all_comb);
	}
}

bool metal_binding_site::check_for_coplanarity(vector<Point3D> points, vector<double> & planar_det){
	double det=compute_determinant(points.at(0), points.at(1), points.at(2), points.at(3));

	cout<<"\nDeterminant Initial Comb: "<<det;
	if (abs(det)<det_threshold){
		planar_det.push_back(abs(det));
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_for_coplanarity(vector<Point3D> points){
	double det=compute_determinant(points.at(0), points.at(1), points.at(2), points.at(3));
	cout<<"\nDeterminant Further: "<<det;
	if (points.size()==4){
		if (abs(det)<det_threshold){
			return true;
		}else{
			return false;
		}
	}else{
		double det1=compute_determinant(points.at(0), points.at(1), points.at(2), points.at(points.size()-1));
		double det2=compute_determinant(points.at(0), points.at(3), points.at(2), points.at(points.size()-1));
		//cout<<"\n\nDeterminant Further Det1: "<<det1<<", Det2: "<<det2;

		if ((abs(det1)<det_threshold) && (abs(det2)<det_threshold)){
			return true;
		}else{
			return false;
		}
	}
}

double metal_binding_site::compute_determinant(Point3D p1, Point3D p2, Point3D p3){
	Matrix<float, 3, 3> det_matrix;
	double det;
	det_matrix<<p1.X(), p2.X(), p3.X(),
				p1.Y(), p2.Y(), p3.Y(),
				p1.Z(), p2.Z(), p3.Z();

	det=det_matrix.determinant();

	return det;
}

double metal_binding_site::compute_determinant(Point3D p1, Point3D p2, Point3D p3, Point3D p4){
	Matrix<float, 4, 4> det_matrix;
	double det;
	det_matrix<<p1.X(), p1.Y(), p1.Z(), 1,
				p2.X(), p2.Y(), p2.Z(), 1,
				p3.X(), p3.Y(), p3.Z(), 1,
				p4.X(), p4.Y(), p4.Z(), 1;

	det=det_matrix.determinant();

	return det;
}

void metal_binding_site::find_non_planar_atom_indices(vector<int> _planar_atom_nos, int* pt){
	int count=0, y;
	int no_of_non_planar_atoms=_inner_coord_no-_planar_atom_nos.size();
	//cout<<"\n\nIn Find_Non_planar, No. of Non_Planar Atoms: "<<no_of_non_planar_atoms<<",  Coordination: "<<_inner_coord_no;

	for (int x=0; x<_inner_coord_no && count<no_of_non_planar_atoms; x++){
		for (y=0; y<_planar_atom_nos.size(); y++){
			if(x!=_planar_atom_nos.at(y)){
				continue;
			}else{
				break;
			}
		}

		if (y==_planar_atom_nos.size()){
			pt[count]=x;
			count++;
			_non_planar_atoms.push_back(_inner_coord_sphere.at(x).get_atom_direct_contact().get_serial_no());
		}
	}
}

bool metal_binding_site::is_same_side(int* pt, vector<int> _planar_atom_nos){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	int no_of_non_planar_atoms=_inner_coord_no-_planar_atom_nos.size();
	vector<Point3D> planar_points;
	vector<Point3D> non_planar_points;
	vector<Point3D> sub_points;

	for (int x=0; x<no_of_non_planar_atoms; x++){
		atom atmnp=_inner_coord_sphere.at(pt[x]).get_atom_direct_contact();
		non_planar_points.push_back(Point3D(atmnp.X(), atmnp.Y(), atmnp.Z()));
	}

	for (auto pli: _planar_atom_nos){
		atom atmp=_inner_coord_sphere.at(pli).get_atom_direct_contact();
		planar_points.push_back(Point3D(atmp.X(), atmp.Y(), atmp.Z()));
	}

	if (non_planar_points.size()==2){
		bool is_same=check_side(planar_points, non_planar_points);
		if (is_same==true){
			return true;
		}else{
			return false;
		}
	}else{
		bool is_coplanar=false;
		if (planar_points.size()==3){
			is_coplanar=true;
		}else if (planar_points.size()>3){
			sub_points.push_back(planar_points.at(0));
			sub_points.push_back(planar_points.at(1));
			sub_points.push_back(planar_points.at(2));
			sub_points.push_back(planar_points.at(3));
			sub_points.push_back(pm);
			for (int x=3; x<planar_points.size(); x++){
				if (x!=3){
					sub_points.at(3)=planar_points.at(x);
				}

				is_coplanar=check_for_coplanarity(sub_points);

				if (is_coplanar==false){
					break;
				}
			}
		}

		if (is_coplanar==true){
			bool is_same=check_side(planar_points, non_planar_points);

			if (is_same==true){
				return true;
			}else{
				return false;
			}
		}else{
			bool intersected=false;
			int count=0;
			for (int x=0; x<planar_points.size(); x++){
				for (int y=0; y<planar_points.size(); y++){
					if (x<y){
						for (int z=0; z<no_of_non_planar_atoms; z++){
							for (int t=0; t<planar_points.size(); t++){
								double angle_planar_non_planar=Point3D::angle_deg_triangle_law(&non_planar_points.at(z), &pm, &planar_points.at(t));

								if (angle_planar_non_planar<180.0){
									intersected=check_intersection_with_plane(planar_points.at(x), planar_points.at(y), pm, non_planar_points.at(z), planar_points.at(t));
									if (intersected==true){
										count++;
										break;
									}
								}else if (angle_planar_non_planar==180.0){
									count++;
									break;
								}
							}
						}
					}

					if (count==no_of_non_planar_atoms){
						break;
					}
				}

				if (count==no_of_non_planar_atoms){
					break;
				}
			}

			if (count==no_of_non_planar_atoms){
				return true;
			}else{
				return false;
			}
		}
	}
}

bool metal_binding_site::check_side(vector<Point3D> planar_points, vector<Point3D> non_planar_points){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	bool intersected=false;
	for (int x=0; x<planar_points.size(); x++){
		for (int y=0; y<planar_points.size(); y++){
			if (x<y){
				for (int z=0; z<non_planar_points.size()-1; z++){
					intersected=check_intersection_with_plane(planar_points.at(x), planar_points.at(y), pm, non_planar_points.at(z), non_planar_points.at(z+1));
					if (intersected==true){
						break;
					}
				}
			}

			if (intersected==true){
				break;
			}
		}

		if (intersected==true){
			break;
		}
	}

	if (intersected==false){
		intersected=check_for_hemispheric(planar_points, non_planar_points);
		if (intersected==true){
			return false;
		}else{
			return true;
		}
	}else{
		return false;
	}
}

bool metal_binding_site::is_square(Point3D p1, Point3D p2, Point3D p3, Point3D p4){
	double dist1=p1.dist(&p2);
	double dist2=p2.dist(&p3);
	double dist3=p3.dist(&p4);
	double dist4=p4.dist(&p1);

	if ((dist1==dist2) && (dist1==dist3) && (dist1==dist4)){
		double angle1=Point3D::angle_deg_triangle_law(&p1, &p2, &p3);
		double angle2=Point3D::angle_deg_triangle_law(&p2, &p3, &p4);
		double angle3=Point3D::angle_deg_triangle_law(&p3, &p4, &p1);
		double angle4=Point3D::angle_deg_triangle_law(&p4, &p1, &p2);

		if (((angle1==90.0)) && (angle1==angle2) && ((angle1==angle3)) && (angle1==angle4)){
			return true;
		}else{
			return false;
		}
	}else{
		return false;
	}
}

bool metal_binding_site::is_cube(int* pt, vector<int> _planar_atom_nos){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	bool square_planar=is_square(p1, p2, p3, p4);

	if (square_planar==false){
		return false;
	}

	vector<Point3D> non_planar_points;
	int no_of_non_planar=_inner_coord_no-_planar_atom_nos.size();
	for (int x=0; x<no_of_non_planar; x++){
		atom atmnp=_inner_coord_sphere.at(pt[x]).get_atom_direct_contact();
		Point3D ptt=Point3D(atmnp.X(), atmnp.Y(), atmnp.Z());
		non_planar_points.push_back(ptt);
	}

	bool is_coplanar1=check_for_coplanarity(non_planar_points);

	if (is_coplanar1==false){
		return false;
	}

	bool square_non_planar=is_square(non_planar_points.at(0), non_planar_points.at(1), non_planar_points.at(2), non_planar_points.at(3));

	if (square_non_planar==false){
		return false;
	}

	double distp1p2=p1.dist(&p2);
	double distnp1np2=non_planar_points.at(0).dist(&non_planar_points.at(1));

	if (distp1p2!=distnp1np2){
		return false;
	}

	double distp1np1=p1.dist(&non_planar_points.at(0));
	double distp1np2=p1.dist(&non_planar_points.at(1));
	double distp1np3=p1.dist(&non_planar_points.at(2));
	double distp1np4=p1.dist(&non_planar_points.at(3));

	if ((distp1np1!=distp1p2) && (distp1np2!=distp1p2) && (distp1np3!=distp1p2) && (distp1np4!=distp1p2)){
		return false;
	}

	double angle1=Point3D::angle_deg_triangle_law(&p1, &non_planar_points.at(0), &non_planar_points.at(1));
	double angle2=Point3D::angle_deg_triangle_law(&p1, &non_planar_points.at(1), &non_planar_points.at(2));
	double angle3=Point3D::angle_deg_triangle_law(&p1, &non_planar_points.at(2), &non_planar_points.at(3));
	double angle4=Point3D::angle_deg_triangle_law(&p1, &non_planar_points.at(3), &non_planar_points.at(0));

	if (angle1==90.0){
		double angle1a=Point3D::angle_deg_triangle_law(&p2, &non_planar_points.at(1), &non_planar_points.at(2));
		if (angle1a==90.0){
			double angle1a1=Point3D::angle_deg_triangle_law(&p3, &non_planar_points.at(2), &non_planar_points.at(3));
			if (angle1a1==90.0){
				double angle1a1b=Point3D::angle_deg_triangle_law(&p4, &non_planar_points.at(3), &non_planar_points.at(0));
				if (angle1a1b==90.0){
					return true;
				}else{
					return false;
				}
			}else{
				return false;
			}
		}else{
			return false;
		}
	}else if (angle2==90.0){
		double angle2a=Point3D::angle_deg_triangle_law(&p2, &non_planar_points.at(2), &non_planar_points.at(3));
		if (angle2a==90.0){
			double angle2a2=Point3D::angle_deg_triangle_law(&p3, &non_planar_points.at(3), &non_planar_points.at(0));
			if (angle2a2==90.0){
				double angle2a2b=Point3D::angle_deg_triangle_law(&p4, &non_planar_points.at(0), &non_planar_points.at(1));
				if (angle2a2b==90.0){
					return true;
				}else{
					return false;
				}
			}else{
				return false;
			}
		}else{
			return false;
		}
	}else if (angle3==90.0){
		double angle3a=Point3D::angle_deg_triangle_law(&p2, &non_planar_points.at(3), &non_planar_points.at(0));
		if (angle3a==90.0){
			double angle3a3=Point3D::angle_deg_triangle_law(&p3, &non_planar_points.at(0), &non_planar_points.at(1));
			if (angle3a3==90.0){
				double angle3a3b=Point3D::angle_deg_triangle_law(&p4, &non_planar_points.at(1), &non_planar_points.at(2));
				if (angle3a3b==90.0){
					return true;
				}else{
					return false;
				}
			}else{
				return false;
			}
		}else{
			return false;
		}
	}else if (angle4==90.0){
		double angle4a=Point3D::angle_deg_triangle_law(&p2, &non_planar_points.at(0), &non_planar_points.at(1));
		if (angle4a==90.0){
			double angle4a4=Point3D::angle_deg_triangle_law(&p3, &non_planar_points.at(1), &non_planar_points.at(2));
			if (angle4a4==90.0){
				double angle4a4b=Point3D::angle_deg_triangle_law(&p4, &non_planar_points.at(2), &non_planar_points.at(3));
				if (angle4a4b==90.0){
					return true;
				}else{
					return false;
				}
			}else{
				return false;
			}
		}else{
			return false;
		}
	}else{
		return false;
	}
}

bool metal_binding_site::check_for_trigonal_bipyramidal(int* pt, vector<int> _planar_atom_nos){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());

	atom atmnp1=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
	atom atmnp2=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();
	Point3D np1=Point3D(atmnp1.X(), atmnp1.Y(), atmnp1.Z());
	Point3D np2=Point3D(atmnp2.X(), atmnp2.Y(), atmnp2.Z());

	bool is_intersected=check_intersection_with_plane(p1, p2, p3, np1, np2);

	if (is_intersected==true){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_for_hemispheric(vector<Point3D> planar_points, vector<Point3D> non_planar_points){
	int no_of_non_planar=_inner_coord_no-planar_points.size();
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	bool intersected=false;
	for (int x=0; x<no_of_non_planar; x++){
		for (int y=0; y<no_of_non_planar; y++){
			if (x<y){
				for (int z=0; z<planar_points.size()-1; z++){
					intersected=check_intersection_with_plane(non_planar_points.at(x), non_planar_points.at(y), pm, planar_points.at(z), planar_points.at(z+1));
					if (intersected==true){
						break;
					}
				}
			}

			if (intersected==true){
				break;
			}
		}

		if (intersected==true){
			break;
		}
	}

	if (intersected==true){
		return false;
	}else{
		return true;
	}
}

bool metal_binding_site::is_square_prism(int* pt, vector<int> _planar_atom_nos){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	bool square_planar=is_square(p1, p2, p3, p4);

	if (square_planar==false){
		return false;
	}

	atom atma=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
	atom atmb=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();
	atom atmc=_inner_coord_sphere.at(pt[2]).get_atom_direct_contact();
	atom atmd=_inner_coord_sphere.at(pt[3]).get_atom_direct_contact();
	Point3D np1=Point3D(atma.X(), atma.Y(), atma.Z());
	Point3D np2=Point3D(atmb.X(), atmb.Y(), atmb.Z());
	Point3D np3=Point3D(atmc.X(), atmc.Y(), atmc.Z());
	Point3D np4=Point3D(atmd.X(), atmd.Y(), atmd.Z());

	vector<Point3D> points{np1, np2, np3, np4};
	bool is_coplanar=check_for_coplanarity(points);

	if (is_coplanar==false){
		return false;
	}

	bool square_non_planar=is_square(np1, np2, np3, np4);

	if (square_non_planar==false){
		return false;
	}

	double distp1p2=p1.dist(&p2);
	double distnp1np2=np1.dist(&np2);

	if (distp1p2!=distnp1np2){
		return false;
	}

	double angle1=Point3D::angle_deg_triangle_law(&p1, &np1, &np2);
	double angle2=Point3D::angle_deg_triangle_law(&p1, &np2, &np3);
	double angle3=Point3D::angle_deg_triangle_law(&p1, &np3, &np4);
	double angle4=Point3D::angle_deg_triangle_law(&p1, &np4, &np1);

	if (angle1==90.0){
		bool ok=check_angles_and_distances(p1, np1, p2, np2, p3, np3, p4, np4);
		if (ok==true){
			return true;
		}else{
			return false;
		}
	}else if (angle2==90.0){
		bool ok=check_angles_and_distances(p1, np2, p2, np3, p3, np4, p4, np1);
		if (ok==true){
			return true;
		}else{
			return false;
		}
	}else if (angle3==90.0){
		bool ok=check_angles_and_distances(p1, np3, p2, np4, p3, np1, p4, np2);
		if (ok==true){
			return true;
		}else{
			return false;
		}
	}else if (angle4==90.0){
		bool ok=check_angles_and_distances(p1, np4, p2, np1, p3, np2, p4, np3);
		if (ok==true){
			return true;
		}else{
			return false;
		}
	}else{
		return false;
	}
}

bool metal_binding_site::check_angles_and_distances(Point3D p1, Point3D np1, Point3D p2, Point3D np2, Point3D p3, Point3D np3, Point3D p4, Point3D np4){
	double distp1np1=p1.dist(&np1);
	double distp2np2=p2.dist(&np2);
	double distp3np3=p3.dist(&np3);
	double distp4np4=p4.dist(&np4);

	if ((distp1np1==distp2np2) && (distp1np1==distp3np3) && (distp1np1==distp4np4)){
		double angle1=Point3D::angle_deg_triangle_law(&p2, &np2, &np3);
		if (angle1==90.0){
			double angle1a=Point3D::angle_deg_triangle_law(&p3, &np3, &np4);
			if (angle1a==90.0){
				double angle1ab=Point3D::angle_deg_triangle_law(&p4, &np4, &np1);
				if (angle1ab==90.0){
					return true;
				}else{
					return false;
				}
			}else{
				return false;
			}
		}else{
			return false;
		}
	}else{
		return false;
	}
}

bool metal_binding_site::test_for_perfect_seesaw(){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	bool is_seesaw=false;
	bool check_seesaw=false;
	double angle;
	atom axial1, axial2, equatorial1, equatorial2;
	Point3D axial_pt1, axial_pt2;
	Point3D equatorial_pt1, equatorial_pt2;

	for (int x=0; x<_inner_coord_sphere.size()-1; x++){
		atom atm1=_inner_coord_sphere.at(x).get_atom_direct_contact();
		Point3D p1=Point3D(atm1.X(), atm1.Y(), atm1.Z());

		for (int y=(x+1); y<_inner_coord_sphere.size(); y++){
			atom atm2=_inner_coord_sphere.at(y).get_atom_direct_contact();
			Point3D p2=Point3D(atm2.X(), atm2.Y(), atm2.Z());
			angle=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);

			if (angle==180.0){
				check_seesaw=true;
				axial_pt1=p1;
				axial_pt2=p2;
				axial1=atm1;
				axial2=atm2;
				int index;

				for (int z=0; z<_inner_coord_sphere.size(); z++){
					if (z!=x && z!=y){
						atom atm3=_inner_coord_sphere.at(z).get_atom_direct_contact();
						equatorial_pt1=Point3D(atm3.X(), atm3.Y(), atm3.Z());
						equatorial1=atm3;
						index=z;
						break;
					}
				}

				for (int z=0; z<_inner_coord_sphere.size(); z++){
					if (z!=x && z!=y && z!=index){
						atom atm4=_inner_coord_sphere.at(z).get_atom_direct_contact();
						equatorial_pt2=Point3D(atm4.X(), atm4.Y(), atm4.Z());
						equatorial2=atm4;
						break;
					}
				}
			}

			if (check_seesaw==true){
				break;
			}
		}

		if (check_seesaw==true){
			break;
		}
	}

	if (check_seesaw==true){
		double equatorial_angle=Point3D::angle_deg_triangle_law(&equatorial_pt1, &pm, &equatorial_pt2);
		if (equatorial_angle>90.0 && equatorial_angle<120.0){
			double equatorial_axial_angle1=Point3D::angle_deg_triangle_law(&equatorial_pt1, &pm, &axial_pt1);
			if (equatorial_axial_angle1==90.0){
				double equatorial_axial_angle2=Point3D::angle_deg_triangle_law(&equatorial_pt1, &pm, &axial_pt2);
				double equatorial_axial_angle3=Point3D::angle_deg_triangle_law(&equatorial_pt2, &pm, &axial_pt1);
				double equatorial_axial_angle4=Point3D::angle_deg_triangle_law(&equatorial_pt2, &pm, &axial_pt2);
				if ((equatorial_axial_angle1==equatorial_axial_angle2) && (equatorial_axial_angle1==equatorial_axial_angle3) && (equatorial_axial_angle1==equatorial_axial_angle4)){
					is_seesaw=true;
					_planar_atoms.push_back(equatorial1.get_serial_no());
					_planar_atoms.push_back(equatorial2.get_serial_no());
					_non_planar_atoms.push_back(axial1.get_serial_no());
					_non_planar_atoms.push_back(axial2.get_serial_no());
				}
			}
		}
	}

	return is_seesaw;
}

bool metal_binding_site::test_for_seesaw(){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	bool is_seesaw=false;
	bool check_seesaw=false;
	double axial_angle;
	atom axial1, axial2, equatorial1, equatorial2;
	Point3D axial_pt1, axial_pt2;
	Point3D equatorial_pt1, equatorial_pt2;

	for (int x=0; x<_inner_coord_sphere.size()-1; x++){
		atom atm1=_inner_coord_sphere.at(x).get_atom_direct_contact();
		Point3D p1=Point3D(atm1.X(), atm1.Y(), atm1.Z());

		for (int y=(x+1); y<_inner_coord_sphere.size(); y++){
			atom atm2=_inner_coord_sphere.at(y).get_atom_direct_contact();
			Point3D p2=Point3D(atm2.X(), atm2.Y(), atm2.Z());
			axial_angle=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);

			if (axial_angle>170.0 && axial_angle<190.0){
				check_seesaw=true;
				axial_pt1=p1;
				axial_pt2=p2;
				axial1=atm1;
				axial2=atm2;

				int index;
				for (int z=0; z<_inner_coord_sphere.size(); z++){
					if (z!=x && z!=y){
						atom atm3=_inner_coord_sphere.at(z).get_atom_direct_contact();
						equatorial_pt1=Point3D(atm3.X(), atm3.Y(), atm3.Z());
						equatorial1=atm3;
						index=z;
						break;
					}
				}

				for (int z=0; z<_inner_coord_sphere.size(); z++){
					if (z!=x && z!=y && z!=index){
						atom atm4=_inner_coord_sphere.at(z).get_atom_direct_contact();
						equatorial_pt2=Point3D(atm4.X(), atm4.Y(), atm4.Z());
						equatorial2=atm4;
						break;
					}
				}
			}

			if (check_seesaw==true){
				break;
			}
		}

		if (check_seesaw==true){
			break;
		}
	}

	if (check_seesaw==true){
		double equatorial_angle=Point3D::angle_deg_triangle_law(&equatorial_pt1, &pm, &equatorial_pt2);
		if (equatorial_angle>90.0 && equatorial_angle<120.0){
			double equatorial_axial_angle1=Point3D::angle_deg_triangle_law(&equatorial_pt1, &pm, &axial_pt1);
			if (equatorial_axial_angle1==(axial_angle/2.0)){
				double equatorial_axial_angle2=Point3D::angle_deg_triangle_law(&equatorial_pt1, &pm, &axial_pt2);
				double equatorial_axial_angle3=Point3D::angle_deg_triangle_law(&equatorial_pt2, &pm, &axial_pt1);
				double equatorial_axial_angle4=Point3D::angle_deg_triangle_law(&equatorial_pt2, &pm, &axial_pt2);

				if ((equatorial_axial_angle1==equatorial_axial_angle2) && (equatorial_axial_angle1==equatorial_axial_angle3) && (equatorial_axial_angle1==equatorial_axial_angle4)){
					is_seesaw=true;
					_planar_atoms.push_back(equatorial1.get_serial_no());
					_planar_atoms.push_back(equatorial2.get_serial_no());
					_non_planar_atoms.push_back(axial1.get_serial_no());
					_non_planar_atoms.push_back(axial2.get_serial_no());
				}
			}
		}
	}

	return is_seesaw;
}

bool metal_binding_site::check_for_octahedral(int* pt, vector<int> _planar_atom_nos){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	vector<Point3D> points1{_inner_coord_sphere.at(pt[0]).get_atom_direct_contact(), _inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact(), _inner_coord_sphere.at(pt[1]).get_atom_direct_contact(), _inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact()};
	bool is_coplanar1=check_for_coplanarity(points1);
	points1.push_back(pm);
	bool is_coplanar2=check_for_coplanarity(points1);

	vector<Point3D> points2{_inner_coord_sphere.at(pt[0]).get_atom_direct_contact(), _inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact(), _inner_coord_sphere.at(pt[1]).get_atom_direct_contact(), _inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact()};
	bool is_coplanar3=check_for_coplanarity(points2);
	points2.push_back(pm);
	bool is_coplanar4=check_for_coplanarity(points2);

	if ((is_coplanar1==true) && (is_coplanar2==true) && (is_coplanar3==true) && (is_coplanar4==true)){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_for_capped_octahedral(int npl, int* pt, vector<int> _planar_atom_nos){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	atom atmnp=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
	atom atmnq=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();
	Point3D np1=Point3D(atmnp.X(), atmnp.Y(), atmnp.Z());
	Point3D np2=Point3D(atmnq.X(), atmnq.Y(), atmnq.Z());
	atom atmnpl=_inner_coord_sphere.at(npl).get_atom_direct_contact();

	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	Point3D np=Point3D(atmnpl.X(), atmnpl.Y(), atmnpl.Z());
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());

	bool intersected1=check_intersection_with_plane(p1, p2, np1, pm, np);
	bool intersected2=check_intersection_with_plane(p1, p2, np2, pm, np);
	bool intersected3=check_intersection_with_plane(p2, p3, np1, pm, np);
	bool intersected4=check_intersection_with_plane(p2, p3, np2, pm, np);
	bool intersected5=check_intersection_with_plane(p3, p4, np1, pm, np);
	bool intersected6=check_intersection_with_plane(p3, p4, np2, pm, np);
	bool intersected7=check_intersection_with_plane(p4, p1, np1, pm, np);
	bool intersected8=check_intersection_with_plane(p4, p1, np2, pm, np);

	if ((intersected1==true) || (intersected2==true) || (intersected3==true) || (intersected4==true) || (intersected5==true) || (intersected6==true) || (intersected7==true) || (intersected8==true)){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_for_trigonal_prismatic(int* pt, vector<int> _planar_atom_nos){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());

	bool planar_side1=check_side(p1, p2, p3, p4);

	if (planar_side1==true){
		atom atmnp=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
		atom atmnq=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();
		Point3D np1=Point3D(atmnp.X(), atmnp.Y(), atmnp.Z());
		Point3D np2=Point3D(atmnq.X(), atmnq.Y(), atmnq.Z());
		vector<Point3D> points1{p1, p2, np1, np2};
		bool is_coplanar1=check_for_coplanarity(points1);
		vector<Point3D> points2{p3, p4, np1, np2};
		bool is_coplanar2=check_for_coplanarity(points2);

		if ((is_coplanar1==true) && (is_coplanar2==true)){
			bool planar_side2=check_side(p1, p2, np1, np2);
			bool planar_side3=check_side(p3, p4, np1, np2);
				if ((planar_side2==true) && (planar_side3==true)){
					return true;
				}else{
					return false;
				}
		}else{
			vector<Point3D> points3{p1, p4, np1, np2};
			bool is_coplanar3=check_for_coplanarity(points3);
			vector<Point3D> points4{p2, p3, np1, np2};
			bool is_coplanar4=check_for_coplanarity(points4);

			if ((is_coplanar3==true) && (is_coplanar4==true)){
				bool planar_side2=check_side(p1, p4, np1, np2);
				bool planar_side3=check_side(p2, p3, np1, np2);
				if ((planar_side2==true) && (planar_side3==true)){
					return true;
				}else{
					return false;
				}
			}else{
				return false;
			}
		}
	}else{
		return false;
	}
}

bool metal_binding_site::check_side(Point3D p1, Point3D p2, Point3D p3, Point3D p4){
	bool parallel1=is_parallel(p1, p2, p3, p4);
	bool parallel2=is_parallel(p1, p4, p2, p3);

	if ((parallel1==true) && (parallel2==true)){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::is_parallel(Point3D p1, Point3D p2, Point3D np1, Point3D np2){
	double r_planar_line=sqrt(((p2.X()-p1.X())*(p2.X()-p1.X()))+((p2.Y()-p1.Y())*(p2.Y()-p1.Y()))+((p2.Z()-p1.Z())*(p2.Z()-p1.Z())));
	double dcl_pl=((p2.X()-p1.X())/r_planar_line);
	double dcm_pl=((p2.Y()-p1.Y())/r_planar_line);
	double dcn_pl=((p2.Z()-p1.Z())/r_planar_line);

	double r_non_planar_line=sqrt(((np2.X()-np1.X())*(np2.X()-np1.X()))+((np2.Y()-np1.Y())*(np2.Y()-np1.Y()))+((np2.Z()-np1.Z())*(np2.Z()-np1.Z())));
	double dcl_npl=((np2.X()-np1.X())/r_non_planar_line);
	double dcm_npl=((np2.Y()-np1.Y())/r_non_planar_line);
	double dcn_npl=((np2.Z()-np1.Z())/r_non_planar_line);

	if (dcl_pl==dcl_npl && dcm_pl==dcm_npl && dcn_pl==dcn_npl){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_for_capped_trigonal_prismatic(int npl, vector<int> _planar_atom_nos){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	atom atmnp=_inner_coord_sphere.at(npl).get_atom_direct_contact();

	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	Point3D np=Point3D(atmnp.X(), atmnp.Y(), atmnp.Z());
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());

	bool intersected=check_intersection_with_plane(p1, p2, p3, pm, np);
	if (intersected){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_position_for_tetragonal_face_bicapped(vector<int> non_pl, vector<int> _planar_atom_nos, int* pt){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	atom atmn1=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
	atom atmn2=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();

	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	Point3D np1=Point3D(atmn1.X(), atmn1.Y(), atmn1.Z());
	Point3D np2=Point3D(atmn2.X(), atmn2.Y(), atmn2.Z());

	bool ok1=check_cap_sides(p1, p2, np1, p3, p4, np1, non_pl);
	bool ok2=check_cap_sides(p1, p4, np1, p2, p3, np1, non_pl);

	if ((ok1==true) || (ok2==true)){
		return true;
	}else{
		bool ok3=check_cap_sides(p1, p2, np1, p1, p2, p4, non_pl);
		bool ok4=check_cap_sides(p1, p4, np1, p1, p4, p2, non_pl);

		if ((ok3==true) || (ok4==true)){
			return true;
		}else{
			bool ok5=check_cap_sides(p3, p4, np1, p3, p4, p2, non_pl);
			bool ok6=check_cap_sides(p2, p3, np1, p2, p3, p1, non_pl);

			if ((ok5==true) || (ok6==true)){
				return true;
			}else{
				return false;
			}
		}
	}
}

bool metal_binding_site::check_cap_sides(Point3D p1, Point3D p2, Point3D p3, Point3D pp1, Point3D pp2, Point3D pp3, vector<int> non_pl){
	atom atmn1=_inner_coord_sphere.at(non_pl.at(0)).get_atom_direct_contact();
	atom atmn2=_inner_coord_sphere.at(non_pl.at(3)).get_atom_direct_contact();
	Point3D np1=Point3D(atmn1.X(), atmn1.Y(), atmn1.Z());
	Point3D np2=Point3D(atmn2.X(), atmn2.Y(), atmn2.Z());
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());

	bool intersected1=check_intersection_with_plane(p1, p2, p3, pm, np1);
	bool intersected2=check_intersection_with_plane(pp1, pp2, pp3, pm, np2);

	if (intersected1==true && intersected2==true){
		return true;
	}else{
		intersected1=check_intersection_with_plane(p1, p2, p3, pm, np2);
		intersected2=check_intersection_with_plane(p1, p2, p3, pm, np1);

		if (intersected1==true && intersected2==true){
			return true;
		}else{
			return false;
		}
	}

}

bool metal_binding_site::check_intersection_with_plane(Point3D p1, Point3D p2, Point3D p3, Point3D np1, Point3D np2){
	Point3D p12=p1.minus(&p2);
	Point3D p13=p1.minus(&p3);
	Point3D normal_to_plane=p12.cross(&p13);

	Point3D direction_vector_np1_np2=Point3D((np2.X()-np1.X()), (np2.Y()-np1.Y()), (np2.Z()-np1.Z()));
	double dot=direction_vector_np1_np2.dot(&normal_to_plane);

	if (dot!=0){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_position_for_triangular_face_bicapped(vector<int> non_pl, vector<int> _planar_atom_nos, int* pt){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	atom atmn1=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
	atom atmn2=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();

	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	Point3D np1=Point3D(atmn1.X(), atmn1.Y(), atmn1.Z());
	Point3D np2=Point3D(atmn2.X(), atmn2.Y(), atmn2.Z());

	bool ok1=check_cap_sides(p1, p4, np1, p2, p3, np2, non_pl);
	bool ok2=check_cap_sides(p1, p4, np2, p2, p3, np1, non_pl);

	if ((ok1==true) && (ok2==true)){
		return true;
	}else{
		bool ok3=check_cap_sides(p1, p2, np1, p3, p4, np2, non_pl);
		bool ok4=check_cap_sides(p1, p2, np2, p3, p4, np1, non_pl);

		if ((ok3==true) && (ok4==true)){
			return true;
		}else{
			return false;
		}
	}
}

bool metal_binding_site::check_position_for_tricapped(vector<int> non_pl, vector<int> _planar_atom_nos, int* pt){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	atom atmn1=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
	atom atmn2=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();

	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	Point3D np1=Point3D(atmn1.X(), atmn1.Y(), atmn1.Z());
	Point3D np2=Point3D(atmn2.X(), atmn2.Y(), atmn2.Z());
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());

	bool intersected1=check_cap_sides(p1, p2, np1, p3, p4, np1, non_pl);
	bool intersected2=check_cap_sides(p1, p4, np1, p2, p3, np1, non_pl);

	atom atmnp=_inner_coord_sphere.at(non_pl.at(2)).get_atom_direct_contact();
	Point3D np=Point3D(atmnp.X(), atmnp.Y(), atmnp.Z());
	bool intersected3=check_intersection_with_plane(p1, p2, p3, pm, np);

	if (((intersected1==true) || (intersected2==true)) && (intersected3==true)){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_position_for_square_antiprism(int* pt, vector<int> _planar_atom_nos){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	bool square_planar=is_square(p1, p2, p3, p4);

	if (square_planar==false){
		return false;
	}

	atom atma=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
	atom atmb=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();
	atom atmc=_inner_coord_sphere.at(pt[2]).get_atom_direct_contact();
	atom atmd=_inner_coord_sphere.at(pt[3]).get_atom_direct_contact();
	Point3D np1=Point3D(atma.X(), atma.Y(), atma.Z());
	Point3D np2=Point3D(atmb.X(), atmb.Y(), atmb.Z());
	Point3D np3=Point3D(atmc.X(), atmc.Y(), atmc.Z());
	Point3D np4=Point3D(atmd.X(), atmd.Y(), atmd.Z());

	vector<Point3D> points{np1, np2, np3, np4};
	bool is_coplanar=check_for_coplanarity(points);

	if (is_coplanar==false){
		return false;
	}

	bool square_non_planar=is_square(np1, np2, np3, np4);

	if (square_non_planar==false){
		return false;
	}

	double distp1p2=p1.dist(&p2);
	double distnp1np2=np1.dist(&np2);

	if (distp1p2!=distnp1np2){
		return false;
	}

	bool ok=check_equilateral_faces(p1, p2, p3, p4, np1, np2, np3, np4);
	if (ok==true){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_equilateral_faces(Point3D p1, Point3D p2, Point3D p3, Point3D p4, Point3D np1, Point3D np2, Point3D np3, Point3D np4){
	double angle1=Point3D::angle_deg_triangle_law(&p1, &np1, &p2);
	double angle2=Point3D::angle_deg_triangle_law(&p1, &np2, &p2);
	double angle3=Point3D::angle_deg_triangle_law(&p1, &np3, &p2);
	double angle4=Point3D::angle_deg_triangle_law(&p1, &np4, &p2);

	if (angle1==60.0){
		bool sq_ant_pris=check_tri_angles(p1, p2, p3, p4, np1, np2, np3, np4);
		if (sq_ant_pris==true){
			return true;
		}else{
			return false;
		}
	}else if (angle2==60.0){
		bool sq_ant_pris=check_tri_angles(p1, p2, p3, p4, np2, np3, np4, np1);
		if (sq_ant_pris==true){
			return true;
		}else{
			return false;
		}
	}else if (angle3==60.0){
		bool sq_ant_pris=check_tri_angles(p1, p2, p3, p4, np3, np4, np1, np2);
		if (sq_ant_pris==true){
			return true;
		}else{
			return false;
		}
	}else if (angle4==60.0){
		bool sq_ant_pris=check_tri_angles(p1, p2, p3, p4, np4, np1, np2, np3);
		if (sq_ant_pris==true){
			return true;
		}else{
			return false;
		}
	}else{
		return false;
	}
}

bool metal_binding_site::check_tri_angles(Point3D p1, Point3D p2, Point3D p3, Point3D p4, Point3D np1, Point3D np2, Point3D np3, Point3D np4){
	double angle1=Point3D::angle_deg_triangle_law(&np1, &p1, &p2);
	double angle2=Point3D::angle_deg_triangle_law(&p2, &np1, &np2);
	double angle3=Point3D::angle_deg_triangle_law(&np1, &p2, &np2);
	double angle4=Point3D::angle_deg_triangle_law(&p2, &np2, &p3);
	double angle5=Point3D::angle_deg_triangle_law(&np2, &p2, &p3);
	double angle6=Point3D::angle_deg_triangle_law(&p3, &np2, &np3);
	double angle7=Point3D::angle_deg_triangle_law(&np2, &p3, &np3);
	double angle8=Point3D::angle_deg_triangle_law(&p3, &np3, &p4);
	double angle9=Point3D::angle_deg_triangle_law(&np3, &p3, &p4);
	double angle10=Point3D::angle_deg_triangle_law(&p4, &np3, &np4);
	double angle11=Point3D::angle_deg_triangle_law(&np3, &p4, &np4);
	double angle12=Point3D::angle_deg_triangle_law(&p4, &np4, &p1);
	double angle13=Point3D::angle_deg_triangle_law(&np4, &p4, &p1);
	double angle14=Point3D::angle_deg_triangle_law(&p1, &np4, &np1);
	double angle15=Point3D::angle_deg_triangle_law(&np4, &p1, &np1);

	if ((angle1==60.0) && (angle1==angle2) && (angle1==angle3) && (angle1==angle4) && (angle1==angle5) && (angle1==angle6) && (angle1==angle7) && (angle1==angle8) && (angle1==angle9) && (angle1==angle10) && (angle1==angle11) && (angle1==angle12) && (angle1==angle13) && (angle1==angle14) && (angle1==angle15)){
		return true;
	}else{
		return false;
	}
}

bool metal_binding_site::check_for_capped_square_anti_prism(int* pt, vector<int> _planar_atom_nos, int no_of_non_planar){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());
	bool square_planar=is_square(p1, p2, p3, p4);

	if (square_planar==false){
		return false;
	}

	atom atma=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
	atom atmb=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();
	atom atmc=_inner_coord_sphere.at(pt[2]).get_atom_direct_contact();
	atom atmd=_inner_coord_sphere.at(pt[3]).get_atom_direct_contact();
	Point3D* np=new Point3D[4];
	np[0]=Point3D(atma.X(), atma.Y(), atma.Z());
	np[1]=Point3D(atmb.X(), atmb.Y(), atmb.Z());
	np[2]=Point3D(atmc.X(), atmc.Y(), atmc.Z());
	np[3]=Point3D(atmd.X(), atmd.Y(), atmd.Z());

	vector<Point3D> points{np[0], np[1], np[2], np[3]};
	bool is_coplanar=check_for_coplanarity(points);

	if (is_coplanar==false){
		delete[] np;
		return false;
	}

	bool square_non_planar=is_square(np[0], np[1], np[2], np[3]);

	if (square_non_planar==false){
		delete[] np;
		return false;
	}

	double distp1p2=p1.dist(&p2);
	double distnp1np2=np[0].dist(&np[1]);

	if (distp1p2!=distnp1np2){
		delete[] np;
		return false;
	}

	bool found=false;
	int z;
	for (int x=0; x<(no_of_non_planar-_planar_atom_nos.size()+1); x++){
		for (int y=(x+1); y<(no_of_non_planar-_planar_atom_nos.size()+2); y++){
			bool ok=check_equilateral_faces(p1, p2, p3, p4, np[x], np[y], np[y+1], np[y+2]);
			if (ok==true){
				for (z=0; z<no_of_non_planar; z++){
					if ((z!=x) && (z!=y) && (z!=y+1) && (z!=y+2)){
						found=true;
						break;
					}
				}
			}

			if (found==true){
				break;
			}
		}

		if (found==true){
			break;
		}
	}

	if (found==true){
		Point3D np12=np[0].minus(&np[1]);
		Point3D np13=np[0].minus(&np[2]);

		Point3D normal_to_plane=np12.cross(&np13);

		atom atmnpl=_inner_coord_sphere.at(z).get_atom_direct_contact();
		Point3D npl=Point3D(atmnpl.X(), atmnpl.Y(), atmnpl.Z());
		Point3D direction_vector_pm_npl=Point3D((npl.X()-pm.X()), (npl.Y()-pm.Y()), (npl.Z()-pm.Z()));
		double dot=normal_to_plane.dot(&direction_vector_pm_npl);

		delete[] np;
		if (dot!=0){
			return true;
		}else{
			return false;
		}
	}else{
		delete[] np;
		return false;
	}
}

bool metal_binding_site::check_for_dodecahedron(int* pt, vector<int> _planar_atom_nos){
	atom atmp=_inner_coord_sphere.at(_planar_atom_nos.at(0)).get_atom_direct_contact();
	atom atmq=_inner_coord_sphere.at(_planar_atom_nos.at(1)).get_atom_direct_contact();
	atom atmr=_inner_coord_sphere.at(_planar_atom_nos.at(2)).get_atom_direct_contact();
	atom atms=_inner_coord_sphere.at(_planar_atom_nos.at(3)).get_atom_direct_contact();
	Point3D p1=Point3D(atmp.X(), atmp.Y(), atmp.Z());
	Point3D p2=Point3D(atmq.X(), atmq.Y(), atmq.Z());
	Point3D p3=Point3D(atmr.X(), atmr.Y(), atmr.Z());
	Point3D p4=Point3D(atms.X(), atms.Y(), atms.Z());

	double angle1=Point3D::angle_deg_triangle_law(&p1, &p2, &p4);
	double angle2=Point3D::angle_deg_triangle_law(&p2, &p1, &p4);
	double angle3=Point3D::angle_deg_triangle_law(&p4, &p3, &p2);
	double angle4=Point3D::angle_deg_triangle_law(&p2, &p4, &p3);

	double angle5=Point3D::angle_deg_triangle_law(&p1, &p4, &p3);
	double angle6=Point3D::angle_deg_triangle_law(&p3, &p1, &p4);
	double angle7=Point3D::angle_deg_triangle_law(&p1, &p2, &p3);
	double angle8=Point3D::angle_deg_triangle_law(&p2, &p1, &p3);

	if (((angle1==60.0) && (angle1==angle2) && (angle1==angle3) && (angle1==angle4)) || ((angle5==60.0) && (angle5==angle6) && (angle5==angle7) && (angle5==angle8))){
		atom atma=_inner_coord_sphere.at(pt[0]).get_atom_direct_contact();
		atom atmb=_inner_coord_sphere.at(pt[1]).get_atom_direct_contact();
		atom atmc=_inner_coord_sphere.at(pt[2]).get_atom_direct_contact();
		atom atmd=_inner_coord_sphere.at(pt[3]).get_atom_direct_contact();
		Point3D np1=Point3D(atma.X(), atma.Y(), atma.Z());
		Point3D np2=Point3D(atmb.X(), atmb.Y(), atmb.Z());
		Point3D np3=Point3D(atmc.X(), atmc.Y(), atmc.Z());
		Point3D np4=Point3D(atmd.X(), atmd.Y(), atmd.Z());

		vector<Point3D> points{np1, np2, np3, np4};
		bool is_coplanar=check_for_coplanarity(points);

		if (is_coplanar==false){
			return false;
		}

		double anglen1=Point3D::angle_deg_triangle_law(&np1, &np2, &np4);
		double anglen2=Point3D::angle_deg_triangle_law(&np2, &np1, &np4);
		double anglen3=Point3D::angle_deg_triangle_law(&np4, &np3, &np2);
		double anglen4=Point3D::angle_deg_triangle_law(&np2, &np4, &np3);

		double anglen5=Point3D::angle_deg_triangle_law(&np1, &np4, &np3);
		double anglen6=Point3D::angle_deg_triangle_law(&np3, &np1, &np4);
		double anglen7=Point3D::angle_deg_triangle_law(&np1, &np2, &np3);
		double anglen8=Point3D::angle_deg_triangle_law(&np2, &np1, &np3);

		if (((anglen1==60.0) && (anglen1==anglen2) && (anglen1==anglen3) && (anglen1==anglen4)) || ((anglen5==60.0) && (anglen5==anglen6) && (anglen5==anglen7) && (anglen5==anglen8))){
			bool ok=check_equilateral_faces(p1, p2, p3, p4, np1, np2, np3, np4);
			if (ok==true){
				return true;
			}else{
				return false;
			}
		}else{
			return false;
		}
	}else{
		return false;
	}
}

void metal_binding_site::compute_mean_angle_planar_atoms_with_metal(vector<int> planar_atom_nos){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	double sum_angle=0.0;
	for (int x=0; x<planar_atom_nos.size()-1; x++){
		atom atm1=_inner_coord_sphere.at(planar_atom_nos.at(x)).get_atom_direct_contact();
		Point3D p1=Point3D(atm1.X(), atm1.Y(), atm1.Z());

		atom atm2=_inner_coord_sphere.at(planar_atom_nos.at(x+1)).get_atom_direct_contact();
		Point3D p2=Point3D(atm2.X(), atm2.Y(), atm2.Z());
		double angle=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);
		sum_angle+=angle;
	}

	atom atm1=_inner_coord_sphere.at(planar_atom_nos.size()-1).get_atom_direct_contact();
	Point3D p1=Point3D(atm1.X(), atm1.Y(), atm1.Z());
	atom atm2=_inner_coord_sphere.at(planar_atom_nos.at(0)).get_atom_direct_contact();
	Point3D p2=Point3D(atm2.X(), atm2.Y(), atm2.Z());
	double angle=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);
	sum_angle+=angle;

	_mean_planar_base_atoms_angles_with_metal=sum_angle/(double) planar_atom_nos.size();
}

void metal_binding_site::compute_mean_angle_non_planar_atoms_with_metal(int* pt, int no_of_non_planar_atoms){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	double sum_angle=0.0;

	for (int x=0; x<no_of_non_planar_atoms-1; x++){
		atom atm1=_inner_coord_sphere.at(pt[x]).get_atom_direct_contact();
		Point3D p1=Point3D(atm1.X(), atm1.Y(), atm1.Z());

		for (int y=(x+1); y<no_of_non_planar_atoms; y++){
			atom atm2=_inner_coord_sphere.at(pt[y]).get_atom_direct_contact();
			Point3D p2=Point3D(atm2.X(), atm2.Y(), atm2.Z());
			double angle=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);
			sum_angle+=angle;
		}
	}

	_mean_non_planar_atoms_angles_with_metal=sum_angle/(double) no_of_non_planar_atoms;
}

void metal_binding_site::compute_mean_angle_planar_and_non_planar_atoms_with_metal(int* pt, int no_of_non_planar_atoms, vector<int> _planar_atom_nos){
	Point3D pm=Point3D(_metal_atom.X(), _metal_atom.Y(), _metal_atom.Z());
	double sum_angle=0.0;
	int total_no_of_angles=no_of_non_planar_atoms*_planar_atom_nos.size();

	for (int x=0; x<_planar_atom_nos.size(); x++){
		atom atm1=_inner_coord_sphere.at(_planar_atom_nos.at(x)).get_atom_direct_contact();
		Point3D p1=Point3D(atm1.X(), atm1.Y(), atm1.Z());

		for (int y=0; y<no_of_non_planar_atoms; y++){
			atom atm2=_inner_coord_sphere.at(pt[y]).get_atom_direct_contact();
			Point3D p2=Point3D(atm2.X(), atm2.Y(), atm2.Z());
			double angle=Point3D::angle_deg_triangle_law(&p1, &pm, &p2);
			//cout<<"\nAngle: "<<angle;
			sum_angle+=angle;
		}
	}

	_mean_planar_and_non_planar_atoms_angles_with_metal=sum_angle/(double) total_no_of_angles;
}

void metal_binding_site::compute_atom_metal_atom_angles(){
	double sum=0.0;
	int no_of_sides=0;

	if (_inner_coord_sphere.size()>1){
		for (int i=0; i<_inner_coord_sphere.size()-1; i++){
			atom at1=_inner_coord_sphere.at(i).get_atom_direct_contact();
			for (int j=(i+1); j<_inner_coord_sphere.size(); j++){
				atom at2=_inner_coord_sphere.at(j).get_atom_direct_contact();
				double angle=Point3D::angle_deg_triangle_law(&at1, &_metal_atom, &at2);
				no_of_sides++;
				sum+=angle;
			}
		}

		_mean_angle_atom_metal_atom=(sum/no_of_sides);
	}
}

void metal_binding_site::compute_atom_atom_atom_angles(){
	double sum=0.0;
	int no_of_sides=0;

	if (_inner_coord_sphere.size()>1){
		for (int i=0; i<_inner_coord_sphere.size()-1; i++){
			atom at1=_inner_coord_sphere.at(i).get_atom_direct_contact();
			for (int j=(i+1); j<_inner_coord_sphere.size(); j++){
				atom at2=_inner_coord_sphere.at(j).get_atom_direct_contact();
				double angle=Point3D::angle_deg_triangle_law(&at1, &_metal_atom, &at2);
				no_of_sides++;
				sum+=angle;
			}
		}

		_mean_angle_atom_metal_atom=(sum/no_of_sides);
	}
}

void metal_binding_site::compute_inter_atomic_distances(){
	double sum=0.0;
	int no_of_sides=0;
	if (_inner_coord_sphere.size()>1){
		for (int i=0; i<_inner_coord_sphere.size()-1; i++){
			atom at1=_inner_coord_sphere.at(i).get_atom_direct_contact();
			for (int j=(i+1); j<_inner_coord_sphere.size(); j++){
				atom at2=_inner_coord_sphere.at(j).get_atom_direct_contact();
				double dist=at1.dist(&at2);
				no_of_sides++;
				sum+=dist;
			}
		}

		_mean_inter_atomic_distance=(sum/no_of_sides);
	}
}

void metal_binding_site::compute_atom_to_metal_distances(){
	double sum_dist=0.0;
	for (auto datm: _inner_coord_sphere){
		sum_dist+=datm.get_distance_direct_contact();
	}

	_mean_distance_direct_ligands_from_metal=(sum_dist/(double ) _inner_coord_no);
}
