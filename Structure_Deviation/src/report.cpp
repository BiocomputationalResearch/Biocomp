#include "biomolecule.h"

void protein::generate_general_report(string accn){
	FILE *fp;
	string report_file_name=accn+"_general_report.mbsi";
	fp=fopen(report_file_name.c_str(), "w");

	if (fp!=NULL){
		fprintf(fp, "HEADER     PROGRAM NAME: StructureDeviation  VERSION: 1.1");
		fprintf(fp, "\nTITLE      General Report for Statistics about Protein Atoms");

		fprintf(fp, "\nProtein Name: %s", accn.c_str());
		fprintf(fp, "\nREMARK     1");
		fprintf(fp, "\nREMARK     1  Total No. of Metal Ions Present: %lu", _metal_list.size());

		if (_metal_list.size()==0){
			fprintf(fp, "\nREMARK     1  No Metal Ion Present in Protein %s!", accn.c_str());
		}else{
			fprintf(fp, "\nREMARK     1\tName & No. of Each Type of Metal Atom:");
			for (int p=0; p<_unique_metal_list.size(); p++){
				fprintf(fp, "\nREMARK     1\t\tNo. of Each Metal Site Present %s: %d", _unique_metal_list.at(p).c_str(), _count_unique_metal.at(p));
			}

			vector<vector<metal_binding_site> > _all_metal_sites=get_all_metal_sites();
			int remark=2;
			for (auto msites: _all_metal_sites){
				fprintf(fp, "\nREMARK     %d", remark);
				string site_name=msites.at(0).get_metal_atom().get_atom_symbol();
				fprintf(fp, "\nREMARK     %d  No. of %s Sites Present: %lu", remark, site_name.c_str(), msites.size());
				for (int p=0; p<msites.size(); p++){
					fprintf(fp, "\nREMARK     %d", remark);
					fprintf(fp, "\nREMARK     %d %s Binding Site %d", remark, site_name.c_str(), (p+1));

					atom metal_atom=msites.at(p).get_metal_atom();
					fprintf(fp, "\nREMARK     %d\tSr. No.\tChain Name\tResidue Name\tResidue No.\tAtom Name", remark);
					fprintf(fp, "\nREMARK     %d\t%ld\t\t\t%s\t\t\t%s\t\t\t%ld\t\t\t%s", remark, metal_atom.get_serial_no(), metal_atom.get_chain_name().c_str(), metal_atom.get_residue_name().c_str(), metal_atom.get_residue_no(), metal_atom.get_atom_name().c_str());
					fprintf(fp, "\nREMARK     %d\tCoordination Number %s: %d", remark, site_name.c_str(), msites.at(p).get_inner_coord_no());

					if (msites.at(p).get_inner_coord_no()>0){
						fprintf(fp, "\nREMARK     %d\t%s Binding Site of %s", remark, site_name.c_str(), msites.at(p).get_site_type().c_str());
						string shape=msites.at(p).get_shape();
						fprintf(fp, "\nREMARK     %d\tShape of the %s Binding Site: %s", remark, site_name.c_str(), shape.c_str());
						fprintf(fp, "\nREMARK     %d", remark);

						fprintf(fp, "\nREMARK     %d", remark);
						if (msites.at(p).get_water_molecules().size()!=0){
							fprintf(fp, "\nREMARK     %d\tNo. of Water Groups Present in %s Site, Shape %s: %lu", remark, site_name.c_str(), shape.c_str(), msites.at(p).get_water_molecules().size());
						}else{
							fprintf(fp, "\nREMARK     %d\tNo Water Molecule Present in this %s Site, Shape %s", remark, site_name.c_str(), shape.c_str());
						}

						fprintf(fp, "\nREMARK     %d", remark);
						fprintf(fp, "\nREMARK     %d\tInner Coordinating Sphere for %s_%ld:", remark, site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
						fprintf(fp, "\nREMARK     %d\tDirect Ligand Atoms for %s_%ld:", remark, site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
						fprintf(fp, "\nREMARK     %d\tMolecule Type\t\tSr. No.\tChain Name\tResidue Name\t\tResidue Type\tResidue No.\tAtom Name\t\tDistance from %s", remark, site_name.c_str());
						for (int q=0; q<msites.at(p).get_inner_coord_no(); q++){
							direct_ligand dl=msites.at(p).get_inner_coord_sphere().at(q);
							atom atm=dl.get_atom_direct_contact();
							fprintf(fp, "\nREMARK     %d \t%s\t\t%ld\t\t\t%s\t\t\t%s\t\t\t%s\t\t\t%ld\t\t\t%s\t\t%8.3lf", remark, atm.get_molecule_type().c_str(), atm.get_serial_no(), atm.get_chain_name().c_str(), atm.get_residue_name().c_str(), atm.get_residue_type().c_str(), atm.get_residue_no(), atm.get_atom_name().c_str(), dl.get_distance_direct_contact());
						}

						fprintf(fp, "\n\nREMARK     %d\tInter-Atomic Distances in %s_%ld Site:\n", remark, site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
						compute_inter_atomic_distances(fp, remark, site_name, msites.at(p).get_metal_atom().get_serial_no(), shape, msites.at(p).get_inner_coord_sphere());

						fprintf(fp, "\n\nREMARK     %d\tAngles Between Two Atoms Made with Metal in %s_%ld Site:\n", remark, site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
						compute_atom_to_metal_angles(fp, remark, site_name, msites.at(p).get_metal_atom().get_serial_no(), shape, msites.at(p).get_inner_coord_sphere());

						if (msites.at(p).get_inner_coord_no()>3){
							if (msites.at(p).get_planar_atoms().size()>0){
								fprintf(fp, "\n\nREMARK     %d\tEquatorial (Planar) Atom Nos. in %s_%ld Site: ", remark, site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
								vector<long int> pl=msites.at(p).get_planar_atoms();
								for (auto p: pl){
									fprintf(fp, "%ld, ", p);
								}

								if (msites.at(p).get_non_planar_atoms().size()>0){
									fprintf(fp, "\n\nREMARK     %d\tAxial (Non-Planar) Atom Nos. in %s_%ld Site: ", remark, site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
									vector<long int> npl=msites.at(p).get_non_planar_atoms();
									for (auto np: npl){
										fprintf(fp, "%ld, ", np);
									}
								}else{
									fprintf(fp, "\n\nREMARK     %d\tAll Coordinated Atoms are Equatorial (Planar, %s) in %s_%ld Site: ", remark, shape.c_str(), site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
								}
							}else{
								fprintf(fp, "\n\nREMARK     %d\tNo Atoms More than 3 are coplanar (%s) in %s_%ld Site: ", remark, shape.c_str(), site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
							}
						}else{
							fprintf(fp, "\n\nREMARK     %d\tIn %s, %s_%ld Site has less than 4 Atoms", remark, shape.c_str(), site_name.c_str(), msites.at(p).get_metal_atom().get_serial_no());
						}
					}else{
						fprintf(fp, "\nREMARK     %d\tThis %s Site Does Not Coordinated with Any Atom!", remark, site_name.c_str());
					}

					remark++;
				}
			}

			fclose(fp);
		}
	}else{
		cout<<"\nGeneral File Cannot be opened";
	}
}

void protein::compute_inter_atomic_distances(FILE *fp, int remark, string site_name, long int site_no, string shape, vector<direct_ligand> inn_sph){
	double sum=0.0;
	int no_of_sides=0;
	if (inn_sph.size()>1){
		for (int i=0; i<inn_sph.size()-1; i++){
			atom at1=inn_sph.at(i).get_atom_direct_contact();
			for (int j=(i+1); j<inn_sph.size(); j++){
				atom at2=inn_sph.at(j).get_atom_direct_contact();
				double dist=at1.dist(&at2);
				no_of_sides++;
				sum+=dist;
				fprintf(fp, "\nREMARK     %d\tIn %s (Distance of, Shape %s)_%ld, %s_%ld to %s_%ld: %8.3lf", remark, site_name.c_str(), shape.c_str(), site_no, at1.get_atom_name().c_str(), at1.get_serial_no(), at2.get_atom_name().c_str(), at2.get_serial_no(), dist);
			}
		}

		double mean_dist=(sum/no_of_sides);
		fprintf(fp, "\n\nREMARK     %d\tIn %s (Mean Inter-Atomic Distance, Shape %s)_%ld: %8.3lf", remark, site_name.c_str(), shape.c_str(), site_no, mean_dist);
	}else if (inn_sph.size()==1){
		fprintf(fp, "\nREMARK     %d\t(Shape %s), No Inter-Atomic Distance for Single Coordinated Atom, in %s_%ld", remark, shape.c_str(), site_name.c_str(), site_no);
	}
}

void protein::compute_atom_to_metal_angles(FILE *fp, int remark, string site_name, long int site_no, string shape, vector<direct_ligand> inn_sph){
	double sum=0.0;
	int no_of_sides=0;
	atom parent_metal=inn_sph.at(0).get_parent_metal();

	if (inn_sph.size()>1){
		for (int i=0; i<inn_sph.size()-1; i++){
			atom at1=inn_sph.at(i).get_atom_direct_contact();
			for (int j=(i+1); j<inn_sph.size(); j++){
				atom at2=inn_sph.at(j).get_atom_direct_contact();
				double angle=Point3D::angle_deg_triangle_law(&at1, &parent_metal, &at2);
				no_of_sides++;
				sum+=angle;
				fprintf(fp, "\nREMARK     %d\tIn %s (Angle Between, Shape %s)_%ld, %s_%ld and %s_%ld with %s_%ld: %8.3lf", remark, site_name.c_str(), shape.c_str(), site_no, at1.get_atom_name().c_str(), at1.get_serial_no(), at2.get_atom_name().c_str(), at2.get_serial_no(), parent_metal.get_atom_symbol().c_str(), parent_metal.get_serial_no(), angle);
			}
		}

		double mean_angle=(sum/no_of_sides);
		fprintf(fp, "\n\nREMARK     %d\tIn %s (Mean Angle Atom-%s-Atom, Shape %s)_%ld: %8.3lf", remark, site_name.c_str(), parent_metal.get_atom_symbol().c_str(), shape.c_str(), site_no, mean_angle);
	}else if (inn_sph.size()==1){
		fprintf(fp, "\nREMARK     %d\t(Shape %s), No Atom-%s-Atom Angle for Single Coordinated Atom, in %s_%ld", remark, shape.c_str(), parent_metal.get_atom_symbol().c_str(), site_name.c_str(), site_no);
	}
}
