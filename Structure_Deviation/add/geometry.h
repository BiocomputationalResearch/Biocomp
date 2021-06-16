
#ifndef METAL_DETECTOR_GEOMETRY_H
#define METAL_DETECTOR_GEOMETRY_H

#include <iostream>
#include <cstdio>
#include <cmath>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>
#include <assert.h>
#include <array>
#include <vector>

using namespace std;
using  namespace  Eigen;

namespace geom{
#define PI 3.141592654
    class Point3D{
    protected:
        double x;
        double y;
        double z;
    public:
        Point3D(){
            x=0;
            y=0;
            z=0;
        }

        Point3D(double x, double y, double z);
        double dist(Point3D* p);
        double dist_sqr(Point3D* p);
        Point3D minus(Point3D* p);
        double dot(Point3D *p);
        Point3D cross(Point3D *p);
        Point3D shift_to_origin();
        double rotate_about_x_axis();
        double rotate_about_y_axis();
        double rotate_about_z_axis();
        Point3D reflect_about_x_axis();
        Point3D reflect_about_y_axis();
        Point3D reflect_about_z_axis();
        Point3D reflect_about_xy_plane();
        Point3D reflect_about_yz_plane();
        Point3D reflect_about_zx_plane();

        double X();
        double Y();
        double Z();
        void print();
        static double angle_deg_triangle_law(Point3D* a, Point3D* b, Point3D* c);
        static Point3D centroid(Point3D* a, Point3D* b, Point3D* c, Point3D* d);
    };


    class Sphere:public Point3D{
    protected:
        double r;
    public:
        Sphere(){}
        Sphere(double x, double y, double z, double r);
        double radius();
    };

class triangle{
    	int _no_of_vertices;
    	int _no_of_sides;
    	int _no_of_angles_with_circum_centre;
    	int _no_of_vertex_vertex_angles;
    	Point3D _vertices[3];
    	Point3D _circum_centre;
    	double _circum_radius;
    	Point3D _centre_of_mass;
    	double _circum_centre_vertex_distances[3][2];
    	double _vertex_vertex_distances[3][3];
    	double _circum_centre_vertex_angles_tri[3][3];
    	double _vertex_vertex_angles_tri[3][4];
    	double _direction_cosines[3][3];
    	int _vertex_order[3];

    	void calc_projections_on_standard_triangle(triangle &, double *);
    	void calc_projections_on_standard_triangle(triangle &, double[][3], double *);
    	vector<double> compute_mse(vector<double>, vector<double>);

    public:
    	triangle();
    	void free_memory();
    	triangle (Point3D*);
    	//triangle(atom*);

    	void find_vertex_order();
    	void set_vertex_order(int, int, int);
    	double (*measure_direction_cosines_of_sides())[3];

    	void compute_circum_centre_radius();
    	void compute_centre_of_mass();
    	void measure_circum_centre_to_vertex_distances();
    	void measure_vertex_to_vertex_distances();
    	void measure_angles_vertices_with_circum_centre();
    	void measure_vertex_vertex_angles();

    	void deviation_distances(triangle, double *);
    	void deviation_distances(triangle, double[][3], double *);
    	void deviation_angles_with_circum_centre(triangle, double *);
    	void deviation_angles_vertex_vertex(triangle, double *);
    	//void gen_report_metal_binding_sites(triangle ori_stand, triangle ori_comp, triangle ctdrn, FILE* fp);

    	Point3D* get_vertices();
    	int* get_vertex_order();
    	Point3D get_circum_centre();
    	double get_circum_radius();
    	Point3D get_centre_of_mass();
    	double (*get_circum_centre_to_vertex_distances())[2];
    	double (*get_vertex_vertex_distances())[3];
    	double (*get_angles_vertices_with_circum_centre())[3];
    	double (*get_vertex_vertex_angles())[4];
    	double (*get_direction_cosines())[3];
    	//void print(Point3D* vert);
    };

class tetragonalplanar{
 		int _no_of_vertices;
    	int _no_of_sides;
    	int _no_of_vertex_vertex_angles;
    	Point3D _vertices[4];
    	int _no_of_triangles;
    	double _vertex_vertex_distances[4][3];
    	double _vertex_vertex_angles_tri[4][4];
    	int _vertex_order[4];

    	void measure_vertex_vertex_distances();
    	void measure_vertex_vertex_angles();
    	void compute_vertex_order();
    public:
    	tetragonalplanar();
    	tetragonalplanar (Point3D*);
    	//tetragonalplanar(atom*);
    	int* get_vertex_order();

    	//void deviation_distances(triangle, triangle, tetraplMetalicStructure &, double[], double[]);
    	//void deviation_angles_vertex_vertex(triangle, triangle, tetraplMetalicStructure &, double[], double[]);

    	Point3D* get_vertices();
    	double (*get_vertex_vertex_distances())[3];
    	double (*get_vertex_vertex_angles())[4];
    	//void print(Point3D* vert);
    };

class tetrahedron{
	int _no_of_vertices;
	int _no_of_sides;
	int _no_of_faces;
	int _no_of_angles_with_circum_centre;
	int _no_of_vertex_vertex_angles;
	Point3D _vertices[4];
	Point3D _circum_centre;
	double _circum_radius;
	Point3D _centre_of_mass;
	double _circum_centre_vertex_distances[4][2];
	double _vertex_vertex_distances[6][3];
	double _surface_areas[4][4];
	double _max_surface_area;
	double _circum_centre_vertex_angles_tri[6][3];
	double _vertex_vertex_angles_tri[12][4];
	double _direction_cosines[6][3];
	int _vertex_order[4];

public:
	tetrahedron();
	tetrahedron (Point3D []);

	void set_vertex_order(int [], int);
	void set_vertex_order(int, int, int, int);
	void compute_circum_centre_radius();
	void compute_centre_of_mass();
	void measure_circum_centre_to_vertex_distances();
	void measure_vertex_to_vertex_distances();
	void measure_angles_vertices_with_circum_centre();
	void measure_vertex_vertex_angles();
	void measure_surface_areas_and_order_of_vertices();
	void find_vertex_order(int, int, int);
	void calc_projections_on_standard_tetra(tetrahedron &, double []);
	void calc_projections_on_standard_tetra(tetrahedron &, double[][3], double []);
	//vector<double> compute_mse(vector<double>, vector<double>);
	double (*measure_direction_cosines_of_sides())[3];
	void deviation_distances(geom::tetrahedron, double[]);
	void deviation_distances(geom::tetrahedron, double[][3], double[]);
	void deviation_angles_with_circum_centre(geom::tetrahedron, double[]);
	void deviation_angles_vertex_vertex(geom::tetrahedron, double[]);

	Point3D* get_vertices();
	double (*get_surface_areas())[4];
	int* get_vertex_order();
	Point3D get_circum_centre();
	double get_circum_radius();
	Point3D get_centre_of_mass();
	double (*get_circum_centre_to_vertex_distances())[2];
	double (*get_vertex_vertex_distances())[3];
	double (*get_angles_vertices_with_circum_centre())[3];
	double (*get_vertex_vertex_angles())[4];
};

class tetragonal_pyramidal{
    	int _no_of_vertices;
    	int _no_of_sides;
    	int _no_of_planar_atoms;
    	int _no_of_non_planar_atom;
    	int _no_of_vertex_vertex_angles;
    	int _no_of_tetrahedron;
    	Point3D _planar_vertices[4];
    	Point3D _non_planar_vertex;
    	double _vertex_vertex_distances_pl_pl[4][3];
    	double _vertex_vertex_distances_npl_pl[4][3];
    	double _vertex_vertex_angles_tri_pl_pl[4][4];
    	double _vertex_vertex_angles_tri_npl_pl[12][4];
    	int _pl_vertex_order[4];
    	int _npl_vertex_index;

    	//void compute_circum_centre_radius();
    	void measure_vertex_to_vertex_distances();
    	void measure_vertex_vertex_angles();
    	void find_vertex_order();
    	//void calc_projections_on_standard_trigonal_bipyramid(geom::tetrahedron &, double **);
    	//vector<double> compute_mse(vector<double>, vector<double>);

    public:
    	tetragonal_pyramidal();
    	tetragonal_pyramidal (vector<Point3D> plv, Point3D npv);
    	tetragonal_pyramidal(Point3D*, Point3D);

    	Point3D* get_pl_vertices();
    	Point3D get_npl_vertex();
    	int* get_pl_vertex_order();
    	double (*get_vertex_vertex_distances_pl_pl())[3];
    	double (*get_vertex_vertex_distances_npl_pl())[3];
    	double (*get_vertex_vertex_angles_pl_pl())[4];
    	double (*get_vertex_vertex_angles_npl_pl())[4];
    };

class trigonal_bipyramidal{
    	int _no_of_vertices;
    	int _no_of_sides;
    	int _no_of_planar_atoms;
    	int _no_of_non_planar_atoms;
    	int _no_of_vertex_vertex_angles;
    	int _no_of_tetrahedron;
    	Point3D _planar_vertices[3];
    	Point3D _non_planar_vertices[2];
    	double _vertex_vertex_distances_pl_pl[3][3];
    	double _vertex_vertex_distances_npl_pl[6][3];
    	double _vertex_vertex_angles_tri_pl_pl[3][4];
    	double _vertex_vertex_angles_tri_npl_pl[18][4];
    	int _pl_vertex_order[3];
    	int _npl_vertex_order[2];

    	//void compute_circum_centre_radius();
    	void measure_vertex_to_vertex_distances();
    	void measure_vertex_vertex_angles();
    	void find_vertex_order();
    	//void calc_projections_on_standard_trigonal_bipyramid(geom::tetrahedron &, double **);
    	//vector<double> compute_mse(vector<double>, vector<double>);

    public:
    	trigonal_bipyramidal();
    	trigonal_bipyramidal (vector<Point3D> plv, vector<Point3D> npv);
    	trigonal_bipyramidal(Point3D*, Point3D*);

    	Point3D* get_pl_vertices();
    	Point3D* get_npl_vertices();
    	int* get_pl_vertex_order();
    	int* get_npl_vertex_order();
    	double (*get_vertex_vertex_distances_pl_pl())[3];
    	double (*get_vertex_vertex_distances_npl_pl())[3];
    	double (*get_vertex_vertex_angles_pl_pl())[4];
    	double (*get_vertex_vertex_angles_npl_pl())[4];
    };

class octahedral{
    	int _no_of_vertices;
    	int _no_of_sides;
    	int _no_of_planar_atoms;
    	int _no_of_non_planar_atoms;
    	int _no_of_vertex_vertex_angles;
    	int _no_of_tetrahedron;
    	Point3D _planar_vertices[4];
    	Point3D _non_planar_vertices[2];
    	double _vertex_vertex_distances_pl_pl[4][3];
    	double _vertex_vertex_distances_npl_pl[8][3];
    	double _vertex_vertex_angles_tri_pl_pl[4][4];
    	double _vertex_vertex_angles_tri_npl_pl[24][4];
    	int _pl_vertex_order[4];
    	int _npl_vertex_order[2];

    	//void compute_circum_centre_radius();
    	void measure_vertex_to_vertex_distances();
    	void measure_vertex_vertex_angles();
    	void find_vertex_order();
    	//void calc_projections_on_standard_trigonal_bipyramid(geom::tetrahedron &, double **);
    	//vector<double> compute_mse(vector<double>, vector<double>);

    public:
    	octahedral();
    	octahedral (vector<Point3D> plv, vector<Point3D> npv);
    	octahedral(Point3D*, Point3D*);

     	Point3D* get_pl_vertices();
    	Point3D* get_npl_vertices();
    	int* get_pl_vertex_order();
    	int* get_npl_vertex_order();
    	double (*get_vertex_vertex_distances_pl_pl())[3];
    	double (*get_vertex_vertex_distances_npl_pl())[3];
    	double (*get_vertex_vertex_angles_pl_pl())[4];
    	double (*get_vertex_vertex_angles_npl_pl())[4];
    };
}

#endif //METAL_DETECTOR_GEOMETRY_H
