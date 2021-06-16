#include "geometry.h"

geom::Point3D::Point3D(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

double geom::Point3D::dist(geom::Point3D *p) {
    return sqrt ( (x-p->x) * (x-p->x) + (y-p->y) * (y-p->y) + (z-p->z) * (z-p->z) );
}

/*friend ostream& operator<<(ostream &out, const Point3D *p) {
    out<<"("<<p->x<<","<<p->y<<","<<p->z<<")";
    //cout<<"pointer version\n";
    return out;
}*/

double geom::Point3D::dist_sqr(geom::Point3D *p) {
    return (x-p->x) * (x-p->x) + (y-p->y) * (y-p->y) + (z-p->z) * (z-p->z) ;

}

geom::Point3D geom::Point3D::minus(geom::Point3D* p) {
    return geom::Point3D(p->x - x, p->y - y, p->z - z);
}

double geom::Point3D::dot(geom::Point3D *p) {
    return (x * p->x + y * p->y + z * p->z) ;
}

geom::Point3D geom::Point3D::cross(geom::Point3D *p) {
	double x_comp=y*p->z-z*p->y;
	double y_comp=x*p->z-z*p->x;
	double z_comp=x*p->y-y*p->x;

	Point3D vect=Point3D(x_comp, y_comp, z_comp);
    return vect ;
}

geom::Point3D geom::Point3D::shift_to_origin(){
	double tx=-x;
	double ty=-y;
	double tz=-z;

	return Point3D(tx, ty, tz);
}

double geom::Point3D::rotate_about_x_axis(){
	double rotate_amount=atan(-z / y);
	return rotate_amount;
}

double geom::Point3D::rotate_about_y_axis(){
	double rotate_amount=atan(-z / x);
	return rotate_amount;
}

double geom::Point3D::rotate_about_z_axis(){
	double rotate_amount=atan(-y / x);
	return rotate_amount;
}

geom::Point3D geom::Point3D::reflect_about_x_axis(){
	return Point3D(x, -y, -z);
}

geom::Point3D geom::Point3D::reflect_about_y_axis(){
	return Point3D(-x, y, -z);
}

geom::Point3D geom::Point3D::reflect_about_z_axis(){
	return Point3D(-x, -y, z);
}

geom::Point3D geom::Point3D::reflect_about_xy_plane(){
	return Point3D(x, y, -z);
}

geom::Point3D geom::Point3D::reflect_about_yz_plane(){
	return Point3D(-x, y, z);
}

geom::Point3D geom::Point3D::reflect_about_zx_plane(){
	return Point3D(x, -y, z);
}

double geom::Point3D::X() {
    return x;
}

double geom::Point3D::Y() {
    return y;
}

double geom::Point3D::Z() {
    return z;
}

void geom::Point3D::print() {
    cout<<"("<<x<<","<<y<<","<<z<<")\n";

}

double geom::Point3D::angle_deg_triangle_law(geom::Point3D *a, geom::Point3D *b, geom::Point3D *c) {
    double dist_sqr_ab = a->dist_sqr(b);
    double dist_sqr_bc = b->dist_sqr(c);
    double dist_sqr_ca   = c->dist_sqr(a);
    double dist_ab=a->dist(b);
    double dist_bc=b->dist(c);
    double numerator = (dist_sqr_ab+dist_sqr_bc-dist_sqr_ca);
    double denominator = 2 * dist_ab * dist_bc;
    double angle;
    if ((numerator/denominator)>1.0){
    	angle=acos(1);
    }else{
    	angle = acos(numerator / denominator);
    }

    angle*=(180.0*7.0)/22.0;
    return angle;
}

geom::Point3D geom::Point3D::centroid(geom::Point3D *a, geom::Point3D *b, geom::Point3D *c, geom::Point3D *d) {
	double _cent_x;
    double _cent_y;
    double _cent_z;
    _cent_x=(a->x + b->x + c->x + d->x)/4.0;
    _cent_y=(a->y + b->y + c->y + d->y)/4.0;
    _cent_z=(a->z + b->z + c->z + d->z)/4.0;

    return (Point3D(_cent_x, _cent_y, _cent_z));
};

geom::Sphere::Sphere(double x, double y, double z, double r) : Point3D(x,y,z) {
    this->r = r;
}

double geom::Sphere::radius() {
    return r;
}

geom::triangle::triangle(){}

geom::triangle::triangle (geom::Point3D* new_vertices){
	_no_of_vertices=3;
	assert(_no_of_vertices==3 && "For Triangle No. of Vertices should be 3!!");

	_no_of_sides=3;
	assert(_no_of_sides==3 && "For Triangle No. of Sides should be 3!!");

	_no_of_angles_with_circum_centre=3;
	assert(_no_of_angles_with_circum_centre==3 && "For Triangle No. of Angles of Two Vertices with Circumcentre should be 3!!");

	_no_of_vertex_vertex_angles=3;
	assert(_no_of_vertex_vertex_angles==3 && "For Triangle No. of Angles of Two Vertices with Another Vertex should be 3!!");

	for (int i=0; i<_no_of_vertices; i++){
		_vertices[i]=Point3D(new_vertices[i].X(), new_vertices[i].Y(), new_vertices[i].Z());
	}

	compute_circum_centre_radius();
	measure_vertex_to_vertex_distances();
	measure_vertex_vertex_angles();
}

void geom::triangle::compute_circum_centre_radius(){
	Point3D centre;
	double a=_vertices[1].dist(&_vertices[2]);
	double b=_vertices[2].dist(&_vertices[0]);
	double c=_vertices[0].dist(&_vertices[1]);
	double alpha=a*a*(b*b+c*c-a*a);
	double beta=b*b*(c*c+a*a-b*b);
	double gama=c*c*(a*a+b*b-c*c);
	double sum=alpha+beta+gama;

	double circum_x=(alpha*_vertices[0].X()+beta*_vertices[1].X()+gama*_vertices[2].X())/sum;
	double circum_y=(alpha*_vertices[0].Y()+beta*_vertices[1].Y()+gama*_vertices[2].Y())/sum;
	double circum_z=(alpha*_vertices[0].Z()+beta*_vertices[1].Z()+gama*_vertices[2].Z())/sum;

	_circum_centre=Point3D(circum_x, circum_y, circum_z);
	_circum_radius=(a*b*c)/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));

	cout<<"\n  Circumcentre: ";
	_circum_centre.print();
	cout<<"  Circumradius: "<<_circum_radius<<endl;
}

void geom::triangle::measure_circum_centre_to_vertex_distances(){
	for (int i=0; i<_no_of_vertices; i++){
		_circum_centre_vertex_distances[i][0]=i;
		_circum_centre_vertex_distances[i][1]=_circum_centre.dist(&_vertices[i]);
	}

	cout<<"\n\nCC-Vertex Distances:\n";
	for (int i=0; i<_no_of_vertices; i++){
		cout<<_circum_centre_vertex_distances[i][0]<<"  "<<_circum_centre_vertex_distances[i][1]<<endl;
	}
}

void geom::triangle::compute_centre_of_mass(){
	double x=((_vertices[0].X()+_vertices[1].X()+_vertices[2].X())/3.0);
	double y=((_vertices[0].Y()+_vertices[1].Y()+_vertices[2].Y())/3.0);
	double z=((_vertices[0].Z()+_vertices[1].Z()+_vertices[2].Z())/3.0);
	_centre_of_mass=Point3D(x, y, z);
}

void geom::triangle::measure_vertex_to_vertex_distances(){
	_vertex_vertex_distances[0][0]=0;
	_vertex_vertex_distances[0][1]=1;
	_vertex_vertex_distances[0][2]=_vertices[0].dist(&_vertices[1]);

	_vertex_vertex_distances[1][0]=1;
	_vertex_vertex_distances[1][1]=2;
	_vertex_vertex_distances[1][2]=_vertices[1].dist(&_vertices[2]);

	_vertex_vertex_distances[2][0]=2;
	_vertex_vertex_distances[2][1]=0;
	_vertex_vertex_distances[2][2]=_vertices[2].dist(&_vertices[0]);

	cout<<"\n\nVertex-Vertex Distances:\n";
	for (int i=0; i<_no_of_sides; i++){
		cout<<_vertex_vertex_distances[i][0]<<"  "<<_vertex_vertex_distances[i][1]<<"  "<<_vertex_vertex_distances[i][2]<<endl;
	}
}

void geom::triangle::measure_angles_vertices_with_circum_centre(){
	_circum_centre_vertex_angles_tri[0][0]=0;
	_circum_centre_vertex_angles_tri[0][1]=1;
	_circum_centre_vertex_angles_tri[0][2]=Point3D::angle_deg_triangle_law(&_vertices[0], &_circum_centre, &_vertices[1]);

	_circum_centre_vertex_angles_tri[1][0]=1;
	_circum_centre_vertex_angles_tri[1][1]=2;
	_circum_centre_vertex_angles_tri[1][2]=Point3D::angle_deg_triangle_law(&_vertices[1], &_circum_centre, &_vertices[2]);

	_circum_centre_vertex_angles_tri[2][0]=2;
	_circum_centre_vertex_angles_tri[2][1]=0;
	_circum_centre_vertex_angles_tri[2][2]=Point3D::angle_deg_triangle_law(&_vertices[2], &_circum_centre, &_vertices[0]);

	cout<<"\n\nCC-Vertex Angles:\n";
	for (int i=0; i<_no_of_angles_with_circum_centre; i++){
		cout<<_circum_centre_vertex_angles_tri[i][0]<<"  "<<_circum_centre_vertex_angles_tri[i][1]<<"  "<<_circum_centre_vertex_angles_tri[i][2]<<endl;
	}
}

void geom::triangle::measure_vertex_vertex_angles(){
	_vertex_vertex_angles_tri[0][0]=1;
	_vertex_vertex_angles_tri[0][1]=0;
	_vertex_vertex_angles_tri[0][2]=2;
	_vertex_vertex_angles_tri[0][3]=Point3D::angle_deg_triangle_law(&_vertices[1], &_vertices[0], &_vertices[2]);

	_vertex_vertex_angles_tri[1][0]=0;
	_vertex_vertex_angles_tri[1][1]=1;
	_vertex_vertex_angles_tri[1][2]=2;
	_vertex_vertex_angles_tri[1][3]=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[1], &_vertices[2]);

	_vertex_vertex_angles_tri[2][0]=1;
	_vertex_vertex_angles_tri[2][1]=2;
	_vertex_vertex_angles_tri[2][2]=0;
	_vertex_vertex_angles_tri[2][3]=Point3D::angle_deg_triangle_law(&_vertices[1], &_vertices[2], &_vertices[0]);

	cout<<"\n\nVertex-Vertex Angles:\n";
	for (int i=0; i<_no_of_vertex_vertex_angles; i++){
		cout<<_vertex_vertex_angles_tri[i][0]<<"  "<<_vertex_vertex_angles_tri[i][1]<<"  "<<_vertex_vertex_angles_tri[i][2]<<"  "<<_vertex_vertex_angles_tri[i][3]<<endl;
	}
}

void geom::triangle::set_vertex_order(int p, int q, int r){
	int i;
	_vertex_order[0]=p;
	_vertex_order[1]=q;
	_vertex_order[2]=r;
}

void geom::triangle::find_vertex_order(){
	double dist1=_vertices[0].dist(&_vertices[1]);
	double dist2=_vertices[1].dist(&_vertices[2]);
	double dist3=_vertices[2].dist(&_vertices[0]);
cout<<"\n\nDIST Vertices: "<<dist1<<"  "<<dist2<<"  "<<dist3<<endl;

	if ((dist1>=dist2) && (dist1>=dist3)){
		if (_vertices[0].X()<=_vertices[1].X()){
			_vertex_order[0]=0;
			_vertex_order[1]=1;
		}else{
			_vertex_order[0]=1;
			_vertex_order[1]=0;
		}

		_vertex_order[2]=2;
	}

	if ((dist2>=dist1) && (dist2>=dist3)){
		if (_vertices[1].X()<=_vertices[2].X()){
			_vertex_order[0]=1;
			_vertex_order[1]=2;
		}else{
			_vertex_order[0]=2;
			_vertex_order[1]=1;
		}

		_vertex_order[2]=0;
	}

	if ((dist3>=dist1) && (dist3>=dist2)){
		if (_vertices[2].X()<=_vertices[0].X()){
			_vertex_order[0]=2;
			_vertex_order[1]=0;
		}else{
			_vertex_order[0]=0;
			_vertex_order[1]=2;
		}

		_vertex_order[2]=1;
	}
}

double (*geom::triangle::measure_direction_cosines_of_sides())[3]{
	double x=_vertices[1].X() - _vertices[0].X();
	double y=_vertices[1].Y() - _vertices[0].Y();
	double z=_vertices[1].Z() - _vertices[0].Z();
	double r=sqrt(x*x+y*y+z*z);

	double l=x/r;
	double m=y/r;
	double n=z/r;

	_direction_cosines[0][0]=l;
	_direction_cosines[0][1]=m;
	_direction_cosines[0][2]=n;
	cout<<"\nDC 0: "<<l<<", "<<m<<", "<<n;
	x=_vertices[2].X() - _vertices[1].X();
	y=_vertices[2].Y() - _vertices[1].Y();
	z=_vertices[2].Z() - _vertices[1].Z();
	r=sqrt(x*x+y*y+z*z);

	l=x/r;
	m=y/r;
	n=z/r;

	_direction_cosines[1][0]=l;
	_direction_cosines[1][1]=m;
	_direction_cosines[1][2]=n;
	cout<<"\nDC 1: "<<l<<", "<<m<<", "<<n;
	x=_vertices[0].X() - _vertices[2].X();
	y=_vertices[0].Y() - _vertices[2].Y();
	z=_vertices[0].Z() - _vertices[2].Z();
	r=sqrt(x*x+y*y+z*z);

	l=x/r;
	m=y/r;
	n=z/r;

	_direction_cosines[2][0]=l;
	_direction_cosines[2][1]=m;
	_direction_cosines[2][2]=n;
	cout<<"\nDC 2: "<<l<<", "<<m<<", "<<n;
	return _direction_cosines;
}

/*void geom::triangle::calc_projections_on_standard_triangle(geom::triangle & ctdrn, double* proj){
	if (_vertices[0].dist(&ctdrn.get_vertices()[0]) == _vertices[1].dist(&ctdrn.get_vertices()[1])){
		proj[0]=ctdrn.get_vertices()[0].dist(&ctdrn.get_vertices()[1]);
	}else{
		double prx=	_direction_cosines[0][0]*(ctdrn.get_vertices()[1].X()-ctdrn.get_vertices()[0].X());
		double pry= _direction_cosines[0][1]*(ctdrn.get_vertices()[1].Y()-ctdrn.get_vertices()[0].Y());
		double prz=	_direction_cosines[0][2]*(ctdrn.get_vertices()[1].Z()-ctdrn.get_vertices()[0].Z());

		proj[0]=(prx+pry+prz);
	}

	if (_vertices[1].dist(&ctdrn.get_vertices()[1]) == _vertices[2].dist(&ctdrn.get_vertices()[2])){
		proj[1]=ctdrn.get_vertices()[1].dist(&ctdrn.get_vertices()[2]);
	}else{
		double prx= _direction_cosines[1][0]*(ctdrn.get_vertices()[2].X()-ctdrn.get_vertices()[1].X());
		double pry= _direction_cosines[1][1]*(ctdrn.get_vertices()[2].Y()-ctdrn.get_vertices()[1].Y());
		double prz= _direction_cosines[1][2]*(ctdrn.get_vertices()[2].Z()-ctdrn.get_vertices()[1].Z());

		proj[1]=(prx+pry+prz);
	}

	if (_vertices[0].dist(&ctdrn.get_vertices()[0]) == _vertices[2].dist(&ctdrn.get_vertices()[2])){
		proj[2]=ctdrn.get_vertices()[0].dist(&ctdrn.get_vertices()[2]);
	}else{
		double prx= _direction_cosines[2][0]*(ctdrn.get_vertices()[0].X()-ctdrn.get_vertices()[2].X());
		double pry= _direction_cosines[2][1]*(ctdrn.get_vertices()[0].Y()-ctdrn.get_vertices()[2].Y());
		double prz= _direction_cosines[2][2]*(ctdrn.get_vertices()[0].Z()-ctdrn.get_vertices()[2].Z());

		proj[2]=(prx+pry+prz);
	}

	cout<<"\n Projections on: "<<endl;

	for (int i=0; i<_no_of_sides; i++){
		cout<<"Side "<<(i+1)<<": "<<proj[i]<<endl;
	}
}

void geom::triangle::deviation_distances(geom::triangle ctdrn, double* deviate_dist){
	double* proj=new double[_no_of_sides];
	calc_projections_on_standard_triangle(ctdrn, proj);

	cout<<"\nVert-Vert Side Deviations: \n";
	for (int i=0; i<_no_of_sides; i++){
		double diff=_vertex_vertex_distances[i][2]-proj[i];
		deviate_dist[i]=diff;
		cout<<"Side "<<(i+1)<<": "<<deviate_dist[i]<<endl;
	}

	delete[] proj;
}*/

void geom::triangle::calc_projections_on_standard_triangle(geom::triangle & ctdrn, double dc[][3], double* proj){
	if (_vertices[0].dist(&ctdrn.get_vertices()[0])==0 && _vertices[1].dist(&ctdrn.get_vertices()[1])==0){
		proj[0]=ctdrn.get_vertices()[0].dist(&ctdrn.get_vertices()[1]);
	}else{
		double prx=	dc[0][0]*(ctdrn.get_vertices()[1].X()-ctdrn.get_vertices()[0].X());
		double pry= dc[0][1]*(ctdrn.get_vertices()[1].Y()-ctdrn.get_vertices()[0].Y());
		double prz=	dc[0][2]*(ctdrn.get_vertices()[1].Z()-ctdrn.get_vertices()[0].Z());

		proj[0]=(prx+pry+prz);
	}

	if (_vertices[1].dist(&ctdrn.get_vertices()[1])==0 && _vertices[2].dist(&ctdrn.get_vertices()[2])==0){
		proj[1]=ctdrn.get_vertices()[1].dist(&ctdrn.get_vertices()[2]);
	}else{
		double prx= dc[1][0]*(ctdrn.get_vertices()[2].X()-ctdrn.get_vertices()[1].X());
		double pry= dc[1][1]*(ctdrn.get_vertices()[2].Y()-ctdrn.get_vertices()[1].Y());
		double prz= dc[1][2]*(ctdrn.get_vertices()[2].Z()-ctdrn.get_vertices()[1].Z());

		proj[1]=(prx+pry+prz);
	}

	if (_vertices[0].dist(&ctdrn.get_vertices()[0])==0 && _vertices[2].dist(&ctdrn.get_vertices()[2])==0){
		proj[2]=ctdrn.get_vertices()[0].dist(&ctdrn.get_vertices()[2]);
	}else{
		double prx= dc[2][0]*(ctdrn.get_vertices()[0].X()-ctdrn.get_vertices()[2].X());
		double pry= dc[2][1]*(ctdrn.get_vertices()[0].Y()-ctdrn.get_vertices()[2].Y());
		double prz= dc[2][2]*(ctdrn.get_vertices()[0].Z()-ctdrn.get_vertices()[2].Z());

		proj[2]=(prx+pry+prz);
	}

	cout<<"\n Projections on: "<<endl;

	for (int i=0; i<_no_of_sides; i++){
		cout<<"Side "<<(i+1)<<": "<<proj[i]<<endl;
	}
}

void geom::triangle::deviation_distances(geom::triangle ctdrn, double dc[][3], double* deviate_dist){
	double* proj=new double[_no_of_sides];
	calc_projections_on_standard_triangle(ctdrn, dc, proj);

	cout<<"\nVert-Vert Side Deviations: \n";
	for (int i=0; i<_no_of_sides; i++){
		double diff=_vertex_vertex_distances[i][2]-proj[i];
		deviate_dist[i]=diff;
		cout<<"Side "<<(i+1)<<": "<<deviate_dist[i]<<endl;
	}

	delete[] proj;
}

void geom::triangle::deviation_angles_with_circum_centre(geom::triangle ctdrn, double* deviate_angle){
	cout<<"\nAngle Deviations CC:\n";
	for (int i=0; i<_no_of_angles_with_circum_centre; i++){
		double diff=get_angles_vertices_with_circum_centre()[i][2] - ctdrn.get_angles_vertices_with_circum_centre()[i][2];
		deviate_angle[i]=diff;
		cout<<"CC-Vert Angle "<<(i+1)<<": "<<deviate_angle[i]<<endl;
	}
}

void geom::triangle::deviation_angles_vertex_vertex(triangle comp, double* deviate_angle){
	cout<<"\nAngle Deviations Vert-Vert:\n";
	for (int i=0; i<_no_of_vertex_vertex_angles; i++){
		double diff=get_vertex_vertex_angles()[i][3] - comp.get_vertex_vertex_angles()[i][3];
		deviate_angle[i]=diff;
		cout<<"Vert-Vert Angle "<<(i+1)<<": "<<deviate_angle[i]<<"  Stand:"<<get_vertex_vertex_angles()[i][3]<<"  Comp:"<<comp.get_vertex_vertex_angles()[i][3]<<endl;
	}
}

geom::Point3D* geom::triangle::get_vertices(){
	return _vertices;
}

int* geom::triangle::get_vertex_order(){
	return _vertex_order;
}

geom::Point3D geom::triangle::get_circum_centre(){
	return _circum_centre;
}

double geom::triangle::get_circum_radius(){
	return _circum_radius;
}

geom::Point3D geom::triangle::get_centre_of_mass(){
	return _centre_of_mass;
}

double (*geom::triangle::get_circum_centre_to_vertex_distances())[2]{
	return _circum_centre_vertex_distances;
}

double (*geom::triangle::get_vertex_vertex_distances())[3]{
	return _vertex_vertex_distances;
}

double (*geom::triangle::get_angles_vertices_with_circum_centre())[3]{
	return _circum_centre_vertex_angles_tri;
}

double (*geom::triangle::get_vertex_vertex_angles())[4]{
	return _vertex_vertex_angles_tri;
}

double (*geom::triangle::get_direction_cosines())[3]{
	return _direction_cosines;
}

geom::tetragonalplanar::tetragonalplanar(){}

geom::tetragonalplanar::tetragonalplanar (geom::Point3D* new_vertices){
	_no_of_vertices=4;
	assert(_no_of_vertices==4 && "For Tetragonal Planar No. of Vertices should be 4!!");

	_no_of_sides=4;
	assert(_no_of_sides==4 && "For Tetragonal Planar No. of Sides should be 4!!");

	_no_of_vertex_vertex_angles=4;
	assert(_no_of_vertex_vertex_angles==4 && "For Tetragonal Planar No. of Angles of Two Vertices with Another Vertex should be 4!!");

	//_vertices=new Point3D[_no_of_vertices];
	for (int i=0; i<_no_of_vertices; i++){
		_vertices[i]=Point3D(new_vertices[i].X(), new_vertices[i].Y(), new_vertices[i].Z());
	}

	//_vertex_order=new int[_no_of_vertices];

	compute_vertex_order();
	measure_vertex_vertex_distances();
	measure_vertex_vertex_angles();
}

void geom::tetragonalplanar::measure_vertex_vertex_distances(){
	_vertex_vertex_distances[0][0]=_vertex_order[0];
	_vertex_vertex_distances[0][1]=_vertex_order[1];
	_vertex_vertex_distances[0][2]=_vertices[_vertex_order[0]].dist(&_vertices[_vertex_order[1]]);

	_vertex_vertex_distances[1][0]=_vertex_order[1];
	_vertex_vertex_distances[1][1]=_vertex_order[2];
	_vertex_vertex_distances[1][2]=_vertices[_vertex_order[1]].dist(&_vertices[_vertex_order[2]]);

	_vertex_vertex_distances[2][0]=_vertex_order[2];
	_vertex_vertex_distances[2][1]=_vertex_order[3];
	_vertex_vertex_distances[2][2]=_vertices[_vertex_order[2]].dist(&_vertices[_vertex_order[3]]);

	_vertex_vertex_distances[3][0]=_vertex_order[3];
	_vertex_vertex_distances[3][1]=_vertex_order[0];
	_vertex_vertex_distances[3][2]=_vertices[_vertex_order[3]].dist(&_vertices[_vertex_order[0]]);

	cout<<"\n\nVertex-Vertex Distances:\n";
	for (int i=0; i<_no_of_sides; i++){
		cout<<_vertex_vertex_distances[i][0]<<"  "<<_vertex_vertex_distances[i][1]<<"  "<<_vertex_vertex_distances[i][2]<<endl;
	}
}

void geom::tetragonalplanar::measure_vertex_vertex_angles(){
	_vertex_vertex_angles_tri[0][0]=_vertex_order[1];
	_vertex_vertex_angles_tri[0][1]=_vertex_order[0];
	_vertex_vertex_angles_tri[0][2]=_vertex_order[3];
	_vertex_vertex_angles_tri[0][3]=Point3D::angle_deg_triangle_law(&_vertices[_vertex_order[1]], &_vertices[_vertex_order[0]], &_vertices[_vertex_order[3]]);

	_vertex_vertex_angles_tri[1][0]=_vertex_order[0];
	_vertex_vertex_angles_tri[1][1]=_vertex_order[1];
	_vertex_vertex_angles_tri[1][2]=_vertex_order[2];
	_vertex_vertex_angles_tri[1][3]=Point3D::angle_deg_triangle_law(&_vertices[_vertex_order[0]], &_vertices[_vertex_order[1]], &_vertices[_vertex_order[2]]);

	_vertex_vertex_angles_tri[2][0]=_vertex_order[1];
	_vertex_vertex_angles_tri[2][1]=_vertex_order[2];
	_vertex_vertex_angles_tri[2][2]=_vertex_order[3];
	_vertex_vertex_angles_tri[2][3]=Point3D::angle_deg_triangle_law(&_vertices[_vertex_order[1]], &_vertices[_vertex_order[2]], &_vertices[_vertex_order[3]]);

	_vertex_vertex_angles_tri[3][0]=_vertex_order[0];
	_vertex_vertex_angles_tri[3][1]=_vertex_order[3];
	_vertex_vertex_angles_tri[3][2]=_vertex_order[2];
	_vertex_vertex_angles_tri[3][3]=Point3D::angle_deg_triangle_law(&_vertices[_vertex_order[0]], &_vertices[_vertex_order[3]], &_vertices[_vertex_order[2]]);

	cout<<"\n\nVertex-Vertex Angles:\n";
	for (int i=0; i<_no_of_vertex_vertex_angles; i++){
		cout<<_vertex_vertex_angles_tri[i][0]<<"  "<<_vertex_vertex_angles_tri[i][1]<<"  "<<_vertex_vertex_angles_tri[i][2]<<"  "<<_vertex_vertex_angles_tri[i][3]<<endl;
	}
}

void geom::tetragonalplanar::compute_vertex_order(){
	double diff1[6], diff2[6];
	int v_order[6][4];
	double ang1=Point3D::angle_deg_triangle_law(&_vertices[3], &_vertices[0], &_vertices[1]);
	double ang2=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[1], &_vertices[2]);
	double ang3=Point3D::angle_deg_triangle_law(&_vertices[1], &_vertices[2], &_vertices[3]);
	double ang4=Point3D::angle_deg_triangle_law(&_vertices[2], &_vertices[3], &_vertices[0]);
	diff1[0]=180.0-(ang1+ang3);
	diff2[0]=180.0-(ang2+ang4);
	cout<<"\n 1st: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[0]<<", "<<diff2[0];

	v_order[0][0]=0;
	v_order[0][1]=1;
	v_order[0][2]=2;
	v_order[0][3]=3;

	ang1=Point3D::angle_deg_triangle_law(&_vertices[2], &_vertices[0], &_vertices[1]);
	ang2=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[1], &_vertices[3]);
	ang3=Point3D::angle_deg_triangle_law(&_vertices[1], &_vertices[3], &_vertices[2]);
	ang4=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[2], &_vertices[3]);
	diff1[1]=180.0-(ang1+ang3);
	diff2[1]=180.0-(ang2+ang4);
	cout<<"\n 2nd: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[1]<<", "<<diff2[1];

	v_order[1][0]=0;
	v_order[1][1]=1;
	v_order[1][2]=3;
	v_order[1][3]=2;

	ang1=Point3D::angle_deg_triangle_law(&_vertices[2], &_vertices[0], &_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[2], &_vertices[1]);
	ang3=Point3D::angle_deg_triangle_law(&_vertices[2], &_vertices[1], &_vertices[3]);
	ang4=Point3D::angle_deg_triangle_law(&_vertices[1], &_vertices[3], &_vertices[0]);
	diff1[2]=180.0-(ang1+ang3);
	diff2[2]=180.0-(ang2+ang4);
	cout<<"\n 3rd: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[2]<<", "<<diff2[2];

	v_order[2][0]=0;
	v_order[2][1]=2;
	v_order[2][2]=1;
	v_order[2][3]=3;

	ang1=Point3D::angle_deg_triangle_law(&_vertices[1], &_vertices[0], &_vertices[2]);
	ang2=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[2], &_vertices[3]);
	ang3=Point3D::angle_deg_triangle_law(&_vertices[2], &_vertices[3], &_vertices[1]);
	ang4=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[1], &_vertices[3]);
	diff1[3]=180.0-(ang1+ang3);
	diff2[3]=180.0-(ang2+ang4);
	cout<<"\n 4th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[3]<<", "<<diff2[3];

	v_order[3][0]=0;
	v_order[3][1]=2;
	v_order[3][2]=3;
	v_order[3][3]=1;

	ang1=Point3D::angle_deg_triangle_law(&_vertices[2], &_vertices[0], &_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[3], &_vertices[1]);
	ang3=Point3D::angle_deg_triangle_law(&_vertices[3], &_vertices[1], &_vertices[2]);
	ang4=Point3D::angle_deg_triangle_law(&_vertices[1], &_vertices[2], &_vertices[0]);
	diff1[4]=180.0-(ang1+ang3);
	diff2[4]=180.0-(ang2+ang4);
	cout<<"\n 5th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[4]<<", "<<diff2[4];

	v_order[4][0]=0;
	v_order[4][1]=3;
	v_order[4][2]=1;
	v_order[4][3]=2;

	ang1=Point3D::angle_deg_triangle_law(&_vertices[1], &_vertices[0], &_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[3], &_vertices[2]);
	ang3=Point3D::angle_deg_triangle_law(&_vertices[3], &_vertices[2], &_vertices[1]);
	ang4=Point3D::angle_deg_triangle_law(&_vertices[0], &_vertices[1], &_vertices[2]);
	diff1[5]=180.0-(ang1+ang3);
	diff2[5]=180.0-(ang2+ang4);
	cout<<"\n 6th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[5]<<", "<<diff2[5];

	v_order[5][0]=0;
	v_order[5][1]=3;
	v_order[5][2]=2;
	v_order[5][3]=1;

	double diff1_min=abs(diff1[0]), diff2_min=abs(diff2[0]);
	int min_index1=0, min_index2=0;
	for (int i=1; i<6; i++){
		if (abs(diff1[i])<abs(diff1_min)){
			diff1_min=abs(diff1[i]);
			min_index1=i;
		}

		if (abs(diff2[i])<abs(diff2_min)){
			diff2_min=abs(diff2[i]);
			min_index2=i;
		}
	}

	if (abs(diff1_min)<=abs(diff2_min)){
		_vertex_order[0]=v_order[min_index1][0];
		_vertex_order[1]=v_order[min_index1][1];
		_vertex_order[2]=v_order[min_index1][2];
		_vertex_order[3]=v_order[min_index1][3];
		cout<<"\nSelected Order Index: "<<min_index1<<"; Diff Min: "<<diff1_min;
	}else{
		_vertex_order[0]=v_order[min_index2][0];
		_vertex_order[1]=v_order[min_index2][1];
		_vertex_order[2]=v_order[min_index2][2];
		_vertex_order[3]=v_order[min_index2][3];
		cout<<"\nSelected Order Index: "<<min_index2<<"; Diff Min: "<<diff2_min;
	}
}

geom::Point3D* geom::tetragonalplanar::get_vertices(){
	return _vertices;
}

int* geom::tetragonalplanar::get_vertex_order(){
	return _vertex_order;
}

double (*geom::tetragonalplanar::get_vertex_vertex_distances())[3]{
	return _vertex_vertex_distances;
}

double (*geom::tetragonalplanar::get_vertex_vertex_angles())[4]{
	return _vertex_vertex_angles_tri;
}

geom::tetrahedron::tetrahedron(){}

geom::tetrahedron::tetrahedron (geom::Point3D new_vertices[]){
	_no_of_vertices=4;
	assert(_no_of_vertices==4 && "For Tetrahedron No. of Vertices should be 4!!");

	_no_of_sides=6;
	assert(_no_of_sides==6 && "For Tetrahedron No. of Sides should be 6!!");

	_no_of_faces=4;
	assert(_no_of_faces==4 && "For Tetrahedron No. of Faces should be 4!!");

	_no_of_angles_with_circum_centre=6;
	assert(_no_of_angles_with_circum_centre==6 && "For Tetrahedron No. of Angles of Two Vertices with Circumcentre should be 6!!");

	_no_of_vertex_vertex_angles=12;
	assert(_no_of_vertex_vertex_angles==12 && "For Tetrahedron No. of Angles of Two Vertices with Another Vertex should be 12!!");

	for (int i=0; i<_no_of_vertices; i++){
		_vertices[i]=Point3D(new_vertices[i].X(), new_vertices[i].Y(), new_vertices[i].Z());
		cout<<"\nVERTICES: "<<_vertices[i].X()<<", "<<_vertices[i].Y()<<", "<<_vertices[i].Z();
	}

	compute_circum_centre_radius();
	measure_vertex_to_vertex_distances();
	measure_vertex_vertex_angles();
}

void geom::tetrahedron::set_vertex_order(int vertex_order[], int v){
	int i;
	for (i=0; i<3; i++){
		_vertex_order[i]=vertex_order[i];
	}

	_vertex_order[i]=v;
}

void geom::tetrahedron::set_vertex_order(int v1, int v2, int v3, int v4){
	_vertex_order[0]=v1;
	_vertex_order[1]=v2;
	_vertex_order[2]=v3;
	_vertex_order[3]=v4;
}

void geom::tetrahedron::compute_circum_centre_radius(){
	Point3D centre;
	Matrix<float, 4, 4> alpha;
	Matrix<float, 4, 4> gamma;
	Matrix<float, 4, 4> mx;
	Matrix<float, 4, 4> my;
	Matrix<float, 4, 4> mz;
	double aa, ga, dx, dy, dz;
	double cx, cy, cz;
	double sum[_no_of_vertices];
	double radius;

	for (int i=0; i<_no_of_vertices; i++){
		sum[i]=_vertices[i].X()*_vertices[i].X()+_vertices[i].Y()*_vertices[i].Y()+_vertices[i].Z()*_vertices[i].Z();
	}

	alpha<<_vertices[0].X(), _vertices[0].Y(), _vertices[0].Z(), 1,
			_vertices[1].X(), _vertices[1].Y(), _vertices[1].Z(), 1,
			_vertices[2].X(), _vertices[2].Y(), _vertices[2].Z(), 1,
			_vertices[3].X(), _vertices[3].Y(), _vertices[3].Z(), 1;

	gamma<<sum[0], _vertices[0].X(), _vertices[0].Y(), _vertices[0].Z(),
			sum[1], _vertices[1].X(), _vertices[1].Y(), _vertices[1].Z(),
			sum[2], _vertices[2].X(), _vertices[2].Y(), _vertices[2].Z(),
			sum[3], _vertices[3].X(), _vertices[3].Y(), _vertices[3].Z();

	mx<<sum[0], _vertices[0].Y(), _vertices[0].Z(), 1,
		sum[1], _vertices[1].Y(), _vertices[1].Z(), 1,
		sum[2], _vertices[2].Y(), _vertices[2].Z(), 1,
		sum[3], _vertices[3].Y(), _vertices[3].Z(), 1;

	my<<sum[0], _vertices[0].X(), _vertices[0].Z(), 1,
		sum[1], _vertices[1].X(), _vertices[1].Z(), 1,
		sum[2], _vertices[2].X(), _vertices[2].Z(), 1,
		sum[3], _vertices[3].X(), _vertices[3].Z(), 1;

	mz<<sum[0], _vertices[0].X(), _vertices[0].Y(), 1,
		sum[1], _vertices[1].X(), _vertices[1].Y(), 1,
		sum[2], _vertices[2].X(), _vertices[2].Y(), 1,
		sum[3], _vertices[3].X(), _vertices[3].Y(), 1;

	aa=alpha.determinant();
	ga=gamma.determinant();
	dx=mx.determinant();
	dy=-my.determinant();
	dz=mz.determinant();

	cx=dx/(2*aa);
	cy=dy/(2*aa);
	cz=dz/(2*aa);

	_circum_centre=Point3D(cx, cy, cz);
	_circum_radius=sqrt(dx*dx+dy*dy+dz*dz-4*aa*ga)/(2*abs(aa));

	cout<<"\n  Circumcentre: ";
	_circum_centre.print();
	cout<<"  Circumradius: "<<_circum_radius<<endl;
}

void geom::tetrahedron::measure_circum_centre_to_vertex_distances(){
	for (int i=0; i<_no_of_vertices; i++){
		_circum_centre_vertex_distances[i][0]=i;
		_circum_centre_vertex_distances[i][1]=_circum_centre.dist(&_vertices[i]);
	}

	cout<<"\n\nCC-Vertex Distances:\n";
	for (int i=0; i<_no_of_vertices; i++){
		cout<<_circum_centre_vertex_distances[i][0]<<"  "<<_circum_centre_vertex_distances[i][1]<<endl;
	}
}

void geom::tetrahedron::compute_centre_of_mass(){
	double x=((_vertices[0].X()+_vertices[1].X()+_vertices[2].X()+_vertices[3].X())/4.0);
	double y=((_vertices[0].Y()+_vertices[1].Y()+_vertices[2].Y()+_vertices[3].Y())/4.0);
	double z=((_vertices[0].Z()+_vertices[1].Z()+_vertices[2].Z()+_vertices[3].Z())/4.0);
	_centre_of_mass=Point3D(x, y, z);
}

void geom::tetrahedron::measure_vertex_to_vertex_distances(){
	int k=0;
	for (int i=0; i<_no_of_vertices-1; i++){
		for (int j=i+1; j<_no_of_vertices; j++){
			_vertex_vertex_distances[k][0]=i;
			_vertex_vertex_distances[k][1]=j;
			_vertex_vertex_distances[k][2]=_vertices[i].dist(&_vertices[j]);
			k++;
		}
	}

	cout<<"\n\nVertex-Vertex Distances:\n";
	for (int i=0; i<k; i++){
		cout<<_vertex_vertex_distances[i][0]<<"  "<<_vertex_vertex_distances[i][1]<<"  "<<_vertex_vertex_distances[i][2]<<endl;
	}
}

void geom::tetrahedron::measure_angles_vertices_with_circum_centre(){
	int k=0;
	for (int i=0; i<_no_of_vertices-1; i++){
		for (int j=i+1; j<_no_of_vertices; j++){
			_circum_centre_vertex_angles_tri[k][0]=i;
			_circum_centre_vertex_angles_tri[k][1]=j;
			double ang=Point3D::angle_deg_triangle_law(&_vertices[i], &_circum_centre, &_vertices[j]);
			_circum_centre_vertex_angles_tri[k][2]=ang;
			k++;
		}
	}

	cout<<"\n\nCC-Vertex Angles:\n";
	for (int i=0; i<_no_of_vertices; i++){
		cout<<_circum_centre_vertex_angles_tri[i][0]<<"  "<<_circum_centre_vertex_angles_tri[i][1]<<"  "<<_circum_centre_vertex_angles_tri[i][2]<<endl;
	}
}

void geom::tetrahedron::measure_vertex_vertex_angles(){
	int p=0;
	for (int i=0; i<_no_of_vertices; i++){
		for (int j=0; j<_no_of_vertices; j++){
			for (int k=0; k<_no_of_vertices; k++){
				if (i!=j && j<k && i!=k){
					_vertex_vertex_angles_tri[p][0]=j;
					_vertex_vertex_angles_tri[p][1]=i;
					_vertex_vertex_angles_tri[p][2]=k;
					double ang=Point3D::angle_deg_triangle_law(&_vertices[j], &_vertices[i], &_vertices[k]);
					_vertex_vertex_angles_tri[p][3]=ang;
					p++;
				}
			}
		}
	}

	cout<<"\n\nVertex-Vertex Angles:\n";
	for (int i=0; i<p; i++){
		cout<<_vertex_vertex_angles_tri[i][0]<<"  "<<_vertex_vertex_angles_tri[i][1]<<"  "<<_vertex_vertex_angles_tri[i][2]<<"  "<<_vertex_vertex_angles_tri[i][3]<<endl;
	}
}

void geom::tetrahedron::measure_surface_areas_and_order_of_vertices(){
	double a, b, c, s, ar;
	int p=0, base_index1, base_index2, base_index3;

	for (int i=0; i<_no_of_vertices-2; i++){
		for (int j=i+1; j<_no_of_vertices-1; j++){
			a=_vertices[i].dist(&_vertices[j]);

			for (int k=j+1; k<_no_of_vertices; k++){
				_surface_areas[p][0]=i;
				_surface_areas[p][1]=j;
				_surface_areas[p][2]=k;
				b=_vertices[j].dist(&_vertices[k]);
				c=_vertices[k].dist(&_vertices[i]);
				s=(a+b+c)/2.0;
				ar=sqrt(s*(s-a)*(s-b)*(s-c));
				_surface_areas[p][3]=ar;
				p++;

				if (ar>_max_surface_area){
					_max_surface_area=ar;
					base_index1=i;
					base_index2=j;
					base_index3=k;
				}
			}

		}
	}

	find_vertex_order(base_index1, base_index2, base_index3);
	cout<<"\n\nMAX SURFACE AREA: "<<_max_surface_area<<"  BASE INDEX: "<<base_index1<<"  "<<base_index2<<"  "<<base_index3;
}

void geom::tetrahedron::find_vertex_order(int _base_vertex1, int _base_vertex2, int _base_vertex3){
	double dist1=_vertices[_base_vertex1].dist(&_vertices[_base_vertex2]);
	double dist2=_vertices[_base_vertex2].dist(&_vertices[_base_vertex3]);
	double dist3=_vertices[_base_vertex1].dist(&_vertices[_base_vertex3]);

	//Order of vertices according to the sides of base with max. surface area
	if ((dist1>=dist2) && (dist1>=dist3)){
		if (_vertices[_base_vertex1].X()<=_vertices[_base_vertex2].X()){
			_vertex_order[0]=_base_vertex1;
			_vertex_order[1]=_base_vertex2;
		}else{
			_vertex_order[0]=_base_vertex2;
			_vertex_order[1]=_base_vertex1;
		}

		_vertex_order[2]=_base_vertex3;
	}

	if ((dist2>=dist1) && (dist2>=dist3)){
		if (_vertices[_base_vertex2].X()<=_vertices[_base_vertex3].X()){
			_vertex_order[0]=_base_vertex2;
			_vertex_order[1]=_base_vertex3;
		}else{
			_vertex_order[0]=_base_vertex3;
			_vertex_order[1]=_base_vertex2;
		}

		_vertex_order[2]=_base_vertex1;
	}

	if ((dist3>=dist1) && (dist3>=dist2)){
		if (_vertices[_base_vertex1].X()<=_vertices[_base_vertex3].X()){
			_vertex_order[0]=_base_vertex1;
			_vertex_order[1]=_base_vertex3;
		}else{
			_vertex_order[0]=_base_vertex3;
			_vertex_order[1]=_base_vertex1;
		}

		_vertex_order[2]=_base_vertex2;
	}

	//Check Fourth Vertex
	for (int i=0; i<_no_of_vertices; i++){
		if ((i!=_base_vertex1) && (i!=_base_vertex2) && (i!=_base_vertex3)){
			_vertex_order[3]=i;
			break;
		}
	}
}

double (*geom::tetrahedron::measure_direction_cosines_of_sides())[3]{
	double x, y, z;
	double l, m, n;
	double r;
	int k=0;

	for (int i=0; i<_no_of_vertices-1; i++){
		for (int j=i+1; j<_no_of_vertices; j++){
			x=_vertices[j].X() - _vertices[i].X();
			y=_vertices[j].Y() - _vertices[i].Y();
			z=_vertices[j].Z() - _vertices[i].Z();
			r=sqrt(x*x+y*y+z*z);

			l=x/r;
			m=y/r;
			n=z/r;

			_direction_cosines[k][0]=l;
			_direction_cosines[k][1]=m;
			_direction_cosines[k][2]=n;
			k++;
		}
	}

	cout<<"\n Direction Cosine: "<<endl;

	for (int i=0; i<_no_of_sides; i++){
		cout<<"Side "<<(i+1)<<" DC: "<<_direction_cosines[i][0]<<", "<<_direction_cosines[i][1]<<", "<<_direction_cosines[i][2]<<endl;
	}

	return _direction_cosines;
}

void geom::tetrahedron::calc_projections_on_standard_tetra(geom::tetrahedron & ctdrn, double proj[]){
	int k=0;

	for (int i=0; i<(_no_of_vertices-1); i++){
		for (int j=i+1; j<_no_of_vertices; j++){
			if (_vertices[i].dist(&ctdrn.get_vertices()[i])==0 && _vertices[j].dist(&ctdrn.get_vertices()[j])==0){
				proj[k]=ctdrn.get_vertices()[i].dist(&ctdrn.get_vertices()[j]);
			}else{
				double prx=	_direction_cosines[k][0]*(ctdrn._vertices[j].X()-ctdrn._vertices[i].X());
				double pry= _direction_cosines[k][1]*(ctdrn._vertices[j].Y()-ctdrn._vertices[i].Y());
				double prz=	_direction_cosines[k][2]*(ctdrn._vertices[j].Z()-ctdrn._vertices[i].Z());

				proj[k]=(prx+pry+prz);
			}

			k++;
		}
	}

	cout<<"\n Projections on: "<<endl;

	for (int i=0; i<_no_of_sides; i++){
		cout<<"Side "<<(i+1)<<" Proj: "<<proj[i]<<endl;
	}
}

void geom::tetrahedron::deviation_distances(geom::tetrahedron ctdrn, double deviate_dist[]){
	double proj[_no_of_sides];
	calc_projections_on_standard_tetra(ctdrn, proj);

	cout<<"\nVert-Vert Side Deviations: \n";
	for (int i=0; i<_no_of_sides; i++){
		double diff=_vertex_vertex_distances[i][2]-proj[i];
		deviate_dist[i]=diff;
		cout<<"Side "<<(i+1)<<": "<<deviate_dist[i]<<endl;
	}
}

void geom::tetrahedron::calc_projections_on_standard_tetra(geom::tetrahedron & ctdrn, double dc[][3], double proj[]){
	int k=0;

	for (int i=0; i<(_no_of_vertices-1); i++){
		for (int j=(i+1); j<_no_of_vertices; j++){
			cout<<"\nProj Vertices 1st End: ("<<_vertices[i].X()<<", "<<_vertices[i].Y()<<", "<<_vertices[i].Z()<<"); ("<<ctdrn.get_vertices()[i].X()<<", "<<ctdrn.get_vertices()[i].Y()<<", "<<ctdrn.get_vertices()[i].Z()<<")";
			cout<<"\nProj Vertices 2nd End: ("<<_vertices[j].X()<<", "<<_vertices[j].Y()<<", "<<_vertices[j].Z()<<"); ("<<ctdrn.get_vertices()[j].X()<<", "<<ctdrn.get_vertices()[j].Y()<<", "<<ctdrn.get_vertices()[j].Z()<<")";
			cout<<"\nDist Vert Proj: "<<_vertices[i].dist(&ctdrn.get_vertices()[i])<<", "<<_vertices[j].dist(&ctdrn.get_vertices()[j]);
			if (_vertices[i].dist(&ctdrn.get_vertices()[i])==0 && _vertices[j].dist(&ctdrn.get_vertices()[j])==0){
				proj[k]=ctdrn.get_vertices()[i].dist(&ctdrn.get_vertices()[j]);
				cout<<"\nProj 1st Tetra";
			}else{
				double prx=	dc[k][0]*(ctdrn._vertices[j].X()-ctdrn._vertices[i].X());
				double pry= dc[k][1]*(ctdrn._vertices[j].Y()-ctdrn._vertices[i].Y());
				double prz=	dc[k][2]*(ctdrn._vertices[j].Z()-ctdrn._vertices[i].Z());

				proj[k]=(prx+pry+prz);
				cout<<"\nProj 2nd Tetra";
			}

			k++;
		}
	}

	cout<<"\n Projections on: "<<endl;

	for (int i=0; i<_no_of_sides; i++){
		cout<<"Side "<<(i+1)<<" Proj: "<<proj[i]<<endl;
	}
}

void geom::tetrahedron::deviation_distances(geom::tetrahedron ctdrn, double dc[][3], double deviate_dist[]){
	double proj[_no_of_sides];
	calc_projections_on_standard_tetra(ctdrn, dc, proj);

	cout<<"\nVert-Vert Side Deviations: \n";
	for (int i=0; i<_no_of_sides; i++){
		double diff=_vertex_vertex_distances[i][2]-proj[i];
		deviate_dist[i]=diff;
		cout<<"Side "<<(i+1)<<": "<<deviate_dist[i]<<"; dist: "<<_vertex_vertex_distances[i][2]<<", Proj: "<<proj[i]<<endl;
	}
}

void geom::tetrahedron::deviation_angles_with_circum_centre(geom::tetrahedron ctdrn, double deviate_angle[]){
	cout<<"\nAngle Deviations CC:\n";
	for (int i=0; i<_no_of_angles_with_circum_centre; i++){
		double diff=get_angles_vertices_with_circum_centre()[i][2] - ctdrn.get_angles_vertices_with_circum_centre()[i][2];
		deviate_angle[i]=diff;
		cout<<"CC-Vert Angle "<<(i+1)<<": "<<deviate_angle[i]<<endl;
	}
}

void geom::tetrahedron::deviation_angles_vertex_vertex(geom::tetrahedron ctdrn, double deviate_angle[]){
	cout<<"\nAngle Deviations Vert-Vert:\n";
	for (int i=0; i<_no_of_vertex_vertex_angles; i++){
		double diff=_vertex_vertex_angles_tri[i][3] - ctdrn.get_vertex_vertex_angles()[i][3];
		deviate_angle[i]=diff;
		cout<<"Vert-Vert Angle "<<(i+1)<<": "<<deviate_angle[i]<<"  Stand:"<<get_vertex_vertex_angles()[i][3]<<"  Comp:"<<ctdrn.get_vertex_vertex_angles()[i][3]<<endl;
	}
}

geom::Point3D* geom::tetrahedron::get_vertices(){
	return _vertices;
}

double (*geom::tetrahedron::get_surface_areas())[4]{
	return _surface_areas;
}

int* geom::tetrahedron::get_vertex_order(){
	return _vertex_order;
}

geom::Point3D geom::tetrahedron::get_circum_centre(){
	return _circum_centre;
}

double geom::tetrahedron::get_circum_radius(){
	return _circum_radius;
}

geom::Point3D geom::tetrahedron::get_centre_of_mass(){
	return _centre_of_mass;
}

double (*geom::tetrahedron::get_circum_centre_to_vertex_distances())[2]{
	return _circum_centre_vertex_distances;
}

double (*geom::tetrahedron::get_vertex_vertex_distances())[3]{
	return _vertex_vertex_distances;
}

double (*geom::tetrahedron::get_angles_vertices_with_circum_centre())[3]{
	return _circum_centre_vertex_angles_tri;
}

double (*geom::tetrahedron::get_vertex_vertex_angles())[4]{
	return _vertex_vertex_angles_tri;
}


geom::trigonal_bipyramidal::trigonal_bipyramidal(){}

geom::trigonal_bipyramidal::trigonal_bipyramidal (vector<Point3D> plv, vector<Point3D> npv){
	_no_of_vertices=5;
	assert(_no_of_vertices==5 && "For Trigonal Bipyramidal No. of Vertices should be 5!!");

	_no_of_sides=9;
	assert(_no_of_sides==9 && "For Trigonal Bipyramidal No. of Sides should be 9!!");

	_no_of_planar_atoms=3;
	assert(_no_of_planar_atoms==3 && "For Trigonal Bipyramidal No. of Planar Atoms should be 3!!");

	_no_of_non_planar_atoms=2;
	assert(_no_of_non_planar_atoms==2 && "For Trigonal Bipyramidal No. of Non Planar Atoms should be 2!!");


	_no_of_vertex_vertex_angles=21;
	assert(_no_of_vertex_vertex_angles==21 && "For Trigonal Bipyramidal No. of Angles of Two Vertices with Another Vertex should be 21!!");

	_no_of_tetrahedron=2;

	for (int i=0; i<plv.size(); i++){
		_planar_vertices[i]=plv.at(i);
	}

	for (int i=0; i<npv.size(); i++){
		_non_planar_vertices[i]=npv.at(i);
	}

	measure_vertex_to_vertex_distances();
	measure_vertex_vertex_angles();
	find_vertex_order();
}

geom::trigonal_bipyramidal::trigonal_bipyramidal (Point3D* tribypyr_vert_pl, Point3D* tribypyr_vert_npl){
	_no_of_vertices=5;
	assert(_no_of_vertices==5 && "For Trigonal Bipyramidal No. of Vertices should be 5!!");

	_no_of_sides=9;
	assert(_no_of_sides==9 && "For Trigonal Bipyramidal No. of Sides should be 9!!");

	_no_of_planar_atoms=3;
	assert(_no_of_planar_atoms==3 && "For Trigonal Bipyramidal No. of Planar Atoms should be 3!!");

	_no_of_non_planar_atoms=2;
	assert(_no_of_non_planar_atoms==2 && "For Trigonal Bipyramidal No. of Non Planar Atoms should be 2!!");

	_no_of_vertex_vertex_angles=21;
	assert(_no_of_vertex_vertex_angles==21 && "For Trigonal Bipyramidal No. of Angles of Two Vertices with Another Vertex should be 21!!");

	_no_of_tetrahedron=2;

	for (int i=0; i<_no_of_planar_atoms; i++){
		_planar_vertices[i]=tribypyr_vert_pl[i];
		//cout<<"\nPl Vert: "<<_planar_vertices[i].X()<<", "<<_planar_vertices[i].Y()<<", "<<_planar_vertices[i].Z();
	}

	for (int i=0; i<_no_of_non_planar_atoms; i++){
		_non_planar_vertices[i]=tribypyr_vert_npl[i];
		//cout<<"\nNPl Vert: "<<_non_planar_vertices[i].X()<<", "<<_non_planar_vertices[i].Y()<<", "<<_non_planar_vertices[i].Z();
	}

	measure_vertex_to_vertex_distances();
	measure_vertex_vertex_angles();
	find_vertex_order();
}

void geom::trigonal_bipyramidal::measure_vertex_to_vertex_distances(){
	int k=0;
	//In Planar Vertices
	for (int i=0; i<_no_of_planar_atoms-1; i++){
		for (int j=i+1; j<_no_of_planar_atoms; j++){
			_vertex_vertex_distances_pl_pl[k][0]=i;
			_vertex_vertex_distances_pl_pl[k][1]=j;
			_vertex_vertex_distances_pl_pl[k][2]=_planar_vertices[i].dist(&_planar_vertices[j]);
			k++;
		}
	}

	//Between Planar and Non Planar Vertices
	int l=0;
	for (int i=0; i<_no_of_non_planar_atoms; i++){
		for (int j=0; j<_no_of_planar_atoms; j++){
			_vertex_vertex_distances_npl_pl[l][0]=i;
			_vertex_vertex_distances_npl_pl[l][1]=j;
			_vertex_vertex_distances_npl_pl[l][2]=_non_planar_vertices[i].dist(&_planar_vertices[j]);
			l++;
		}
	}

	cout<<"\n\nVertex-Vertex Distances Pl-Pl:\n";
	for (int i=0; i<3; i++){
		cout<<_vertex_vertex_distances_pl_pl[i][0]<<"  "<<_vertex_vertex_distances_pl_pl[i][1]<<"  "<<_vertex_vertex_distances_pl_pl[i][2]<<endl;
	}

	cout<<"\n\nVertex-Vertex Distances Npl-Pl:\n";
	for (int i=0; i<6; i++){
		cout<<_vertex_vertex_distances_npl_pl[i][0]<<"  "<<_vertex_vertex_distances_npl_pl[i][1]<<"  "<<_vertex_vertex_distances_npl_pl[i][2]<<endl;
	}
}

void geom::trigonal_bipyramidal::measure_vertex_vertex_angles(){
	int p=0;
	for (int i=0; i<_no_of_planar_atoms; i++){
		for (int j=0; j<_no_of_planar_atoms; j++){
			for (int k=0; k<_no_of_planar_atoms; k++){
				if (i!=j && j<k && i!=k){
					_vertex_vertex_angles_tri_pl_pl[p][0]=j;
					_vertex_vertex_angles_tri_pl_pl[p][1]=i;
					_vertex_vertex_angles_tri_pl_pl[p][2]=k;
					double ang=Point3D::angle_deg_triangle_law(&_planar_vertices[j], &_planar_vertices[i], &_planar_vertices[k]);
					_vertex_vertex_angles_tri_pl_pl[p][3]=ang;
					p++;
				}
			}
		}
	}

	p=0;
	Point3D test_vertices1[4]={_planar_vertices[0], _planar_vertices[1], _planar_vertices[2], _non_planar_vertices[0]};

	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			for (int k=0; k<4; k++){
				if (i!=j && j<k && i!=k){
					if ((i==1 && j==0 && k==2) || (i==2 && j==0 && k==1) || (i==0 && j==1 && k==2)){
						continue;
					}else{
						_vertex_vertex_angles_tri_npl_pl[p][0]=j;
						_vertex_vertex_angles_tri_npl_pl[p][1]=i;
						_vertex_vertex_angles_tri_npl_pl[p][2]=k;
						double ang=Point3D::angle_deg_triangle_law(&test_vertices1[j], &test_vertices1[i], &test_vertices1[k]);
						_vertex_vertex_angles_tri_npl_pl[p][3]=ang;
						p++;
					}
				}
			}
		}
	}

	Point3D test_vertices2[4]={_planar_vertices[0], _planar_vertices[1], _planar_vertices[2], _non_planar_vertices[1]};

	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			for (int k=0; k<4; k++){
				if (i!=j && j<k && i!=k){
					if ((i==1 && j==0 && k==2) || (i==2 && j==0 && k==1) || (i==0 && j==1 && k==2)){
						continue;
					}else{
						_vertex_vertex_angles_tri_npl_pl[p][0]=j;
						_vertex_vertex_angles_tri_npl_pl[p][1]=i;
						_vertex_vertex_angles_tri_npl_pl[p][2]=k;
						double ang=Point3D::angle_deg_triangle_law(&test_vertices2[j], &test_vertices2[i], &test_vertices2[k]);
						_vertex_vertex_angles_tri_npl_pl[p][3]=ang;
						p++;
					}
				}
			}
		}
	}

	cout<<"\n\nVertex-Vertex Angles Pl-Pl:\n";
	for (int i=0; i<_no_of_planar_atoms; i++){
		cout<<_vertex_vertex_angles_tri_pl_pl[i][0]<<"  "<<_vertex_vertex_angles_tri_pl_pl[i][1]<<"  "<<_vertex_vertex_angles_tri_pl_pl[i][2]<<"  "<<_vertex_vertex_angles_tri_pl_pl[i][3]<<endl;
	}

	cout<<"\n\nVertex-Vertex Angles Npl-Pl:\n";
	for (int i=0; i<p; i++){
		cout<<_vertex_vertex_angles_tri_npl_pl[i][0]<<"  "<<_vertex_vertex_angles_tri_npl_pl[i][1]<<"  "<<_vertex_vertex_angles_tri_npl_pl[i][2]<<"  "<<_vertex_vertex_angles_tri_npl_pl[i][3]<<endl;
	}
}

void geom::trigonal_bipyramidal::find_vertex_order(){
	double dist1=_planar_vertices[0].dist(&_planar_vertices[1]);
	double dist2=_planar_vertices[1].dist(&_planar_vertices[2]);
	double dist3=_planar_vertices[2].dist(&_planar_vertices[0]);

	cout<<"\n\nDIST Vertices: "<<dist1<<"  "<<dist2<<"  "<<dist3<<endl;

	if ((dist1>=dist2) && (dist1>=dist3)){
		if (_planar_vertices[0].X()<=_planar_vertices[1].X()){
			_pl_vertex_order[0]=0;
			_pl_vertex_order[1]=1;
		}else{
			_pl_vertex_order[0]=1;
			_pl_vertex_order[1]=0;
		}

		_pl_vertex_order[2]=2;
	}

	if ((dist2>=dist1) && (dist2>=dist3)){
		if (_planar_vertices[1].X()<=_planar_vertices[2].X()){
			_pl_vertex_order[0]=1;
			_pl_vertex_order[1]=2;
		}else{
			_pl_vertex_order[0]=2;
			_pl_vertex_order[1]=1;
		}

		_pl_vertex_order[2]=0;
	}

	if ((dist3>=dist1) && (dist3>=dist2)){
		if (_planar_vertices[2].X()<=_planar_vertices[0].X()){
			_pl_vertex_order[0]=2;
			_pl_vertex_order[1]=0;
		}else{
			_pl_vertex_order[0]=0;
			_pl_vertex_order[1]=2;
		}

		_pl_vertex_order[2]=1;
	}

	Point3D p12=_planar_vertices[1].minus(&_planar_vertices[0]);
	Point3D p13=_planar_vertices[2].minus(&_planar_vertices[0]);
	Point3D normal_to_plane=p12.cross(&p13);

	double product1=_non_planar_vertices[0].X()*normal_to_plane.X()+_non_planar_vertices[0].Y()*normal_to_plane.Y()+_non_planar_vertices[0].Z()*normal_to_plane.Z();
	double product2=_planar_vertices[0].X()*normal_to_plane.X()+_planar_vertices[0].Y()*normal_to_plane.Y()+_planar_vertices[0].Z()*normal_to_plane.Z();
	double product3=normal_to_plane.X()*normal_to_plane.X()+normal_to_plane.Y()*normal_to_plane.Y()+normal_to_plane.Z()*normal_to_plane.Z();

	double dist4=abs((product1-product2)/product3);

	double product4=_non_planar_vertices[1].X()*normal_to_plane.X()+_non_planar_vertices[1].Y()*normal_to_plane.Y()+_non_planar_vertices[1].Z()*normal_to_plane.Z();

	double dist5=abs((product4-product2)/product3);

	if (dist4>=dist5){
		_npl_vertex_order[0]=0;
		_npl_vertex_order[1]=1;
	}else{
		_npl_vertex_order[0]=1;
		_npl_vertex_order[1]=0;
	}
}

geom::Point3D* geom::trigonal_bipyramidal::get_pl_vertices(){
	return _planar_vertices;
}

geom::Point3D* geom::trigonal_bipyramidal::get_npl_vertices(){
	return _non_planar_vertices;
}

int* geom::trigonal_bipyramidal::get_pl_vertex_order(){
	return _pl_vertex_order;
}

int* geom::trigonal_bipyramidal::get_npl_vertex_order(){
	return _npl_vertex_order;
}

double (*geom::trigonal_bipyramidal::get_vertex_vertex_distances_pl_pl())[3]{
	return _vertex_vertex_distances_pl_pl;
}

double (*geom::trigonal_bipyramidal::get_vertex_vertex_distances_npl_pl())[3]{
	return _vertex_vertex_distances_npl_pl;
}

double (*geom::trigonal_bipyramidal::get_vertex_vertex_angles_pl_pl())[4]{
	return _vertex_vertex_angles_tri_pl_pl;
}

double (*geom::trigonal_bipyramidal::get_vertex_vertex_angles_npl_pl())[4]{
	return _vertex_vertex_angles_tri_npl_pl;
}

geom::tetragonal_pyramidal::tetragonal_pyramidal(){}

geom::tetragonal_pyramidal::tetragonal_pyramidal (vector<Point3D> plv, Point3D npv){
	_no_of_vertices=5;
	assert(_no_of_vertices==5 && "For Tetragonal Pyramidal No. of Vertices should be 5!!");

	_no_of_sides=8;
	assert(_no_of_sides==8 && "For Tetragonal Pyramidal No. of Sides should be 8!!");

	_no_of_planar_atoms=4;
	assert(_no_of_planar_atoms==4 && "For Tetragonal Pyramidal No. of Planar Atoms should be 4!!");

	_no_of_non_planar_atom=1;
	assert(_no_of_non_planar_atom==1 && "For Tetragonal Pyramidal No. of Non Planar Atom should be 1!!");

	_no_of_vertex_vertex_angles=16;
	assert(_no_of_vertex_vertex_angles==16 && "For Tetragonal Pyramidal No. of Angles of Two Vertices with Another Vertex should be 16!!");

	_no_of_tetrahedron=2;

	for (int i=0; i<plv.size(); i++){
		_planar_vertices[i]=plv.at(i);
	}

	_non_planar_vertex=npv;
	_npl_vertex_index=0;

	find_vertex_order();
	measure_vertex_to_vertex_distances();
	measure_vertex_vertex_angles();
}

geom::tetragonal_pyramidal::tetragonal_pyramidal (Point3D tetrapyr_vert_pl[], Point3D tetrapyr_vert_npl){
	_no_of_vertices=5;
	assert(_no_of_vertices==5 && "For Tetragonal Pyramidal No. of Vertices should be 5!!");

	_no_of_sides=8;
	assert(_no_of_sides==8 && "For Tetragonal Pyramidal No. of Sides should be 8!!");

	_no_of_planar_atoms=4;
	assert(_no_of_planar_atoms==4 && "For Tetragonal Pyramidal No. of Planar Atoms should be 4!!");

	_no_of_non_planar_atom=1;
	assert(_no_of_non_planar_atom==1 && "For Tetragonal Pyramidal No. of Non Planar Atom should be 1!!");

	_no_of_vertex_vertex_angles=16;
	assert(_no_of_vertex_vertex_angles==16 && "For Tetragonal Pyramidal No. of Angles of Two Vertices with Another Vertex should be 16!!");

	_no_of_tetrahedron=2;

	for (int i=0; i<_no_of_planar_atoms; i++){
		_planar_vertices[i]=tetrapyr_vert_pl[i];
		//cout<<"\nPl Vert: "<<_planar_vertices[i].X()<<", "<<_planar_vertices[i].Y()<<", "<<_planar_vertices[i].Z();
	}

	_non_planar_vertex=tetrapyr_vert_npl;
	_npl_vertex_index=0;

	find_vertex_order();
	measure_vertex_to_vertex_distances();
	measure_vertex_vertex_angles();
}

void geom::tetragonal_pyramidal::measure_vertex_to_vertex_distances(){
	int k=0;
	//In Planar Vertices
	_vertex_vertex_distances_pl_pl[0][0]=_pl_vertex_order[0];
	_vertex_vertex_distances_pl_pl[0][1]=_pl_vertex_order[1];
	_vertex_vertex_distances_pl_pl[0][2]=_planar_vertices[_pl_vertex_order[0]].dist(&_planar_vertices[_pl_vertex_order[1]]);

	_vertex_vertex_distances_pl_pl[1][0]=_pl_vertex_order[1];
	_vertex_vertex_distances_pl_pl[1][1]=_pl_vertex_order[2];
	_vertex_vertex_distances_pl_pl[1][2]=_planar_vertices[_pl_vertex_order[1]].dist(&_planar_vertices[_pl_vertex_order[2]]);

	_vertex_vertex_distances_pl_pl[2][0]=_pl_vertex_order[2];
	_vertex_vertex_distances_pl_pl[2][1]=_pl_vertex_order[3];
	_vertex_vertex_distances_pl_pl[2][2]=_planar_vertices[_pl_vertex_order[2]].dist(&_planar_vertices[_pl_vertex_order[3]]);

	_vertex_vertex_distances_pl_pl[3][0]=_pl_vertex_order[3];
	_vertex_vertex_distances_pl_pl[3][1]=_pl_vertex_order[0];
	_vertex_vertex_distances_pl_pl[3][2]=_planar_vertices[_pl_vertex_order[3]].dist(&_planar_vertices[_pl_vertex_order[0]]);

	//Between Planar and Non Planar Vertices
	int l=0;
	for (int j=0; j<_no_of_planar_atoms; j++){
		_vertex_vertex_distances_npl_pl[l][0]=_npl_vertex_index;
		_vertex_vertex_distances_npl_pl[l][1]=j;
		_vertex_vertex_distances_npl_pl[l][2]=_non_planar_vertex.dist(&_planar_vertices[j]);
		l++;
	}

	cout<<"\n\nVertex-Vertex Distances Pl-Pl:\n";
	for (int i=0; i<4; i++){
		cout<<_vertex_vertex_distances_pl_pl[i][0]<<"  "<<_vertex_vertex_distances_pl_pl[i][1]<<"  "<<_vertex_vertex_distances_pl_pl[i][2]<<endl;
	}

	cout<<"\n\nVertex-Vertex Distances Npl-Pl:\n";
	for (int i=0; i<4; i++){
		cout<<_vertex_vertex_distances_npl_pl[i][0]<<"  "<<_vertex_vertex_distances_npl_pl[i][1]<<"  "<<_vertex_vertex_distances_npl_pl[i][2]<<endl;
	}
}

void geom::tetragonal_pyramidal::measure_vertex_vertex_angles(){
	cout<<"\nStart Vert Vert Angle";
	_vertex_vertex_angles_tri_pl_pl[0][0]=_pl_vertex_order[1];
	_vertex_vertex_angles_tri_pl_pl[0][1]=_pl_vertex_order[0];
	_vertex_vertex_angles_tri_pl_pl[0][2]=_pl_vertex_order[3];
	_vertex_vertex_angles_tri_pl_pl[0][3]=Point3D::angle_deg_triangle_law(&_planar_vertices[_pl_vertex_order[1]], &_planar_vertices[_pl_vertex_order[0]], &_planar_vertices[_pl_vertex_order[3]]);

	_vertex_vertex_angles_tri_pl_pl[1][0]=_pl_vertex_order[0];
	_vertex_vertex_angles_tri_pl_pl[1][1]=_pl_vertex_order[1];
	_vertex_vertex_angles_tri_pl_pl[1][2]=_pl_vertex_order[2];
	_vertex_vertex_angles_tri_pl_pl[1][3]=Point3D::angle_deg_triangle_law(&_planar_vertices[_pl_vertex_order[0]], &_planar_vertices[_pl_vertex_order[1]], &_planar_vertices[_pl_vertex_order[2]]);

	_vertex_vertex_angles_tri_pl_pl[2][0]=_pl_vertex_order[1];
	_vertex_vertex_angles_tri_pl_pl[2][1]=_pl_vertex_order[2];
	_vertex_vertex_angles_tri_pl_pl[2][2]=_pl_vertex_order[3];
	_vertex_vertex_angles_tri_pl_pl[2][3]=Point3D::angle_deg_triangle_law(&_planar_vertices[_pl_vertex_order[1]], &_planar_vertices[_pl_vertex_order[2]], &_planar_vertices[_pl_vertex_order[3]]);

	_vertex_vertex_angles_tri_pl_pl[3][0]=_pl_vertex_order[0];
	_vertex_vertex_angles_tri_pl_pl[3][1]=_pl_vertex_order[3];
	_vertex_vertex_angles_tri_pl_pl[3][2]=_pl_vertex_order[2];
	_vertex_vertex_angles_tri_pl_pl[3][3]=Point3D::angle_deg_triangle_law(&_planar_vertices[_pl_vertex_order[0]], &_planar_vertices[_pl_vertex_order[3]], &_planar_vertices[_pl_vertex_order[2]]);

	int p=0, i;
	for (i=0; i<(_no_of_planar_atoms-1); i++){
		_vertex_vertex_angles_tri_npl_pl[p][0]=i;
		_vertex_vertex_angles_tri_npl_pl[p][1]=_npl_vertex_index;
		_vertex_vertex_angles_tri_npl_pl[p][2]=(i+1);
		double ang=Point3D::angle_deg_triangle_law(&_planar_vertices[i], &_non_planar_vertex, &_planar_vertices[i+1]);
		_vertex_vertex_angles_tri_npl_pl[p][3]=ang;

		p++;
		_vertex_vertex_angles_tri_npl_pl[p][0]=_npl_vertex_index;
		_vertex_vertex_angles_tri_npl_pl[p][1]=i;
		_vertex_vertex_angles_tri_npl_pl[p][2]=(i+1);
		ang=Point3D::angle_deg_triangle_law(&_non_planar_vertex, &_planar_vertices[i], &_planar_vertices[i+1]);
		_vertex_vertex_angles_tri_npl_pl[p][3]=ang;

		p++;
		_vertex_vertex_angles_tri_npl_pl[p][0]=_npl_vertex_index;
		_vertex_vertex_angles_tri_npl_pl[p][1]=(i+1);
		_vertex_vertex_angles_tri_npl_pl[p][2]=i;
		ang=Point3D::angle_deg_triangle_law(&_non_planar_vertex, &_planar_vertices[i+1], &_planar_vertices[i]);
		_vertex_vertex_angles_tri_npl_pl[p][3]=ang;
		p++;
	}

	_vertex_vertex_angles_tri_npl_pl[p][0]=i;
	_vertex_vertex_angles_tri_npl_pl[p][1]=_npl_vertex_index;
	_vertex_vertex_angles_tri_npl_pl[p][2]=0;
	double ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[i], &_non_planar_vertex, &_planar_vertices[0]);
	_vertex_vertex_angles_tri_npl_pl[p][3]=ang1;

	p++;
	_vertex_vertex_angles_tri_npl_pl[p][0]=_npl_vertex_index;
	_vertex_vertex_angles_tri_npl_pl[p][1]=i;
	_vertex_vertex_angles_tri_npl_pl[p][2]=0;
	ang1=Point3D::angle_deg_triangle_law(&_non_planar_vertex, &_planar_vertices[i], &_planar_vertices[0]);
	_vertex_vertex_angles_tri_npl_pl[p][3]=ang1;

	p++;
	_vertex_vertex_angles_tri_npl_pl[p][0]=_npl_vertex_index;
	_vertex_vertex_angles_tri_npl_pl[p][1]=0;
	_vertex_vertex_angles_tri_npl_pl[p][2]=i;
	ang1=Point3D::angle_deg_triangle_law(&_non_planar_vertex, &_planar_vertices[0], &_planar_vertices[i]);
	_vertex_vertex_angles_tri_npl_pl[p][3]=ang1;

	cout<<"\n\nVertex-Vertex Angles Pl-Pl:\n";
	for (int x=0; x<4; x++){
		cout<<_vertex_vertex_angles_tri_pl_pl[x][0]<<"  "<<_vertex_vertex_angles_tri_pl_pl[x][1]<<"  "<<_vertex_vertex_angles_tri_pl_pl[x][2]<<"  "<<_vertex_vertex_angles_tri_pl_pl[x][3]<<endl;
	}

	cout<<"\n\nVertex-Vertex Angles Npl-Pl:\n";
	for (int x=0; x<=p; x++){
		cout<<_vertex_vertex_angles_tri_npl_pl[x][0]<<"  "<<_vertex_vertex_angles_tri_npl_pl[x][1]<<"  "<<_vertex_vertex_angles_tri_npl_pl[x][2]<<"  "<<_vertex_vertex_angles_tri_npl_pl[x][3]<<endl;
	}
	cout<<"\nEnd Vert Vert Angle";
}

void geom::tetragonal_pyramidal::find_vertex_order(){
	cout<<"\nOK Order Started";
	double diff1[6], diff2[6];
	int v_order[6][4];
	double ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[3], &_planar_vertices[0], &_planar_vertices[1]);
	double ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[1], &_planar_vertices[2]);
	double ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[2], &_planar_vertices[3]);
	double ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[3], &_planar_vertices[0]);
	diff1[0]=180.0-(ang1+ang3);
	diff2[0]=180.0-(ang2+ang4);
	cout<<"\n 1st: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[0]<<", "<<diff2[0];

	v_order[0][0]=0;
	v_order[0][1]=1;
	v_order[0][2]=2;
	v_order[0][3]=3;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[0], &_planar_vertices[1]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[1], &_planar_vertices[3]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[3], &_planar_vertices[2]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[2], &_planar_vertices[3]);
	diff1[1]=180.0-(ang1+ang3);
	diff2[1]=180.0-(ang2+ang4);
	cout<<"\n 2nd: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[1]<<", "<<diff2[1];

	v_order[1][0]=0;
	v_order[1][1]=1;
	v_order[1][2]=3;
	v_order[1][3]=2;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[0], &_planar_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[2], &_planar_vertices[1]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[1], &_planar_vertices[3]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[3], &_planar_vertices[0]);
	diff1[2]=180.0-(ang1+ang3);
	diff2[2]=180.0-(ang2+ang4);
	cout<<"\n 3rd: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[2]<<", "<<diff2[2];

	v_order[2][0]=0;
	v_order[2][1]=2;
	v_order[2][2]=1;
	v_order[2][3]=3;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[0], &_planar_vertices[2]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[2], &_planar_vertices[3]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[3], &_planar_vertices[1]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[1], &_planar_vertices[3]);
	diff1[3]=180.0-(ang1+ang3);
	diff2[3]=180.0-(ang2+ang4);
	cout<<"\n 4th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[3]<<", "<<diff2[3];

	v_order[3][0]=0;
	v_order[3][1]=2;
	v_order[3][2]=3;
	v_order[3][3]=1;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[0], &_planar_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[3], &_planar_vertices[1]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[3], &_planar_vertices[1], &_planar_vertices[2]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[2], &_planar_vertices[0]);
	diff1[4]=180.0-(ang1+ang3);
	diff2[4]=180.0-(ang2+ang4);
	cout<<"\n 5th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[4]<<", "<<diff2[4];

	v_order[4][0]=0;
	v_order[4][1]=3;
	v_order[4][2]=1;
	v_order[4][3]=2;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[0], &_planar_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[3], &_planar_vertices[2]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[3], &_planar_vertices[2], &_planar_vertices[1]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[1], &_planar_vertices[2]);
	diff1[5]=180.0-(ang1+ang3);
	diff2[5]=180.0-(ang2+ang4);
	cout<<"\n 6th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[5]<<", "<<diff2[5];

	v_order[5][0]=0;
	v_order[5][1]=3;
	v_order[5][2]=2;
	v_order[5][3]=1;

	double diff1_min=abs(diff1[0]), diff2_min=abs(diff2[0]);
	int min_index1=0, min_index2=0;
	for (int i=1; i<6; i++){
		if (abs(diff1[i])<abs(diff1_min)){
			diff1_min=abs(diff1[i]);
			min_index1=i;
		}

		if (abs(diff2[i])<abs(diff2_min)){
			diff2_min=abs(diff2[i]);
			min_index2=i;
		}
	}

	if (abs(diff1_min)<=abs(diff2_min)){
		_pl_vertex_order[0]=v_order[min_index1][0];
		_pl_vertex_order[1]=v_order[min_index1][1];
		_pl_vertex_order[2]=v_order[min_index1][2];
		_pl_vertex_order[3]=v_order[min_index1][3];
		cout<<"\nSelected Order Index: "<<min_index1<<"; Diff Min: "<<diff1_min;
	}else{
		_pl_vertex_order[0]=v_order[min_index2][0];
		_pl_vertex_order[1]=v_order[min_index2][1];
		_pl_vertex_order[2]=v_order[min_index2][2];
		_pl_vertex_order[3]=v_order[min_index2][3];
		cout<<"\nSelected Order Index: "<<min_index2<<"; Diff Min: "<<diff2_min;
	}

	cout<<"\nOK Order End: "<<_pl_vertex_order[0]<<", "<<_pl_vertex_order[1]<<", "<<_pl_vertex_order[2]<<", "<<_pl_vertex_order[3];
}

geom::Point3D* geom::tetragonal_pyramidal::get_pl_vertices(){
	return _planar_vertices;
}

geom::Point3D geom::tetragonal_pyramidal::get_npl_vertex(){
	return _non_planar_vertex;
}

int* geom::tetragonal_pyramidal::get_pl_vertex_order(){
	return _pl_vertex_order;
}

double (*geom::tetragonal_pyramidal::get_vertex_vertex_distances_pl_pl())[3]{
	return _vertex_vertex_distances_pl_pl;
}

double (*geom::tetragonal_pyramidal::get_vertex_vertex_distances_npl_pl())[3]{
	return _vertex_vertex_distances_npl_pl;
}

double (*geom::tetragonal_pyramidal::get_vertex_vertex_angles_pl_pl())[4]{
	return _vertex_vertex_angles_tri_pl_pl;
}

double (*geom::tetragonal_pyramidal::get_vertex_vertex_angles_npl_pl())[4]{
	return _vertex_vertex_angles_tri_npl_pl;
}

geom::octahedral::octahedral(){}

geom::octahedral::octahedral (vector<Point3D> plv, vector<Point3D> npv){
	_no_of_vertices=6;
	assert(_no_of_vertices==6 && "For Octahedral No. of Vertices should be 6!!");

	_no_of_sides=12;
	assert(_no_of_sides==12 && "For Octahedral No. of Sides should be 12!!");

	_no_of_planar_atoms=4;
	assert(_no_of_planar_atoms==4 && "For Octahedral No. of Planar Atoms should be 4!!");

	_no_of_non_planar_atoms=2;
	assert(_no_of_non_planar_atoms==2 && "For Octahedral No. of Non Planar Atom should be 2!!");

	_no_of_vertex_vertex_angles=28;
	assert(_no_of_vertex_vertex_angles==28 && "For Octahedral No. of Angles of Two Vertices with Another Vertex should be 28!!");

	_no_of_tetrahedron=4;

	for (int i=0; i<plv.size(); i++){
		_planar_vertices[i]=plv.at(i);
	}

	for (int i=0; i<plv.size(); i++){
		_non_planar_vertices[i]=npv.at(i);
	}

	find_vertex_order();
	measure_vertex_to_vertex_distances();
	measure_vertex_vertex_angles();
}

geom::octahedral::octahedral (Point3D octa_vert_pl[], Point3D octa_vert_npl[]){
	_no_of_vertices=6;
	assert(_no_of_vertices==6 && "For Octahedral No. of Vertices should be 6!!");

	_no_of_sides=12;
	assert(_no_of_sides==12 && "For Octahedral No. of Sides should be 12!!");

	_no_of_planar_atoms=4;
	assert(_no_of_planar_atoms==4 && "For Octahedral No. of Planar Atoms should be 4!!");

	_no_of_non_planar_atoms=2;
	assert(_no_of_non_planar_atoms==2 && "For Octahedral No. of Non Planar Atom should be 2!!");

	_no_of_vertex_vertex_angles=28;
	assert(_no_of_vertex_vertex_angles==28 && "For Octahedral No. of Angles of Two Vertices with Another Vertex should be 28!!");

	_no_of_tetrahedron=4;

	for (int i=0; i<_no_of_planar_atoms; i++){
		_planar_vertices[i]=octa_vert_pl[i];
		//cout<<"\nPl Vert: "<<_planar_vertices[i].X()<<", "<<_planar_vertices[i].Y()<<", "<<_planar_vertices[i].Z();
	}

	for (int i=0; i<_no_of_planar_atoms; i++){
		_non_planar_vertices[i]=octa_vert_npl[i];
		//cout<<"\nNpl Vert: "<<_non_planar_vertices[i].X()<<", "<<_non_planar_vertices[i].Y()<<", "<<_non_planar_vertices[i].Z();
	}

	find_vertex_order();
	measure_vertex_to_vertex_distances();
	measure_vertex_vertex_angles();
}

void geom::octahedral::measure_vertex_to_vertex_distances(){
	int k=0;
	//In Planar Vertices
	_vertex_vertex_distances_pl_pl[0][0]=_pl_vertex_order[0];
	_vertex_vertex_distances_pl_pl[0][1]=_pl_vertex_order[1];
	_vertex_vertex_distances_pl_pl[0][2]=_planar_vertices[_pl_vertex_order[0]].dist(&_planar_vertices[_pl_vertex_order[1]]);

	_vertex_vertex_distances_pl_pl[1][0]=_pl_vertex_order[1];
	_vertex_vertex_distances_pl_pl[1][1]=_pl_vertex_order[2];
	_vertex_vertex_distances_pl_pl[1][2]=_planar_vertices[_pl_vertex_order[1]].dist(&_planar_vertices[_pl_vertex_order[2]]);

	_vertex_vertex_distances_pl_pl[2][0]=_pl_vertex_order[2];
	_vertex_vertex_distances_pl_pl[2][1]=_pl_vertex_order[3];
	_vertex_vertex_distances_pl_pl[2][2]=_planar_vertices[_pl_vertex_order[2]].dist(&_planar_vertices[_pl_vertex_order[3]]);

	_vertex_vertex_distances_pl_pl[3][0]=_pl_vertex_order[3];
	_vertex_vertex_distances_pl_pl[3][1]=_pl_vertex_order[0];
	_vertex_vertex_distances_pl_pl[3][2]=_planar_vertices[_pl_vertex_order[3]].dist(&_planar_vertices[_pl_vertex_order[0]]);

	//Between Planar and Non Planar Vertices
	int l=0;
	for (int i=0; i<_no_of_non_planar_atoms; i++){
		for (int j=0; j<_no_of_planar_atoms; j++){
			_vertex_vertex_distances_npl_pl[l][0]=i;
			_vertex_vertex_distances_npl_pl[l][1]=j;
			_vertex_vertex_distances_npl_pl[l][2]=_non_planar_vertices[i].dist(&_planar_vertices[j]);
			l++;
		}
	}

	cout<<"\n\nVertex-Vertex Distances Pl-Pl:\n";
	for (int i=0; i<4; i++){
		cout<<_vertex_vertex_distances_pl_pl[i][0]<<"  "<<_vertex_vertex_distances_pl_pl[i][1]<<"  "<<_vertex_vertex_distances_pl_pl[i][2]<<endl;
	}

	cout<<"\n\nVertex-Vertex Distances Npl-Pl:\n";
	for (int i=0; i<8; i++){
		cout<<_vertex_vertex_distances_npl_pl[i][0]<<"  "<<_vertex_vertex_distances_npl_pl[i][1]<<"  "<<_vertex_vertex_distances_npl_pl[i][2]<<endl;
	}
}

void geom::octahedral::measure_vertex_vertex_angles(){
	_vertex_vertex_angles_tri_pl_pl[0][0]=_pl_vertex_order[1];
	_vertex_vertex_angles_tri_pl_pl[0][1]=_pl_vertex_order[0];
	_vertex_vertex_angles_tri_pl_pl[0][2]=_pl_vertex_order[3];
	_vertex_vertex_angles_tri_pl_pl[0][3]=Point3D::angle_deg_triangle_law(&_planar_vertices[_pl_vertex_order[1]], &_planar_vertices[_pl_vertex_order[0]], &_planar_vertices[_pl_vertex_order[3]]);

	_vertex_vertex_angles_tri_pl_pl[1][0]=_pl_vertex_order[0];
	_vertex_vertex_angles_tri_pl_pl[1][1]=_pl_vertex_order[1];
	_vertex_vertex_angles_tri_pl_pl[1][2]=_pl_vertex_order[2];
	_vertex_vertex_angles_tri_pl_pl[1][3]=Point3D::angle_deg_triangle_law(&_planar_vertices[_pl_vertex_order[0]], &_planar_vertices[_pl_vertex_order[1]], &_planar_vertices[_pl_vertex_order[2]]);

	_vertex_vertex_angles_tri_pl_pl[2][0]=_pl_vertex_order[1];
	_vertex_vertex_angles_tri_pl_pl[2][1]=_pl_vertex_order[2];
	_vertex_vertex_angles_tri_pl_pl[2][2]=_pl_vertex_order[3];
	_vertex_vertex_angles_tri_pl_pl[2][3]=Point3D::angle_deg_triangle_law(&_planar_vertices[_pl_vertex_order[1]], &_planar_vertices[_pl_vertex_order[2]], &_planar_vertices[_pl_vertex_order[3]]);

	_vertex_vertex_angles_tri_pl_pl[3][0]=_pl_vertex_order[0];
	_vertex_vertex_angles_tri_pl_pl[3][1]=_pl_vertex_order[3];
	_vertex_vertex_angles_tri_pl_pl[3][2]=_pl_vertex_order[2];
	_vertex_vertex_angles_tri_pl_pl[3][3]=Point3D::angle_deg_triangle_law(&_planar_vertices[_pl_vertex_order[0]], &_planar_vertices[_pl_vertex_order[3]], &_planar_vertices[_pl_vertex_order[2]]);


	int p=0;
	Point3D test_vertices1[4]={_planar_vertices[_pl_vertex_order[0]], _planar_vertices[_pl_vertex_order[1]], _planar_vertices[_pl_vertex_order[3]], _non_planar_vertices[_npl_vertex_order[0]]};

	int k=3;
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			for (int k=0; k<4; k++){
				if (i!=j && j<k && i!=k){
					if ((i==1 && j==0 && k==2) || (i==2 && j==0 && k==1) || (i==0 && j==1 && k==2) || (i==1 && j==2 && k==3) || (i==2 && j==1 && k==3) || (i==3 && j==1 && k==2)){
						continue;
					}else{
						_vertex_vertex_angles_tri_npl_pl[p][0]=j;
						_vertex_vertex_angles_tri_npl_pl[p][1]=i;
						_vertex_vertex_angles_tri_npl_pl[p][2]=k;
						double ang=Point3D::angle_deg_triangle_law(&test_vertices1[j], &test_vertices1[i], &test_vertices1[k]);
						_vertex_vertex_angles_tri_npl_pl[p][3]=ang;
						p++;
					}
				}
			}
		}
	}

	Point3D test_vertices2[4]={_planar_vertices[_pl_vertex_order[3]], _planar_vertices[_pl_vertex_order[1]], _planar_vertices[_pl_vertex_order[2]], _non_planar_vertices[_npl_vertex_order[0]]};

	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			for (int k=0; k<4; k++){
				if (i!=j && j<k && i!=k){
					if ((i==1 && j==0 && k==2) || (i==2 && j==0 && k==1) || (i==0 && j==1 && k==2) || (i==0 && j==1 && k==3) || (i==1 && j==0 && k==3) || (i==3 && j==0 && k==1)){
						continue;
					}else{
						_vertex_vertex_angles_tri_npl_pl[p][0]=j;
						_vertex_vertex_angles_tri_npl_pl[p][1]=i;
						_vertex_vertex_angles_tri_npl_pl[p][2]=k;
						double ang=Point3D::angle_deg_triangle_law(&test_vertices2[j], &test_vertices2[i], &test_vertices2[k]);
						_vertex_vertex_angles_tri_npl_pl[p][3]=ang;
						p++;
					}
				}
			}
		}
	}

	Point3D test_vertices3[4]={_planar_vertices[_pl_vertex_order[0]], _planar_vertices[_pl_vertex_order[1]], _planar_vertices[_pl_vertex_order[3]], _non_planar_vertices[_npl_vertex_order[1]]};

	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			for (int k=0; k<4; k++){
				if (i!=j && j<k && i!=k){
					if ((i==1 && j==0 && k==2) || (i==2 && j==0 && k==1) || (i==0 && j==1 && k==2) || (i==1 && j==2 && k==3) || (i==2 && j==1 && k==3) || (i==3 && j==1 && k==2)){
						continue;
					}else{
						_vertex_vertex_angles_tri_npl_pl[p][0]=j;
						_vertex_vertex_angles_tri_npl_pl[p][1]=i;
						_vertex_vertex_angles_tri_npl_pl[p][2]=k;
						double ang=Point3D::angle_deg_triangle_law(&test_vertices3[j], &test_vertices3[i], &test_vertices3[k]);
						_vertex_vertex_angles_tri_npl_pl[p][3]=ang;
						p++;
					}
				}
			}
		}
	}

	Point3D test_vertices4[4]={_planar_vertices[_pl_vertex_order[3]], _planar_vertices[_pl_vertex_order[1]], _planar_vertices[_pl_vertex_order[2]], _non_planar_vertices[_npl_vertex_order[1]]};

	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			for (int k=0; k<4; k++){
				if (i!=j && j<k && i!=k){
					if ((i==1 && j==0 && k==2) || (i==2 && j==0 && k==1) || (i==0 && j==1 && k==2) || (i==0 && j==1 && k==3) || (i==1 && j==0 && k==3) || (i==3 && j==0 && k==1)){
						continue;
					}else{
						_vertex_vertex_angles_tri_npl_pl[p][0]=j;
						_vertex_vertex_angles_tri_npl_pl[p][1]=i;
						_vertex_vertex_angles_tri_npl_pl[p][2]=k;
						double ang=Point3D::angle_deg_triangle_law(&test_vertices4[j], &test_vertices4[i], &test_vertices4[k]);
						_vertex_vertex_angles_tri_npl_pl[p][3]=ang;
						p++;
					}
				}
			}
		}
	}

	cout<<"\n\nVertex-Vertex Angles Pl-Pl:\n";
	for (int x=0; x<4; x++){
		cout<<_vertex_vertex_angles_tri_pl_pl[x][0]<<"  "<<_vertex_vertex_angles_tri_pl_pl[x][1]<<"  "<<_vertex_vertex_angles_tri_pl_pl[x][2]<<"  "<<_vertex_vertex_angles_tri_pl_pl[x][3]<<endl;
	}

	cout<<"\n\nVertex-Vertex Angles Npl-Pl:\n";
	for (int x=0; x<p; x++){
		cout<<_vertex_vertex_angles_tri_npl_pl[x][0]<<"  "<<_vertex_vertex_angles_tri_npl_pl[x][1]<<"  "<<_vertex_vertex_angles_tri_npl_pl[x][2]<<"  "<<_vertex_vertex_angles_tri_npl_pl[x][3]<<endl;
	}
}

void geom::octahedral::find_vertex_order(){
	double diff1[6], diff2[6];
	int v_order[6][4];
	double ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[3], &_planar_vertices[0], &_planar_vertices[1]);
	double ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[1], &_planar_vertices[2]);
	double ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[2], &_planar_vertices[3]);
	double ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[3], &_planar_vertices[0]);
	diff1[0]=180.0-(ang1+ang3);
	diff2[0]=180.0-(ang2+ang4);
	cout<<"\n 1st: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[0]<<", "<<diff2[0];

	v_order[0][0]=0;
	v_order[0][1]=1;
	v_order[0][2]=2;
	v_order[0][3]=3;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[0], &_planar_vertices[1]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[1], &_planar_vertices[3]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[3], &_planar_vertices[2]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[2], &_planar_vertices[3]);
	diff1[1]=180.0-(ang1+ang3);
	diff2[1]=180.0-(ang2+ang4);
	cout<<"\n 2nd: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[1]<<", "<<diff2[1];

	v_order[1][0]=0;
	v_order[1][1]=1;
	v_order[1][2]=3;
	v_order[1][3]=2;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[0], &_planar_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[2], &_planar_vertices[1]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[1], &_planar_vertices[3]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[3], &_planar_vertices[0]);
	diff1[2]=180.0-(ang1+ang3);
	diff2[2]=180.0-(ang2+ang4);
	cout<<"\n 3rd: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[2]<<", "<<diff2[2];

	v_order[2][0]=0;
	v_order[2][1]=2;
	v_order[2][2]=1;
	v_order[2][3]=3;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[0], &_planar_vertices[2]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[2], &_planar_vertices[3]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[3], &_planar_vertices[1]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[1], &_planar_vertices[3]);
	diff1[3]=180.0-(ang1+ang3);
	diff2[3]=180.0-(ang2+ang4);
	cout<<"\n 4th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[3]<<", "<<diff2[3];

	v_order[3][0]=0;
	v_order[3][1]=2;
	v_order[3][2]=3;
	v_order[3][3]=1;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[2], &_planar_vertices[0], &_planar_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[3], &_planar_vertices[1]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[3], &_planar_vertices[1], &_planar_vertices[2]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[2], &_planar_vertices[0]);
	diff1[4]=180.0-(ang1+ang3);
	diff2[4]=180.0-(ang2+ang4);
	cout<<"\n 5th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[4]<<", "<<diff2[4];

	v_order[4][0]=0;
	v_order[4][1]=3;
	v_order[4][2]=1;
	v_order[4][3]=2;

	ang1=Point3D::angle_deg_triangle_law(&_planar_vertices[1], &_planar_vertices[0], &_planar_vertices[3]);
	ang2=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[3], &_planar_vertices[2]);
	ang3=Point3D::angle_deg_triangle_law(&_planar_vertices[3], &_planar_vertices[2], &_planar_vertices[1]);
	ang4=Point3D::angle_deg_triangle_law(&_planar_vertices[0], &_planar_vertices[1], &_planar_vertices[2]);
	diff1[5]=180.0-(ang1+ang3);
	diff2[5]=180.0-(ang2+ang4);
	cout<<"\n 6th: "<<(ang1+ang3)<<", "<<(ang2+ang4)<<"; Diff: "<<diff1[5]<<", "<<diff2[5];

	v_order[5][0]=0;
	v_order[5][1]=3;
	v_order[5][2]=2;
	v_order[5][3]=1;

	double diff1_min=abs(diff1[0]), diff2_min=abs(diff2[0]);
	int min_index1=0, min_index2=0;
	for (int i=1; i<6; i++){
		if (abs(diff1[i])<abs(diff1_min)){
			diff1_min=abs(diff1[i]);
			min_index1=i;
		}

		if (abs(diff2[i])<abs(diff2_min)){
			diff2_min=abs(diff2[i]);
			min_index2=i;
		}
	}

	if (abs(diff1_min)<=abs(diff2_min)){
		_pl_vertex_order[0]=v_order[min_index1][0];
		_pl_vertex_order[1]=v_order[min_index1][1];
		_pl_vertex_order[2]=v_order[min_index1][2];
		_pl_vertex_order[3]=v_order[min_index1][3];
		cout<<"\nSelected Order Index: "<<min_index1<<"; Diff Min: "<<diff1_min;
	}else{
		_pl_vertex_order[0]=v_order[min_index2][0];
		_pl_vertex_order[1]=v_order[min_index2][1];
		_pl_vertex_order[2]=v_order[min_index2][2];
		_pl_vertex_order[3]=v_order[min_index2][3];
		cout<<"\nSelected Order Index: "<<min_index2<<"; Diff Min: "<<diff2_min;
	}

	Point3D p12=_planar_vertices[1].minus(&_planar_vertices[0]);
	Point3D p13=_planar_vertices[2].minus(&_planar_vertices[0]);
	Point3D normal_to_plane=p12.cross(&p13);

	double product1=_non_planar_vertices[0].X()*normal_to_plane.X()+_non_planar_vertices[0].Y()*normal_to_plane.Y()+_non_planar_vertices[0].Z()*normal_to_plane.Z();
	double product2=_planar_vertices[0].X()*normal_to_plane.X()+_planar_vertices[0].Y()*normal_to_plane.Y()+_planar_vertices[0].Z()*normal_to_plane.Z();
	double product3=normal_to_plane.X()*normal_to_plane.X()+normal_to_plane.Y()*normal_to_plane.Y()+normal_to_plane.Z()*normal_to_plane.Z();

	double dist4=abs((product1-product2)/product3);

	double product4=_non_planar_vertices[1].X()*normal_to_plane.X()+_non_planar_vertices[1].Y()*normal_to_plane.Y()+_non_planar_vertices[1].Z()*normal_to_plane.Z();

	double dist5=abs((product4-product2)/product3);

	if (dist4>=dist5){
		_npl_vertex_order[0]=0;
		_npl_vertex_order[1]=1;
	}else{
		_npl_vertex_order[0]=1;
		_npl_vertex_order[1]=0;
	}

	cout<<"\nInitial Planar Vertex Order: "<<_pl_vertex_order[0]<<", "<<_pl_vertex_order[1]<<", "<<_pl_vertex_order[2]<<", "<<_pl_vertex_order[3];
	cout<<"\nInitial Non Planar Vertex Order: "<<_npl_vertex_order[0]<<", "<<_npl_vertex_order[1];
}

geom::Point3D* geom::octahedral::get_pl_vertices(){
	return _planar_vertices;
}

geom::Point3D* geom::octahedral::get_npl_vertices(){
	return _non_planar_vertices;
}

int* geom::octahedral::get_pl_vertex_order(){
	return _pl_vertex_order;
}

int* geom::octahedral::get_npl_vertex_order(){
	return _npl_vertex_order;
}

double (*geom::octahedral::get_vertex_vertex_distances_pl_pl())[3]{
	return _vertex_vertex_distances_pl_pl;
}

double (*geom::octahedral::get_vertex_vertex_distances_npl_pl())[3]{
	return _vertex_vertex_distances_npl_pl;
}

double (*geom::octahedral::get_vertex_vertex_angles_pl_pl())[4]{
	return _vertex_vertex_angles_tri_pl_pl;
}

double (*geom::octahedral::get_vertex_vertex_angles_npl_pl())[4]{
	return _vertex_vertex_angles_tri_npl_pl;
}
