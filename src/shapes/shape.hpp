#ifndef SHAPES_SHAPE_HPP
#define SHAPES_SHAPE_HPP

#include <Eigen/Dense>
#include "display/drawable.hpp"
#include "collision/sdf.hpp"

#define DIM 3
#define DEFAULT_COLOR 0.5, 0.5, 0.0
#define DEFAULT_POSITION 0, 0, 0
#define DEFAULT_VELOCITY 0, 0, 0
#define DEFAULT_MASS 1

class Shape : public Drawable, public SDF
{
public:
	Eigen::Vector3d color;
	Eigen::Vector3d x;
	Eigen::Vector3d v;
	unsigned int i;	// Index in the list
	double m;
public:
	Shape() : color(DEFAULT_COLOR), x(DEFAULT_POSITION), v(DEFAULT_VELOCITY), m(DEFAULT_MASS) {}
	~Shape(){}
	void setColor(const Eigen::Vector3d& _color);
	void getColor(Eigen::Vector3d& _color);
	void setPosition(const Eigen::Vector3d& _x);
	void getPosition(Eigen::Vector3d& _x);
	void setVelocity(const Eigen::Vector3d& _v);
	void getVelocity(Eigen::Vector3d& _v);
	void setIndex(const unsigned int _i);
	unsigned int getIndex();
	void setMass(const double _m);
	double getMass();
	
	virtual void draw() = 0;

	virtual double signedDistance(Eigen::Vector3d x) = 0;
	virtual Eigen::Vector3d normal(Eigen::Vector3d x) = 0;
};

#define DEFAULT_RADIUS 0.8

class Sphere : public Shape
{
private:
	double r;
public:
	Sphere(double _r = DEFAULT_RADIUS) : r(_r){}
	// Sphere(Eigen::Vector3d _c = DEFAULT_CENTER, double _r = DEFAULT_RADIUS, Eigen::Vector3d _color) : c(_c), r(_r){}
	~Sphere() {}

	void draw();

	double signedDistance(Eigen::Vector3d x);
	Eigen::Vector3d normal(Eigen::Vector3d x);
	
};

#define DEFAULT_V0 Eigen::Vector3d(-1, 0, -1)
#define DEFAULT_V1 Eigen::Vector3d(1, 0, -1)
#define DEFAULT_V2 Eigen::Vector3d(1, 0, 1)

class Plane : public Shape
{
private:
	
public:
	Eigen::Vector3d v0, v1, v2;	// Only 3 vertices are needed to determine a plane
	Plane(Eigen::Vector3d _v0 = DEFAULT_V0, Eigen::Vector3d _v1 = DEFAULT_V1, Eigen::Vector3d _v2 = DEFAULT_V2) : v0(_v0), v1(_v1), v2(_v2) {}
	~Plane() {}
	void draw();
	
	double signedDistance(Eigen::Vector3d x);
	Eigen::Vector3d normal(Eigen::Vector3d x);
};

#endif