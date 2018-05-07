#ifndef FORCES_FORCE_HPP
#define FORCES_FORCE_HPP

#include "shapes/shape.hpp"

typedef Eigen::Triplet<double> T;
class Force
{
public:
	Force() {}
	~Force() {}
	virtual void applyForces(Eigen::VectorXd& f) = 0;
	virtual void applyJacobians(std::vector<T> &triplesJx, std::vector<T> &triplesJv) = 0;
};

#define GRAVITY_ACCELERATION -9.8
class Gravity : public Force
{
public:
	Shape* shape;
	Eigen::Vector3d g;
public:
	Gravity(Shape* _shape) : shape(_shape), g(0, GRAVITY_ACCELERATION, 0) {}
	~Gravity() {}
	void applyForces(Eigen::VectorXd& f);
	void applyJacobians(std::vector<T> &triplesJx, std::vector<T> &triplesJv);
};

class SpringForce : public Force
{
public:
	double ks, kd, l;
	Shape *p0, *p1;
	SpringForce(double _ks, double _kd, double _l, Shape* _p0, Shape* _p1) : ks(_ks), kd(_kd), l(_l), p0(_p0), p1(_p1) {}

	void applyForces(Eigen::VectorXd& f);
	void applyJacobians(std::vector<T> &triplesJx, std::vector<T> &triplesJv);
};

#endif