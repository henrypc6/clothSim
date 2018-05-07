#ifndef COLLISION_SDF_HPP
#define COLLISION_SDF_HPP

#include <Eigen/Sparse>

class SDF
{
public:
	virtual double signedDistance(Eigen::Vector3d x) = 0;
	virtual Eigen::Vector3d normal(Eigen::Vector3d x) = 0;
};

class PlaneSDF : public SDF
{
private:
	Eigen::Vector3d x0, n;
public:
	PlaneSDF(Eigen::Vector3d _x0, Eigen::Vector3d _n) : x0(_x0), n(_n.normalized()) {}
	~PlaneSDF(){}

	double signedDistance(Eigen::Vector3d x);
	Eigen::Vector3d normal(Eigen::Vector3d x);
	
};

#endif