#ifndef SHAPES_PARTICLE_HPP
#define SHAPES_PARTICLE_HPP

#include <Eigen/Sparse>
#include "shape.hpp"

#define DEFAULT_LIFETIME 5

class Particle : public Shape
{
private:
public:
	double lifetime;
	bool subParticle = false;
	Eigen::Vector3d originalColor;
	Eigen::Vector3d n;
	
	Particle(Eigen::Vector3d& _x, double _lifetime = DEFAULT_LIFETIME);
	~Particle(){}
	void draw();
	
	double signedDistance(Eigen::Vector3d x);
	Eigen::Vector3d normal(Eigen::Vector3d x);
};

#endif