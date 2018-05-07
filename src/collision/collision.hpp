#ifndef COLLISION_COLLISION_HPP
#define COLLISION_COLLISION_HPP

#include "sdf.hpp"

class Collision
{
private:
	unsigned int i0, i1;
	Eigen::Vector3d n;
	double d;
public:
	Collision() : i0(-1), i1(-1), d(0) {}
	~Collision() {}
	bool isCollision();
};

#endif