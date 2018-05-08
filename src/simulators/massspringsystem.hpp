#ifndef SIMULATORS_MASSSPRINGSYSTEM_HPP
#define SIMULATORS_MASSSPRINGSYSTEM_HPP

#include "simulator.hpp"
#include "primitives/trimesh.hpp"
#include "bvh/bvh.hpp"

#define N 35
#define L 0.8

class MassSpringSystem : public Simulator
{
public:
	TriMesh* mesh;
	TriMesh* obstacle;
	bvh *mbvh;
	double ks;
	double kd;

	MassSpringSystem(Integrator* _integrator);
	void init();
	std::vector<Drawable*> getObjectDrawables();
	std::vector<Drawable*> getObstacleDrawables();
	DrawableData getObjectData();
	DrawableData getObstacleData();

	int getDOFs();
	void getState(Eigen::VectorXd &q, Eigen::VectorXd &dq);
	void setState(const Eigen::VectorXd &q, const Eigen::VectorXd &dq);
	void getForces(Eigen::VectorXd &f);
	void getInertia(SpMat &M);
	void getJacobians(SpMat &Jx, SpMat &Jv);

	void collisionProjection();
	void advanceStep();
	
};

#endif