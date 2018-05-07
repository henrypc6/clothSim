#ifndef SIMULATORS_RIGIDBODY_HPP
#define SIMULATORS_RIGIDBODY_HPP

#include "simulator.hpp"
#include "display/drawable.hpp"
#include "shapes/shape.hpp"

class RigidBody : public Simulator
{
private:
	std::vector<Shape*> objects;
	std::vector<Shape*> obstacles;
public:
	RigidBody(Integrator* _integrator) : Simulator(_integrator) {}
	~RigidBody();
	
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

	void advanceStep();
	void collisionProjection();
};

#endif