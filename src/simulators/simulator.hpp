#ifndef SIMULATORS_SIMULATOR_HPP
#define SIMULATORS_SIMULATOR_HPP

#include "display/drawable.hpp"
#include "forces/force.hpp"
#include "dynamics.hpp"
#include "integrators/timeintegrator.hpp"
#include <vector>

#define DEFAULT_STEP_TIME 1e-2
#define DEFAULT_END_TIME -1

class Simulator : public Dynamics
{
protected:
	Integrator* integrator;
	std::vector<Force*> forces;

	double dt = DEFAULT_STEP_TIME;
	double endTime = DEFAULT_END_TIME;
	unsigned int frame = 0;

	const char* imagesDir = NULL;

public:
	Simulator(Integrator* _integrator);
	~Simulator(){}

	virtual void init()=0;
	virtual std::vector<Drawable*> getObjectDrawables()= 0;
	virtual std::vector<Drawable*> getObstacleDrawables()= 0;
	virtual DrawableData getObjectData() = 0;
	virtual DrawableData getObstacleData() = 0;

	virtual int getDOFs() = 0;
	virtual void getState(Eigen::VectorXd &q, Eigen::VectorXd &dq) = 0;
	virtual void setState(const Eigen::VectorXd &q, const Eigen::VectorXd &dq) = 0;
	virtual void getForces(Eigen::VectorXd &f) = 0;
	virtual void getInertia(SpMat &M) = 0;
	virtual void getJacobians(SpMat &Jx, SpMat &Jv) = 0;

	virtual void collisionProjection() = 0;
	virtual void advanceStep();

	double getStepTime();
	unsigned int getCurrentFrame();
	const char* getImagesDir();
	void setImagesDir(const char* _imagesDir);
};

#endif