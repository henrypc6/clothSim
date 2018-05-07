#ifndef SIMULATORS_PARTICLESYSTEM_HPP
#define SIMULATORS_PARTICLESYSTEM_HPP

#include <vector>
#include "simulator.hpp"
#include "shapes/particle.hpp"
#include "display/drawable.hpp"

class Generator
{
private:
	Eigen::Vector3d pos;
	Eigen::Vector3d dir;
	Eigen::Vector3d varVel;
	double meanVel;
public:
	Generator(Eigen::Vector3d _pos, Eigen::Vector3d _dir, Eigen::Vector3d _varVel, double _meanVel) : pos(_pos), dir(_dir), varVel(_varVel), meanVel(_meanVel) {}
	Particle* generateParticle();
	
};

class ParticleSystem : public Simulator, public Drawable
{
private:
	int n;
	std::vector<Particle*> particles;
	std::vector<Shape*> obstacles;

public:
	ParticleSystem(Integrator* _integrator);
	~ParticleSystem(){}

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

	// void draw();
	void setupObstacles();
	void setupFire();
	void generateParticles(int num);
	void addParticle(Particle* particle);
	void resetColor();
	void updateLifetime();
	void fireworkExplode();
	void generateExplosion(Eigen::Vector3d x, int num);
};

#endif