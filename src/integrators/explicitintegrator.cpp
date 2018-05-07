#include "explicitintegrator.hpp"

void ExplicitIntegrator::step(double dt) {
	// int d = dynamics->getDOFs();
	// Eigen::VectorXd q(d), dq(d);
	// dynamics->getState(q, dq);
	// Eigen::VectorXd f(d);
	// dynamics->getForces(f);
	// SpMat M;
	// dynamics->getInertia(M);
	// Eigen::SparseLU<SpMat> solver;
	// solver.analyzePattern(M);
	// solver.factorize(M);
	// Eigen::VectorXd dqDelta = solver.solve(f*dt);
	// dq += dqDelta;
	// q += dq*dt;
	// dynamics->setState(q, dq);
}