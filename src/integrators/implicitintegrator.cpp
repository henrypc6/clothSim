#include "implicitintegrator.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "tools/timer.hpp"
#include <iostream>

void ImplicitIntegrator::step(double dt) {
	int d = dynamics->getDOFs();
	Eigen::VectorXd q(d), dq(d);
	dynamics->getState(q, dq);
	Eigen::VectorXd f(d);
	dynamics->getForces(f);
	SpMat M, Jx, Jv;
	Timer t;
	t.tick();
	dynamics->getInertia(M);
	t.tock();
	// std::cout << "getInertia:" << t.last << std::endl;
	dynamics->getJacobians(Jx, Jv);
	t.tock();
	// std::cout << "getJacobians:" << t.last << std::endl;
	
	SpMat A = M/dt - Jx*dt - Jv;
	Eigen::VectorXd b = f - Jv*dq + M*dq/dt;
	Eigen::ConjugateGradient< SpMat > solver;
	solver.setTolerance(1e-3);

	dq = solver.compute(A).solve(b);
	t.tock();
	// std::cout << "solve:" << t.last << std::endl;
	q += dq*dt;
	dynamics->setState(q, dq);
}