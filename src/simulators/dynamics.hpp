#ifndef SIMULATORS_DYNAMICS_HPP
#define SIMULATORS_DYNAMICS_HPP

#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
class Dynamics
{
public:
	virtual int getDOFs() = 0;
	virtual void getState(Eigen::VectorXd &q, Eigen::VectorXd &dq) = 0;
	virtual void setState(const Eigen::VectorXd &q, const Eigen::VectorXd &dq) = 0;
	virtual void getForces(Eigen::VectorXd &f) = 0;
	virtual void getInertia(SpMat &M) = 0;
	virtual void getJacobians(SpMat &Jx, SpMat &Jv) = 0;
	
};

#endif