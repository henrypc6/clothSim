#include "force.hpp"
#include <iostream>

void Gravity::applyForces(Eigen::VectorXd& f) {
	f.segment(DIM*shape->getIndex(), DIM) += shape->getMass()*g;
}

void Gravity::applyJacobians(std::vector<T> &triplesJx, std::vector<T> &triplesJv) {

}

void SpringForce::applyForces(Eigen::VectorXd& f) {
	Eigen::Vector3d dx = (p1->x - p0->x);
	double disX = (dx.norm() - l);
	Eigen::Vector3d dv = (p1->v - p0->v);
	Eigen::Vector3d dir = dx.normalized();
	Eigen::Vector3d force = ks*disX*dir + kd*dv.dot(dir)*dir;
	f.segment(p0->i*3, 3) += force;
	f.segment(p1->i*3, 3) -= force;
}

void SpringForce::applyJacobians(std::vector<T> &triplesJx, std::vector<T> &triplesJv) {
	Eigen::Vector3d dx = (p1->x - p0->x);
	Eigen::Vector3d dir = dx.normalized();
	Eigen::Matrix3d xxt = dir*dir.transpose();
	Eigen::Matrix3d Jxi = ks*((1 - l/dx.norm())*(Eigen::Matrix3d::Identity() - xxt) + xxt);
	Eigen::Matrix3d Jvi = kd*xxt;
	int i = p0->i, j = p1->i;
	for (int row = 0; row < 3; row++) {
		for (int col = 0; col < 3; col++) {
			triplesJx.push_back(T(i*DIM + row, i*DIM + col, -Jxi(row, col)));
			triplesJx.push_back(T(i*DIM + row, j*DIM + col, Jxi(row, col)));
			triplesJx.push_back(T(j*DIM + row, i*DIM + col, Jxi(row, col)));
			triplesJx.push_back(T(j*DIM + row, j*DIM + col, -Jxi(row, col)));

			triplesJv.push_back(T(i*DIM + row, i*DIM + col, -Jvi(row, col)));
			triplesJv.push_back(T(i*DIM + row, j*DIM + col, Jvi(row, col)));
			triplesJv.push_back(T(j*DIM + row, i*DIM + col, Jvi(row, col)));
			triplesJv.push_back(T(j*DIM + row, j*DIM + col, -Jvi(row, col)));
		}
	}
	// triplesJx.push_back(i*DIM)
	// Jx.block(i*DIM, i*DIM, DIM, DIM) += -Jxi;
	// Jx.block(i*DIM, j*DIM, DIM, DIM) += Jxi;
	// Jx.block(j*DIM, i*DIM, DIM, DIM) += Jxi;
	// Jx.block(j*DIM, j*DIM, DIM, DIM) += -Jxi;
	
	// Jv.block(i*DIM, i*DIM, DIM, DIM) += -Jvi;
	// Jv.block(i*DIM, j*DIM, DIM, DIM) += Jvi;
	// Jv.block(j*DIM, i*DIM, DIM, DIM) += Jvi;
	// Jv.block(j*DIM, j*DIM, DIM, DIM) += -Jvi;
}