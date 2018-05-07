#include "sdf.hpp"

double PlaneSDF::signedDistance(Eigen::Vector3d x) {
	return (x - x0).dot(n);
}

Eigen::Vector3d PlaneSDF::normal(Eigen::Vector3d x) {
	return n;
}