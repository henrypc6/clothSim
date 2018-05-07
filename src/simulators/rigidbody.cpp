#include "rigidbody.hpp"

RigidBody::~RigidBody() {
	for (int o = 0; o < objects.size(); o++) {
		delete objects[o];
	}
}

void RigidBody::init() {
	// TODO: modify this to loading from config file
	// shapes.push_back(Sphere());
	Sphere* sphere = new Sphere();
	sphere->setPosition(Eigen::Vector3d(-5, 8, 0));
	sphere->setColor(Eigen::Vector3d(1, 0, 0));
	objects.push_back(sphere);
	sphere->setMass(1);
	sphere->setIndex(0);
	Gravity* gravity = new Gravity(sphere);
	forces.push_back(gravity);
	Eigen::Vector3d v0(-10, 0, 5), v1(0, -3, 5), v2(0, -3, -5);
	obstacles.push_back(new Plane(v0, v1, v2));
	v0 << 0, -3, 5;
	v1 << 10, -3, 5;
	v2 << 10, -3, -5;
	obstacles.push_back(new Plane(v0, v1, v2));
	v0 << -10, -3, 5;
	v1 << -10, -3, -5;
	v2 << -10, 7, -5;
	obstacles.push_back(new Plane(v0, v1, v2));
	v0 << -10, -3, 5;
	v1 << -10, 7, 5;
	v2 << 10, 7, 5;
	obstacles.push_back(new Plane(v0, v1, v2));
	v0 << 10, -3, 5;
	v1 << 10, 7, 5;
	v2 << 10, 7, -5;
	obstacles.push_back(new Plane(v0, v1, v2));
	v0 << 10, -3, -5;
	v1 << 10, 7, -5;
	v2 << -10, 7, -5;
	obstacles.push_back(new Plane(v0, v1, v2));
}

std::vector<Drawable*> RigidBody::getObjectDrawables() {
	std::vector<Drawable*> drawables;
	for (int o = 0; o < objects.size(); o++) {
		drawables.push_back(objects[o]);
	}
	return drawables;
}

std::vector<Drawable*> RigidBody::getObstacleDrawables() {
	std::vector<Drawable*> drawables;
	for (int o = 0; o < obstacles.size(); o++) {
		drawables.push_back(obstacles[o]);
	}
	return drawables;
}

DrawableData RigidBody::getObjectData() {

}

DrawableData RigidBody::getObstacleData() {

}

int RigidBody::getDOFs() {
	return DIM*objects.size();	// Currently only consider the position and velocity of the objects.
}

void RigidBody::getState(Eigen::VectorXd &q, Eigen::VectorXd &dq) {
	assert(q.size() == getDOFs());
	assert(dq.size() == getDOFs());
	for (int o = 0; o < objects.size(); o++) {
		Eigen::Vector3d x;
		objects[o]->getPosition(x);
		q.segment(DIM*o, DIM) = x;
		Eigen::Vector3d v;
		objects[o]->getVelocity(v);
		dq.segment(DIM*o, DIM) = v;
	}
}

void RigidBody::setState(const Eigen::VectorXd &q, const Eigen::VectorXd &dq) {
	assert(q.size() == getDOFs());
	assert(dq.size() == getDOFs());
	for (int o = 0; o < objects.size(); o++) {
		objects[o]->setPosition(q.segment(DIM*o, DIM));
		objects[o]->setVelocity(dq.segment(DIM*o, DIM));
	}
}

void RigidBody::getForces(Eigen::VectorXd &f) {
	assert(f.size() == getDOFs());
	f.setZero();
	for (int ff = 0; ff < forces.size(); ff++) {
		forces[ff]->applyForces(f);
	}
}

void RigidBody::getInertia(SpMat &M) {
    unsigned int d = getDOFs();
	M.resize(d, d);

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(d);
	for (int o = 0; o < objects.size(); o++) {
		for (int i = 0; i < 3; i++) {
			tripletList.push_back(T(o*3 + i, o*3 + i, objects[o]->getMass()));
		}
	}
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}

void RigidBody::getJacobians(SpMat &Jx, SpMat &Jv) {

}

void RigidBody::collisionProjection() {
	for (int o0 = 0; o0 < objects.size(); o0++) {
		Shape* object = objects[o0];
		for (int o1 = 0; o1 < obstacles.size(); o1++) {
			Shape* obstacle = obstacles[o1];
			Eigen::Vector3d center;
			object->getPosition(center);
			Eigen::Vector3d n = obstacle->normal(center);
			Eigen::Vector3d x = center - DEFAULT_RADIUS*n;
			Eigen::Vector3d v;
			object->getVelocity(v);
			// Hard coded, to be modified
			double d = obstacle->signedDistance(x);
			if (d < 0) {
				center += -d*n;
				Eigen::Vector3d vn = v.dot(n)*n;
				Eigen::Vector3d vt = v - vn;
				vn *= -0.3;
				v = vt + vn;
				object->setPosition(center);
				object->setVelocity(v);
			}
		}
	}
}

void RigidBody::advanceStep() {
	integrator->step(dt);
	collisionProjection();
	Simulator::advanceStep();
}