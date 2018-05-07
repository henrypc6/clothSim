#include "massspringsystem.hpp"
#include <iostream>

MassSpringSystem::MassSpringSystem(Integrator* _integrator) : Simulator(_integrator), ks(800), kd(10){

}

void MassSpringSystem::init() {
	mesh = new TriMesh;
	Eigen::Vector3d o(-(N - 1)*L*.5, 8, -(N - 1)*L*.5);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Eigen::Vector3d x = Eigen::Vector3d(j*L, 0, i*L) + o;
			Particle* particle = new Particle(x);
			mesh->vertices.push_back(particle);
			Gravity* gravity = new Gravity(particle);
			forces.push_back(gravity);
			particle->setIndex(mesh->vertices.size() - 1);
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			int idx0 = i*N + j;
			int idx1 = (i + 1)*N + j;
			int idx2 = i*N + j + 1;
			int idx3 = (i + 1)*N + j + 1;
			int idx4 = (i - 1)*N + j + 1;
			if (i < N - 1) {
				mesh->edges.push_back(Eigen::Vector2i(idx0, idx1));
				SpringForce* springForce = new SpringForce(ks, kd, L, mesh->vertices[idx0], mesh->vertices[idx1]);
				forces.push_back(springForce);
			}
			if (j < N - 1) {
				mesh->edges.push_back(Eigen::Vector2i(idx0, idx2));
				SpringForce* springForce = new SpringForce(ks, kd, L, mesh->vertices[idx0], mesh->vertices[idx2]);
				forces.push_back(springForce);
			}
			if (i < N - 1 && j < N - 1) {
				mesh->edges.push_back(Eigen::Vector2i(idx0, idx3));
				SpringForce* springForce = new SpringForce(ks, kd, sqrt(2)*L, mesh->vertices[idx0], mesh->vertices[idx3]);
				forces.push_back(springForce);
				mesh->faces.push_back(Eigen::Vector3i(idx0, idx1, idx3));
				mesh->faces.push_back(Eigen::Vector3i(idx0, idx3, idx2));
			}
			if (i > 0 && j < N - 1) {
				SpringForce* springForce = new SpringForce(ks, kd, sqrt(2)*L, mesh->vertices[idx0], mesh->vertices[idx4]);
				forces.push_back(springForce);
			}
		}
	}
}

std::vector<Drawable*> MassSpringSystem::getObjectDrawables() {

}

std::vector<Drawable*> MassSpringSystem::getObstacleDrawables() {

}

DrawableData MassSpringSystem::getObjectData() {
	DrawableData objData;
	objData.vSize = mesh->vertices.size();
	objData.fSize = mesh->faces.size();
	objData.vertices = new float[objData.vSize*3];
	objData.faces = new int[objData.fSize*3];
	objData.color = new float[objData.vSize*3];
	objData.alpha = new float[objData.vSize];
	objData.normal = new float[objData.vSize*3];
	float alpha = 1;
	int* adjfNum = new int[objData.vSize];
	for (int v = 0; v < mesh->vertices.size(); v++) {
		mesh->vertices[v]->n = Eigen::Vector3d(0, 0, 0);
		adjfNum[v] = 0;
	}
	for (int f = 0; f < mesh->faces.size(); f++) {
			int idx0 = mesh->faces[f][0];
			int idx1 = mesh->faces[f][1];
			int idx2 = mesh->faces[f][2];
			objData.faces[f*3] = idx0;
			objData.faces[f*3 + 1] = idx1;
			objData.faces[f*3 + 2] = idx2;

			Eigen::Vector3d x0 = mesh->vertices[idx0]->x;
			Eigen::Vector3d x1 = mesh->vertices[idx1]->x;
			Eigen::Vector3d x2 = mesh->vertices[idx2]->x;

			Eigen::Vector3d e0 = x1 - x0;
			Eigen::Vector3d e1 = x2 - x0;
			Eigen::Vector3d n = e0.cross(e1).normalized();

			mesh->vertices[idx0]->n += n;
			adjfNum[idx0]++;
			mesh->vertices[idx1]->n += n;
			adjfNum[idx1]++;
			mesh->vertices[idx2]->n += n;
			adjfNum[idx2]++;
	}
	for (int v = 0; v < mesh->vertices.size(); v++) {
		Eigen::Vector3d x = mesh->vertices[v]->x;
		objData.vertices[v*3] = x[0];
		objData.vertices[v*3 + 1] = x[1];
		objData.vertices[v*3 + 2] = x[2];
		Eigen::Vector3d c(0, 0, 1);
		objData.color[v*3] = c[0];
		objData.color[v*3 + 1] = c[1];
		objData.color[v*3 + 2] = c[2];
		objData.alpha[v] = alpha;
		Eigen::Vector3d n = mesh->vertices[v]->n/adjfNum[v];
		objData.normal[v*3] = n[0];
		objData.normal[v*3 + 1] = n[1];
		objData.normal[v*3 + 2] = n[2];
	}
	return objData;
}

DrawableData MassSpringSystem::getObstacleData() {

}

int MassSpringSystem::getDOFs() {
	return DIM*mesh->vertices.size();
}

void MassSpringSystem::getState(Eigen::VectorXd &q, Eigen::VectorXd &dq) {
	assert(q.size() == getDOFs());
	assert(dq.size() == getDOFs());
	for (int v = 0; v < mesh->vertices.size(); v++) {
		Eigen::Vector3d x;
		mesh->vertices[v]->getPosition(x);
		q.segment(DIM*v, DIM) = x;
		Eigen::Vector3d vel;
		mesh->vertices[v]->getVelocity(vel);
		dq.segment(DIM*v, DIM) = vel;
	}
}

void MassSpringSystem::setState(const Eigen::VectorXd &q, const Eigen::VectorXd &dq) {
	assert(q.size() == getDOFs());
	assert(dq.size() == getDOFs());
	for (int v = 0; v < mesh->vertices.size(); v++) {
		mesh->vertices[v]->setPosition(q.segment(DIM*v, DIM));
		mesh->vertices[v]->setVelocity(dq.segment(DIM*v, DIM));
	}
}

void MassSpringSystem::getForces(Eigen::VectorXd &f) {
	assert(f.size() == getDOFs());
	f.setZero();
	for (int ff = 0; ff < forces.size(); ff++) {
		forces[ff]->applyForces(f);
	}
}

void MassSpringSystem::getInertia(SpMat &M) {
	unsigned int d = getDOFs();
	M.resize(d, d);

	std::vector<T> tripletList;
	tripletList.reserve(d);
	for (int v = 0; v < mesh->vertices.size(); v++) {
		for (int i = 0; i < 3; i++) {
			tripletList.push_back(T(v*3 + i, v*3 + i, mesh->vertices[v]->getMass()));
		}
	}
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}

void addBlock(SpMat& M, int row, int col, Eigen::Matrix3d b) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			M.coeffRef(row*3 + i, col*3 + j) += b(i, j);
		}
	}
}

void MassSpringSystem::getJacobians(SpMat &Jx, SpMat &Jv) {
	unsigned int d = getDOFs();
	Jx.resize(d, d);
	Jv.resize(d, d);

	// int springNum = (forces.size() - mesh->vertices.size());
	std::vector<T> triplesJx;
	std::vector<T> triplesJv;
	// triplesJx.reserve(springNum);
	// triplesJv.reserve(springNum);

	for (int f = mesh->vertices.size(); f < forces.size(); f++) {
		forces[f]->applyJacobians(triplesJx, triplesJv);
	}
	Jx.setFromTriplets(triplesJx.begin(), triplesJx.end());
	Jv.setFromTriplets(triplesJv.begin(), triplesJv.end());
}

void MassSpringSystem::collisionProjection() {

}

void MassSpringSystem::advanceStep() {
	integrator->step(dt);
	for (int i = 0; i < N; i++) {
		int idx = i*N + N - 1;
		mesh->vertices[idx]->x = Eigen::Vector3d((N - 1)*L*.5, 8, (i - (N - 1)*.5)*L);
		mesh->vertices[idx]->v = Eigen::Vector3d(0, 0, 0);
	}
	collisionProjection();
	Simulator::advanceStep();
}