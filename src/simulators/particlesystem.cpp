#include "particlesystem.hpp"
#include <iostream>
#include <stdlib.h>

double rand_1() {
	int ri = rand()%10000;
	double rd = ri/5000.0 - 1;
	return rd;
}

Particle* Generator::generateParticle() {

	Particle* particle;

	// double defaultLife = 4;
	// Eigen::Vector3d pt0 = velocity, pt1 = pos;
	// pt0[1] = 0;
	// pt1[1] = 0;

	// double w0 = 1 - pt0.norm()/vt.norm();
	// double w1 = 1 - pt1.norm()/2;
	// particle->lifetime = w0*defaultLife + rand_1()*defaultLife*.3;

	// Fire
	// Eigen::Vector3d velocity = meanVel*dir;
	// // Eigen::Vector3d diff;
	// // for (int i = 0; i < 3; i++) {
	// // 	diff[i] = varVel[i]*rand_1();
	// // }

	// double radius = 2*rand_1();
	// double angle = M_PI*rand_1();
	// // pos[0] = radius*cos(angle);
	// // pos[2] = radius*sin(angle);

	// Eigen::Vector3d vt = varVel;
	// vt[1] = 0;
	// radius = vt.norm()*rand_1();
	// angle = M_PI*rand_1();
	// velocity[0] += radius*cos(angle);
	// velocity[2] += radius*sin(angle);
	// velocity[1] = varVel[1]*(rand_1() + 1)/2;

	// // velocity += diff;
	// Particle* particle = new Particle(pos);
	// particle->setVelocity(velocity);
	// double defaultLife = 4;
	// Eigen::Vector3d pt0 = velocity, pt1 = pos;
	// pt0[1] = 0;
	// pt1[1] = 0;

	// double w0 = 1 - pt0.norm()/vt.norm();
	// double w1 = 1 - pt1.norm()/2;
	// particle->lifetime = w0*defaultLife + rand_1()*defaultLife*.3;
	return particle;
}

ParticleSystem::ParticleSystem(Integrator* _integrator) : Simulator(_integrator) {
	
}

void ParticleSystem::setupObstacles() {
	// Waterfall
	// Eigen::Vector3d v0(-10, 0, 5), v1(0, -3, 5), v2(0, -3, -5);
	// obstacles.push_back(new Plane(v0, v1, v2));
	// v0 << 0, -3, 5;
	// v1 << 10, -3, 5;
	// v2 << 10, -3, -5;
	// obstacles.push_back(new Plane(v0, v1, v2));
	// v0 << -10, 0, 5;
	// v1 << -10, 0, -5;
	// v2 << -10, 10, -5;
	// obstacles.push_back(new Plane(v0, v1, v2));
	// v0 << -10, -3, 5;
	// v1 << -10, 7, 5;
	// v2 << 10, 7, 5;
	// obstacles.push_back(new Plane(v0, v1, v2));
	// v0 << 10, -3, 5;
	// v1 << 10, 7, 5;
	// v2 << 10, 7, -5;
	// obstacles.push_back(new Plane(v0, v1, v2));
	// v0 << 10, -3, -5;
	// v1 << 10, 7, -5;
	// v2 << -10, 7, -5;
	// obstacles.push_back(new Plane(v0, v1, v2));

	// Fire
	Eigen::Vector3d v0(-100, -10, 100), v1(100, -10, 100), v2(100, -10, -100);
	obstacles.push_back(new Plane(v0, v1, v2));
}

void ParticleSystem::setupFire() {

}

void ParticleSystem::generateParticles(int num) {
	// waterfall
	// for (int i = 0; i < num; i++) {
	// 	Particle* particle = Generator(Eigen::Vector3d(0, 10, 0), Eigen::Vector3d(-10, 0, 0), Eigen::Vector3d(2, 2, 2), 1).generateParticle();
	// 	Gravity* gravity = new Gravity(particle);
	// 	forces.push_back(gravity);
	// 	addParticle(particle);
	// 	particle->setIndex(particles.size() - 1);
	// }

	// fire
	// for (int i = 0; i < num; i++) {
	// 	Particle* particle = Generator(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 20, 0), Eigen::Vector3d(2, 2, 2), 1).generateParticle();
	// 	Gravity* gravity = new Gravity(particle);
	// 	// gravity->g = Eigen::Vector3d(0, 1, 0);
	// 	forces.push_back(gravity);
	// 	addParticle(particle);
	// 	particle->setIndex(particles.size() - 1);
	// }

	// Fireworks
	// Eigen::Vector3d velocity = Eigen::Vector3d(0, 20, 0);
	// Eigen::Vector3d pos;
	// double l = 4, h = 0.5;
	// pos[0] = l*rand_1();
	// pos[2] = l*rand_1();
	// pos[1] = -10 + h*(1 + rand_1())/2;

	// Eigen::Vector3d color;
	// int idx = rand()%3;
	// color[idx] = 0;
	// color[(idx + 1)%3] = (rand_1() + 1)/4 + 0.5;
	// color[(idx + 2)%3] = (rand_1() + 1)/4 + 0.5;

	// for (int i = 0; i < 50; i++) {
	// 	Eigen::Vector3d x;
	// 	double theta = M_PI*rand_1(), phi = M_PI*rand_1();
	// 	x[0] = pos[0] + .25*cos(phi)*cos(theta);
	// 	x[1] = pos[1] + .25*sin(phi);
	// 	x[2] = pos[2] + .25*cos(phi)*sin(theta);
	// 	Particle* particle = new Particle(x);
	// 	particle->setVelocity(velocity);

	// 	double w = 1 - (x - pos).norm()/0.25;
	// 	w = w < 0 ? 0 : w;

	// 	particle->originalColor = ((1 + w)*color);
	// 	particle->lifetime = 2.2;
	// 	Gravity* gravity = new Gravity(particle);
	// 	forces.push_back(gravity);
	// 	addParticle(particle);
	// 	particle->setIndex(particles.size() - 1);
	// }

	// Snow
	for (int i = 0; i < num; i++) {

		Eigen::Vector3d x;
		x[0] = 100*rand_1();
		x[2] = 100*rand_1();
		x[1] = 50 + 1*rand_1();
		Particle* particle = new Particle(x);

		Eigen::Vector3d v;
		v[0] = 10 + rand_1();
		v[1] = rand_1();
		v[2] = rand_1();
		particle->setVelocity(v);

		Eigen::Vector3d c;
		c[0] = 0.05*(1 + rand_1());
		c[1] = 0.05*(1 + rand_1());
		c[2] = 0.9 + 0.1*rand_1();
		particle->originalColor = c;
		particle->lifetime = 7;
		Gravity* gravity = new Gravity(particle);
		forces.push_back(gravity);
		addParticle(particle);
		particle->setIndex(particles.size() - 1);
	}
}

void ParticleSystem::init() {
	// TODO: modify this to loading from config file
	// Eigen::Vector3d o(-1, 5, -1);
	// double step = 0.2;
	// for (int i = 0; i < 10; i++) {
	// 	for (int j = 0; j < 10; j++) {
	// 		for (int k = 0; k < 10; k++) {
	// 			Eigen::Vector3d x = o + Eigen::Vector3d(i*step, j*step, k*step);
	// 			Particle *p = new Particle(x);
	// 			p->setVelocity(Eigen::Vector3d(10, 1, 1));
	// 			p->setMass(1);
	// 			p->setIndex(i*100 + j*10 + k);
	// 			Gravity* gravity = new Gravity(p);
	// 			forces.push_back(gravity);
	// 			addParticle(p);
	// 		}
	// 	}
	// }
	// generateParticles(1);

	setupObstacles();
}

std::vector<Drawable*> ParticleSystem::getObjectDrawables() {
	std::vector<Drawable*> drawables;
	for (int p = 0; p < particles.size(); p++) {
		drawables.push_back(particles[p]);
	}
	return drawables;
}

std::vector<Drawable*> ParticleSystem::getObstacleDrawables() {
	std::vector<Drawable*> drawables;
	for (int o = 0; o < obstacles.size(); o++) {
		drawables.push_back(obstacles[o]);
	}
	return drawables;
}

DrawableData ParticleSystem::getObjectData() {
	DrawableData objData;
	objData.vSize = particles.size();
	objData.vertices = new float[objData.vSize*3];
	objData.color = new float[objData.vSize*3];
	objData.alpha = new float[objData.vSize];
	objData.normal = new float[objData.vSize*3];
	float maxAlpha = 0.8, minAlpha = 0.3;
	double maxLifetime = 5;
	for (int p = 0; p < particles.size(); p++) {
		Eigen::Vector3d x;
		particles[p]->getPosition(x);
		objData.vertices[p*3] = x[0];
		objData.vertices[p*3 + 1] = x[1];
		objData.vertices[p*3 + 2] = x[2];
		Eigen::Vector3d c;
		particles[p]->getColor(c);
		objData.color[p*3] = c[0];
		objData.color[p*3 + 1] = c[1];
		objData.color[p*3 + 2] = c[2];
		double lifetime = particles[p]->lifetime;
		double w = lifetime >= maxLifetime ? 1 : lifetime/maxLifetime;
		objData.alpha[p] = w*maxAlpha + (1 - w)*minAlpha;
	}
	return objData;
}

DrawableData ParticleSystem::getObstacleData() {
	DrawableData obsData;
	obsData.vSize = obstacles.size()*6;
	obsData.vertices = new float[obsData.vSize*3];
	obsData.color = new float[obsData.vSize*3];
	obsData.alpha = new float[obsData.vSize];
	double defaultAlpha = 0.6;
	for (int o = 0; o < obstacles.size(); o++) {
		Plane* plane = (Plane*)obstacles[o];
		for (int i = 0; i < 3; i++) {
			obsData.vertices[o*18 + i] = plane->v0[i];
			obsData.vertices[o*18 + 3 + i] = plane->v1[i];
			obsData.vertices[o*18 + 6 + i] = plane->v2[i];

			Eigen::Vector3d v3 = plane->v0 + plane->v2 - plane->v1;
			obsData.vertices[o*18 + 9 + i] = plane->v0[i];
			obsData.vertices[o*18 + 12 + i] = plane->v2[i];
			obsData.vertices[o*18 + 15 + i] = v3[i];
		}
		Eigen::Vector3d c;
		plane->getColor(c);
		for (int v = 0; v < 6; v++) {
			for (int i = 0; i < 3; i++) {
				obsData.color[o*18 + v*3 + i] = c[i];
			}
			obsData.alpha[o*6 + v] = defaultAlpha;
		}
	}
	return obsData;
}

int ParticleSystem::getDOFs() {
	return DIM*particles.size();
}

void ParticleSystem::getState(Eigen::VectorXd &q, Eigen::VectorXd &dq) {
	assert(q.size() == getDOFs());
	assert(dq.size() == getDOFs());
	for (int p = 0; p < particles.size(); p++) {
		Eigen::Vector3d x;
		particles[p]->getPosition(x);
		q.segment(DIM*p, DIM) = x;
		Eigen::Vector3d v;
		particles[p]->getVelocity(v);
		dq.segment(DIM*p, DIM) = v;
	}
}

void ParticleSystem::setState(const Eigen::VectorXd &q, const Eigen::VectorXd &dq) {
	assert(q.size() == getDOFs());
	assert(dq.size() == getDOFs());
	for (int p = 0; p < particles.size(); p++) {
		particles[p]->setPosition(q.segment(DIM*p, DIM));
		particles[p]->setVelocity(dq.segment(DIM*p, DIM));
	}
}

void ParticleSystem::getForces(Eigen::VectorXd &f) {
	assert(f.size() == getDOFs());
	f.setZero();
	for (int ff = 0; ff < forces.size(); ff++) {
		forces[ff]->applyForces(f);
	}

}

void ParticleSystem::getInertia(SpMat &M) {
    unsigned int d = getDOFs();
	M.resize(d, d);

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(d);
	for (int p = 0; p < particles.size(); p++) {
		for (int i = 0; i < 3; i++) {
			tripletList.push_back(T(p*3 + i, p*3 + i, particles[p]->getMass()));
		}
	}
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}

void ParticleSystem::getJacobians(SpMat &Jx, SpMat &Jv) {

}

void ParticleSystem::addParticle(Particle* particle) {
	particles.push_back(particle);
}

void ParticleSystem::resetColor() {
	// double max = 0, min = 1e10;
	Eigen::Vector3d slow(1, 1, 1);
	double ceiling = 50;
	for (int p = 0; p < particles.size(); p++) {
		Particle* part = particles[p];
		Eigen::Vector3d fast = part->originalColor;
		Eigen::Vector3d v;
		part->getVelocity(v);
		double velValue = v.norm();
		double w = velValue > ceiling ? 1 : velValue/ceiling;
		part->setColor(w*fast + (1 - w)*slow);
	}
}

void ParticleSystem::updateLifetime() {
	std::vector<int> toRemove;
	for (int p = 0; p < particles.size(); p++) {
		Particle* particle = particles[p];
		particle->lifetime -= dt;
		if (particle->lifetime < 0) {
			toRemove.push_back(p);
		}
	}
	for (int r = toRemove.size() - 1; r >= 0; r--) {
		particles.erase(particles.begin() + toRemove[r]);
		forces.erase(forces.begin() + toRemove[r]);
	}
	for (int p = 0; p < particles.size(); p++) {
		particles[p]->setIndex(p);
	}
}

void ParticleSystem::fireworkExplode() {
	std::vector<int> toRemove;
	std::vector<Eigen::Vector3d> xs;
	for (int p = 0; p < particles.size(); p++) {
		Particle* particle = particles[p];
		if (!particle->subParticle && particle->lifetime < 1.5*dt) {
			toRemove.push_back(p);
			xs.push_back(particle->x);
		}

	}
	for (int r = toRemove.size() - 1; r >= 0; r--) {
		particles.erase(particles.begin() + toRemove[r]);
		forces.erase(forces.begin() + toRemove[r]);
	}
	for (int i = 0; i < xs.size(); i++) {
		generateExplosion(xs[i], 50);
	}
	for (int p = 0; p < particles.size(); p++) {
		particles[p]->setIndex(p);
	}
}

void ParticleSystem::generateExplosion(Eigen::Vector3d x, int num) {
	int i = rand()%3;
	Eigen::Vector3d c;
	c[i] = 0;
	c[(i + 1)%3] = (rand_1() + 1)/8 + 0.75;
	c[(i + 2)%3] = (rand_1() + 1)/8 + 0.75;
	for (int p = 0; p < num; p++) {
		Particle* particle = new Particle(x);
		double theta = M_PI*rand_1(), phi = M_PI*rand_1();
		double mag = 2;
		Eigen::Vector3d v;
		v[0] = mag*cos(phi)*cos(theta);
		v[2] = mag*cos(phi)*sin(theta);
		v[1] = mag*sin(phi);
		particle->setVelocity(v);
		particle->subParticle = true;
		particle->originalColor = c;
		particle->color = c;
		particle->lifetime = 2;
		Gravity* gravity = new Gravity(particle);
		forces.push_back(gravity);
		addParticle(particle);
	}
}

void ParticleSystem::collisionProjection() {
	for (int p = 0; p < particles.size(); p++) {
		Particle* particle = particles[p];
		for (int o1 = 0; o1 < obstacles.size(); o1++) {
			Shape* obstacle = obstacles[o1];
			Eigen::Vector3d center;
			particle->getPosition(center);
			Eigen::Vector3d n = obstacle->normal(center);
			Eigen::Vector3d x = center - 0.1*n;
			Eigen::Vector3d v;
			particle->getVelocity(v);
			// Hard coded, to be modified
			double d = obstacle->signedDistance(x);
			if (d < 0) {
				center += -d*n;
				Eigen::Vector3d vn = v.dot(n)*n;
				Eigen::Vector3d vt = v - vn;
				vn *= -0;
				v = vt + vn;
				v = Eigen::Vector3d(0, 0, 0);	// For snow demo
				particle->setPosition(center);
				particle->setVelocity(v);
			}
		}
	}
}

void ParticleSystem::advanceStep() {
	// if (frame%50 == 0) {
		generateParticles(100);
	// }
	resetColor();
	// integrator->step(dt);
	for (int p = 0; p < particles.size(); p++) {
		Particle* particle = particles[p];
		Gravity* gravity = (Gravity*)forces[p];
		assert(gravity->shape == particle);
		particle->v += gravity->g*dt;
		particle->x += particle->v*dt;
	}
	updateLifetime();
	// fireworkExplode();
	std::cout << "size:" << particles.size() << std::endl;
	collisionProjection();
	Simulator::advanceStep();
}