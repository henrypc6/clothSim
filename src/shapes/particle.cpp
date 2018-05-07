#include "particle.hpp"

Particle::Particle(Eigen::Vector3d& _x, double _lifetime) {
	x = _x;
	lifetime = _lifetime;
	setColor(Eigen::Vector3d(0.1, 0.2, 0.9));
}

void Particle::draw() {
	glPushMatrix();
	glColor3f(color[0], color[1], color[2]);
	glTranslatef(x[0], x[1], x[2]);
    glutSolidSphere(0.1, 5, 5);
    glPopMatrix();
}

double Particle::signedDistance(Eigen::Vector3d x) {
	
}

Eigen::Vector3d Particle::normal(Eigen::Vector3d x) {

}