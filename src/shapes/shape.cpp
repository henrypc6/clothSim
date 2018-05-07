#include "shape.hpp"

void Shape::setColor(const Eigen::Vector3d& _color) {
	color = _color;
}

void Shape::getColor(Eigen::Vector3d& _color) {
	_color = color;
}

void Shape::setPosition(const Eigen::Vector3d& _x) {
	x = _x;
}

void Shape::getPosition(Eigen::Vector3d& _x) {
	_x = x;
}

void Shape::setVelocity(const Eigen::Vector3d& _v) {
	v = _v;
}

void Shape::getVelocity(Eigen::Vector3d& _v) {
	_v = v;
}

void Shape::setIndex(const unsigned int _i) {
	i = _i;
}

unsigned int Shape::getIndex() {
	return i;
}

void Shape::setMass(const double _m) {
	m = _m;
}

double Shape::getMass() {
	return m;
}

void Sphere::draw() {
	glPushMatrix();
	glColor3f(color[0], color[1], color[2]);
	glTranslatef(x[0], x[1], x[2]);
    glutSolidSphere(r, 50, 50);
	glPopMatrix();
}

double Sphere::signedDistance(Eigen::Vector3d x) {
	return (x - this->x).norm() - r;
}

Eigen::Vector3d Sphere::normal(Eigen::Vector3d x) {
	return (x - this->x).normalized();
}

#include <iostream>
void Plane::draw() {
	Eigen::Vector3d e0 = v1 - v0;
	Eigen::Vector3d e1 = v2 - v1;
	Eigen::Vector3d n = e0.cross(e1).normalized();
	Eigen::Vector3d v3 = v0 + v2 - v1;
	glPushMatrix();
	// glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor4f(color[0], color[1], color[2], 0.6);
	glNormal3f(n[0], n[1], n[2]);
	glBegin(GL_TRIANGLES);
	glVertex3f(v0[0], v0[1], v0[2]);
	glVertex3f(v1[0], v1[1], v1[2]);
	glVertex3f(v2[0], v2[1], v2[2]);

	glVertex3f(v0[0], v0[1], v0[2]);
	glVertex3f(v2[0], v2[1], v2[2]);
	glVertex3f(v3[0], v3[1], v3[2]);
	glEnd();
	glPopMatrix();
}

double Plane::signedDistance(Eigen::Vector3d x) {

	Eigen::Vector3d e0 = v1 - v0;
	Eigen::Vector3d e1 = v2 - v1;
	Eigen::Vector3d n = e0.cross(e1).normalized();
	Eigen::Vector3d v0x = x - v0;
	v0x -= v0x.dot(n)*n;
	Eigen::Vector3d v2x = x - v2;
	v2x -= v2x.dot(n)*n;
	if (v0x.dot(v2x) > 0) {
		return 1;
	}
	return (x - v0).dot(n);
}

Eigen::Vector3d Plane::normal(Eigen::Vector3d x) {
	Eigen::Vector3d e0 = v1 - v0;
	Eigen::Vector3d e1 = v2 - v1;
	Eigen::Vector3d n = e0.cross(e1).normalized();
	return n;
}