#include "simulator.hpp"

Simulator::Simulator(Integrator* _integrator) : integrator(_integrator) {
	assert(_integrator);
	integrator->setDynamics(this);

}

void Simulator::advanceStep() {
	if (endTime > 0 && frame*dt > endTime) {
		exit(0);
	}
	frame++;
}

double Simulator::getStepTime() {
	return dt;
}

unsigned int Simulator::getCurrentFrame() {
	return frame;
}

const char* Simulator::getImagesDir() {
	return imagesDir;
}

#include <iostream>
void Simulator::setImagesDir(const char* _imagesDir) {
	imagesDir = _imagesDir;
}