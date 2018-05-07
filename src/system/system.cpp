#include "system.hpp"
#include "simulators/rigidbody.hpp"
#include "simulators/particlesystem.hpp"
#include "simulators/massspringsystem.hpp"
#include "integrators/explicitintegrator.hpp"
#include "integrators/implicitintegrator.hpp"
#include "display/interface.hpp"
#include <iostream>

System::System(int argc, char const *argv[], SIM_TYPE simType, INTEGRATOR_TYPE integratorType) {
	const char* imagesDir = NULL;
	if (argc == 2) {
		imagesDir = argv[1];
	}
	switch(integratorType) {
		case INTEGRATOR_TYPE::EXPLICIT:
			integrator = new ExplicitIntegrator();
			break;
		case INTEGRATOR_TYPE::IMPLICIT:
			integrator = new ImplicitIntegrator();
			break;
		default:
			std::cerr << "Undefined integration type " << integratorType << std::endl;
			exit(-1);
	}
	switch(simType) {
		case SIM_TYPE::RIGID_BODY:
			simulator = new RigidBody(integrator);
			break;
		case SIM_TYPE::PARTICLE_SYSTEM:
			simulator = new ParticleSystem(integrator);
			break;
		case SIM_TYPE::MASSSPRING_SYSTEM:
			simulator = new MassSpringSystem(integrator);
			break;
		default:
			std::cerr << "Undefined simulation type " << simType << std::endl;
			exit(-1);
	}
	if (imagesDir) {
		simulator->setImagesDir(imagesDir);
	}
}

System::~System() {
	delete simulator;
	delete integrator;
}

void System::init() {
	simulator->init();
}

void System::start() {
	init();
	Interface interface;
	interface.run(simulator);
}