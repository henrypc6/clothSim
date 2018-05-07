#ifndef SYSTEM_SYSTEM_HPP
#define SYSTEM_SYSTEM_HPP

#include "parameters.hpp"
#include "simulators/simulator.hpp"
#include "integrators/timeintegrator.hpp"

class System
{
private:
	Simulator* simulator;
	Integrator* integrator;
public:
	System(int argc, char const *argv[], SIM_TYPE simType, INTEGRATOR_TYPE integratorType);
	~System();
	void init();
	void start();
};

#endif