#ifndef INTEGRATORS_TIMEINTEGRATOR_HPP
#define INTEGRATORS_TIMEINTEGRATOR_HPP

#include "simulators/dynamics.hpp"

class Integrator
{
private:
protected:
	Dynamics* dynamics;
public:
	Integrator(){}
	~Integrator(){}
	void setDynamics(Dynamics* _dynamics) {
		dynamics = _dynamics;
	}
	virtual void step(double dt) = 0;
};

#endif