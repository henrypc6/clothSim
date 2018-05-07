#ifndef INTEGRATORS_EXPLICITINTEGRATOR_HPP
#define INTEGRATORS_EXPLICITINTEGRATOR_HPP

#include "timeintegrator.hpp"

class ExplicitIntegrator : public Integrator
{
public:
	ExplicitIntegrator() {}
	~ExplicitIntegrator(){}
	void step(double dt);
};

#endif