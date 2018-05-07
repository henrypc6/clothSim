#ifndef INTEGRATORS_IMPLICITINTEGRATOR_HPP
#define INTEGRATORS_IMPLICITINTEGRATOR_HPP

#include "timeintegrator.hpp"

class ImplicitIntegrator : public Integrator
{
public:
	ImplicitIntegrator() {}
	~ImplicitIntegrator() {}
	void step(double dt);
	
};

#endif