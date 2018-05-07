#include "system/system.hpp"

int main(int argc, char const *argv[])
{
	System(argc, argv, SIM_TYPE::MASSSPRING_SYSTEM, INTEGRATOR_TYPE::IMPLICIT).start();
	return 0;
}