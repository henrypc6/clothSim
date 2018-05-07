#ifndef TOOLS_TIMER_HPP
#define TOOLS_TIMER_HPP

class Timer
{
public:
	Timer();
	double then;
	double last;
	double total;
	void tick();
	void tock();
	void sleep(double t);
};
#endif