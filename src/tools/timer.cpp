#include "timer.hpp"
#include <thread>
#include <chrono>
using namespace std::chrono;

Timer::Timer() {
	tick();
}

void Timer::tick() {
	then = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
}

void Timer::tock() {
	double now = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	last = (now - then)/1e6;
	total += last;
	then = now;
}

void Timer::sleep(double t) {
	std::this_thread::sleep_for(std::chrono::seconds((long)t));
}