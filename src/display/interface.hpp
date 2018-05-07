#ifndef DISPLAY_INTERFACE_HPP
#define DISPLAY_INTERFACE_HPP

#include <vector>
#include "drawable.hpp"
#include "simulators/simulator.hpp"

#define WIDTH 720
#define HEIGHT 720

class Interface
{
private:
public:
	Interface(int w = WIDTH, int h = HEIGHT);
	~Interface();
	// void display(std::vector<Drawable*> drawables);
	void init();
	void run(Simulator* simulator);
};

void idle();
void display();
void reshape (int w, int h);
void keyboard(unsigned char key, int x, int y);
void special(int key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
// void errorCallback(int error, const char* description);
// void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
// void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
// void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);
// void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);

#endif