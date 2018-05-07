#ifndef DISPLAY_DRAWABLE_HPP
#define DISPLAY_DRAWABLE_HPP

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <Eigen/Dense>
#include <vector>

class Drawable
{
public:
	Drawable(){}
	~Drawable(){}

	virtual void draw() {}
};

class DrawableData
{
private:
public:
	int vSize, fSize;
	float* vertices;
	int* faces;
	float* color;
	float* alpha;
	float* normal;
	DrawableData();
	~DrawableData();
};

#endif