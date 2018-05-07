#ifndef DISPLAY_CAMERA_HPP
#define DISPLAY_CAMERA_HPP

#include <Eigen/Dense>
#include <GL/freeglut.h>
#include <iostream>

#define XMIN -1
#define XMAX 1
#define YMIN -1
#define YMAX 1

class Mat4x4
{
public:

	float m[16];

	Mat4x4(){ // Default: Identity
		m[0] = 1.f;  m[4] = 0.f;  m[8]  = 0.f;  m[12] = 0.f;
		m[1] = 0.f;  m[5] = 1.f;  m[9]  = 0.f;  m[13] = 0.f;
		m[2] = 0.f;  m[6] = 0.f;  m[10] = 1.f;  m[14] = 0.f;
		m[3] = 0.f;  m[7] = 0.f;  m[11] = 0.f;  m[15] = 1.f;
	}

	void makeIdentity(){
		m[0] = 1.f;  m[4] = 0.f;  m[8]  = 0.f;  m[12] = 0.f;
		m[1] = 0.f;  m[5] = 1.f;  m[9]  = 0.f;  m[13] = 0.f;
		m[2] = 0.f;  m[6] = 0.f;  m[10] = 1.f;  m[14] = 0.f;
		m[3] = 0.f;  m[7] = 0.f;  m[11] = 0.f;  m[15] = 1.f;
	}

	void print(){
		std::cout << m[0] << ' ' <<  m[4] << ' ' <<  m[8]  << ' ' <<  m[12] << "\n";
		std::cout << m[1] << ' ' <<   m[5] << ' ' <<  m[9]  << ' ' <<   m[13] << "\n";
		std::cout << m[2] << ' ' <<   m[6] << ' ' <<  m[10] << ' ' <<   m[14] << "\n";
		std::cout << m[3] << ' ' <<   m[7] << ' ' <<  m[11] << ' ' <<   m[15] << "\n";
	}

	void makeScale(float x, float y, float z){
		makeIdentity();
		m[0] = x; m[5] = y; m[10] = z;
	}

	void makeTranslate(float x, float y, float z) {
		makeIdentity();
		m[12] = x; m[13] = y; m[14] = z;
	}

	void makeTranslate(Eigen::Vector3d v) {
		makeIdentity();
		m[12] = v[0];
		m[13] = v[1];
		m[14] = v[2];
	}

	// axis indicates which axis is this rotation with respect with:
	// 1:x 		2:y 	3:z
	// Rotating angle is counter clock wise
	void makeRotate(int axis, float theta) {
		makeIdentity();
		if (axis > 3 || axis < 1) {
			std::cerr << "Axis should be chosen in {1. x-axis, 2. y-axis, 3. z-axis}." << std::endl;
			exit(-1);
		}
		switch(axis) {
			case 1: {
				m[5] = cos(theta);	m[9] = -sin(theta);
				m[6] = sin(theta);	m[10] = cos(theta);
				break;
			}
			case 2: {
				m[0] = cos(theta);	m[8] = sin(theta);
				m[2] = -sin(theta);	m[10] = cos(theta);
				break;
			}
			case 3: {
				m[0] = cos(theta);	m[4] = -sin(theta);
				m[1] = sin(theta);	m[5] = cos(theta);
				break;
			}
		}
	}

	const Mat4x4 operator*(const Mat4x4 &B) {
		Mat4x4 R;
		for (int i = 0; i < 4; i++) {
	        for (int j = 0; j < 4; j++) {
	            R.m[j*4 + i] = 0;
	            for (int k = 0; k < 4; k++) {
	                R.m[j*4 + i] += m[k*4 + i]*B.m[j*4 + k];
	            }
	        }
	    }
	    return R;
	}

	void operator=(const Mat4x4 &A) {
		for (int i = 0; i < 16; i++) {
			m[i] = A.m[i];
		}
	}
};

static inline const Eigen::Vector3d operator*(const Mat4x4 &m, const Eigen::Vector3d &v){
	Eigen::Vector3d r( m.m[0]*v[0]+m.m[4]*v[1]+m.m[8]*v[2],
		m.m[1]*v[0]+m.m[5]*v[1]+m.m[9]*v[2],
		m.m[2]*v[0]+m.m[6]*v[1]+m.m[10]*v[2] );
	return r;
}

class Camera {
private:
    // double distance;
    double fov;
    double near, far;
    int winWidth, winHeight;
    Eigen::Vector3d target;

public:
    // double latitude, longitude;
    Eigen::Vector3d pos;
    Eigen::Vector3d up;
    double scale;
    Mat4x4 model;
    Mat4x4 view;
    Mat4x4 projection;
    Camera(int width, int height);
    void setProjectionMatrix();
    void setViewMatrix();
    void resize(int width, int height);
    void rotateCamera(float horizontal, float vertical);
    void translateCamera(float displacement);
};

#endif