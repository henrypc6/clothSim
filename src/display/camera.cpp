#include "camera.hpp"

Camera::Camera(int width, int height) : winWidth(width), winHeight(height), pos(10, 10, 18), target(0, 0, 0), up(0, 1, 0), fov(45), near(0.05), far(10), scale(1)  {
    
}

void Camera::setProjectionMatrix() {
    double w = 2*near*tan(fov*.5);
    double h = w*winHeight/winWidth;
    double right = w/2, left = -w/2;
    double top = h/2, bottom = - h/2;
    projection.m[0] = 2*near/(right - left);
    projection.m[5] = 2*near/(top - bottom);
    projection.m[8] = (right + left)/(right - left);
    projection.m[9] = (top + bottom)/(top - bottom);
    projection.m[10] = -(far + near)/(far - near);
    projection.m[11] = -1;
    projection.m[14] = -2*far*near/(far - near);
}

void Camera::setViewMatrix() {
    view.makeIdentity();
    Mat4x4 mT, mR;
    mT.makeTranslate(Eigen::Vector3d(0, 0, 0) - pos);    // Translate the camera to the origin
    Eigen::Vector3d viewDir = (target - pos).normalized();
    Eigen::Vector3d n = Eigen::Vector3d(0, 0, 0) - viewDir;    // n is negative the view direction
    n.normalize();
    Eigen::Vector3d u = up.cross(n);
    u.normalize();
    Eigen::Vector3d v = n.cross(u);
    v.normalize();

    mR.makeIdentity();
    mR.m[0] = u[0]; mR.m[4] = u[1]; mR.m[8] = u[2];
    mR.m[1] = v[0]; mR.m[5] = v[1]; mR.m[9] = v[2];
    mR.m[2] = n[0]; mR.m[6] = n[1]; mR.m[10] = n[2];

    view = mR*mT;
}

void Camera::resize(int width, int height) {
    winWidth = width;
    winHeight = height;
}

void Camera::rotateCamera(float horizontal, float vertical) {
    Mat4x4 mT1, mT2;
    Mat4x4 mR1, mR2;
    mT1.makeTranslate(-target[0], -target[1], -target[2]);
    mT2.makeTranslate(target[0], target[1], target[2]);
    mR1.makeRotate(2, horizontal);
    mR2.makeRotate(1, vertical);
    pos = mT2*mR2*(mR1*mT1*pos);

    setViewMatrix();
}

void Camera::translateCamera(float displacement) {
    Eigen::Vector3d viewDir = (target - pos).normalized();
    pos += viewDir*displacement;
    setViewMatrix();
}

// void Camera::apply() {

    // double w = glutGet(GLUT_WINDOW_WIDTH), h = glutGet(GLUT_WINDOW_HEIGHT);
    // double aspect = w/h;
    // glMatrixMode(GL_PROJECTION);
    // glLoadIdentity();
    // gluPerspective(fov, aspect, 0.01*distance, 100*distance);
    // glMatrixMode(GL_MODELVIEW);
    // glLoadIdentity();
    // gluLookAt(25, 9, 2*distance, 0, 0, 0, 0, 1, 0);
    // // glTranslatef(0, 0, -distance);
    // // glTranslatef(-target(0), -target(1), -target(2));

    // glTranslatef(pos[0], pos[1], pos[2]);
    // glRotatef(latitude, 1, 0, 0);
    // glRotatef(-longitude, 0, 1, 0);
    // glScalef(scale, scale, scale);
// }
