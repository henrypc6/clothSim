#include "interface.hpp"
#include "camera.hpp"
#include "screenshot.hpp"
#include "shader.hpp"
#include "tools/timer.hpp"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <iostream>

struct MouseInfo {
    bool leftDragging, middleDraggin;    // is mouse being dragged?
    double xdrag, ydrag; // where the drag began
    double xprev, yprev; // where the mouse was last seen
};


static Camera* camera;
static MouseInfo mouseInfo;
static Simulator* simulator;
static bool play;
static int mainWindow;

static mcl::Shader shader;

namespace Global {
    GLuint vboParticleHandles[5], vboObstacleHandles[3];
    GLuint vaoParticleHandle, vaoObstacleHandle;
    GLuint posParticleBufferHandle, texCoordsBufferHandle,
        colorParticleBufferHandle, alphaParticleBufferHandle,
        normalParticleBufferHandle;
    GLuint faceBufferHandle;
    GLuint posObstacleBufferHandle, colorObstacleBufferHandle,
        alphaObstacleBufferHandle;

    DrawableData drawData;
    DrawableData obstacleData;

    int width, height, nrChannels;
    unsigned char *data;

}

Interface::Interface(int w, int h) {
    int argc = 1;
    char argv0[] = "";
    char *argv = argv0;
    glutInit(&argc, &argv);
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH|GLUT_MULTISAMPLE);
    glutInitWindowSize(w, h);
    ::mainWindow = glutCreateWindow("Animplan");
    ::camera = new Camera(w, h);

    glewExperimental = GL_TRUE;
    glewInit();
    std::stringstream ss;
    ss << "../shader/shader.";
    shader.init_from_files(ss.str() + "vert", ss.str() + "frag");
    init();

    // glEnable(GL_DEPTH_TEST);
    // glClearColor(0.9,0.7,0, 0);

    ::shader.enable();

    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Global::)

    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
}

Interface::~Interface() {
    glBindVertexArray(0);
    ::shader.disable();
    delete ::camera;
    if (::Global::obstacleData.vertices) delete[] ::Global::obstacleData.vertices;
    if (::Global::obstacleData.faces) delete[] ::Global::obstacleData.faces;
    if (::Global::obstacleData.color) delete[] ::Global::obstacleData.color;
    if (::Global::obstacleData.alpha) delete[] ::Global::obstacleData.alpha;
    if (::Global::obstacleData.normal) delete[] ::Global::obstacleData.normal;
    if (::Global::data) delete[] ::Global::data;
}

void Interface::init() {
    ::camera->setViewMatrix();
}

static Timer timer;
void Interface::run(Simulator* simulator) {
    ::simulator = simulator;
    play = false;

    using namespace Global;

    glGenBuffers(5, vboParticleHandles);
    posParticleBufferHandle = vboParticleHandles[0];
    texCoordsBufferHandle = vboParticleHandles[1];
    colorParticleBufferHandle = vboParticleHandles[2];
    alphaParticleBufferHandle = vboParticleHandles[3];
    normalParticleBufferHandle = vboParticleHandles[4];

    glGenBuffers(1, &faceBufferHandle);

    glGenVertexArrays(1, &vaoParticleHandle);
    glBindVertexArray(vaoParticleHandle);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, posParticleBufferHandle);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, texCoordsBufferHandle);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, colorParticleBufferHandle);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL);

    glEnableVertexAttribArray(3);
    glBindBuffer(GL_ARRAY_BUFFER, alphaParticleBufferHandle);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL);

    glEnableVertexAttribArray(4);
    glBindBuffer(GL_ARRAY_BUFFER, normalParticleBufferHandle);
    glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL);

    // data = stbi_load("../parks_rec_photo.0.jpg", &width, &height, &nrChannels, 0);
    data = stbi_load("../cloth.png", &width, &height, &nrChannels, 0);
    glActiveTexture(GL_TEXTURE0);
    GLuint tid;
    glGenTextures(1, &tid);
    glBindTexture(GL_TEXTURE_2D, tid);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0,
                GL_RGB, GL_UNSIGNED_BYTE, data);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glUniform1i(shader.uniform("Tex1"), 0);

    int N = 41;
    float* texCoords = new float[N*N*2];
    float L = 1.0/(N - 1);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            texCoords[2*(i*N + j)] = j*L;
            texCoords[2*(i*N + j) + 1] = i*L;
        }
    }

    glBindBuffer(GL_ARRAY_BUFFER, texCoordsBufferHandle);
    glBufferData(GL_ARRAY_BUFFER, N*N*2*sizeof(float), texCoords, GL_STATIC_DRAW);

    delete[] texCoords;

    // obstacleData = ::simulator->getObstacleData();
    // glGenBuffers(3, vboObstacleHandles);
    // posObstacleBufferHandle = vboObstacleHandles[0];
    // colorObstacleBufferHandle = vboObstacleHandles[1];
    // alphaObstacleBufferHandle = vboObstacleHandles[2];

    // glBindBuffer(GL_ARRAY_BUFFER, posObstacleBufferHandle);
    // glBufferData(GL_ARRAY_BUFFER, obstacleData.size*3*sizeof(float), obstacleData.data, GL_STATIC_DRAW);

    // glBindBuffer(GL_ARRAY_BUFFER, colorObstacleBufferHandle);
    // glBufferData(GL_ARRAY_BUFFER, obstacleData.size*3*sizeof(float), obstacleData.color, GL_STATIC_DRAW);

    // glBindBuffer(GL_ARRAY_BUFFER, alphaObstacleBufferHandle);
    // glBufferData(GL_ARRAY_BUFFER, obstacleData.size*sizeof(float), obstacleData.alpha, GL_STATIC_DRAW);

    // glGenVertexArrays(1, &vaoObstacleHandle);
    // glBindVertexArray(vaoObstacleHandle);

    // glEnableVertexAttribArray(0);
    // glBindBuffer(GL_ARRAY_BUFFER, posObstacleBufferHandle);
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL);

    // glEnableVertexAttribArray(1);
    // glBindBuffer(GL_ARRAY_BUFFER, colorObstacleBufferHandle);
    // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL);

    // glEnableVertexAttribArray(2);
    // glBindBuffer(GL_ARRAY_BUFFER, alphaObstacleBufferHandle);
    // glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL);

    timer.tick();
    glutMainLoop();
}

void idle() {
    if (play) {
        ::simulator->advanceStep();
        glutPostRedisplay();
        timer.tock();

        // long currentTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        // long elapsedTime = (currentTime - lastTime)*1e3;
        double elapsedTime = timer.last;
        double waitTime = ::simulator->getStepTime() - elapsedTime;
        std::cout << "waitTime:" << waitTime << std::endl;

        if (waitTime > 0) {
            timer.sleep(waitTime);
        }
        // lastTime = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
        timer.tock();
    }
}

void display() {
    glClearColor(0.4, 0.3, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_NORMALIZE);
    glUniformMatrix4fv( shader.uniform("model"), 1, GL_FALSE, ::camera->model.m  ); // model transformation
    glUniformMatrix4fv( shader.uniform("view"), 1, GL_FALSE, ::camera->view.m  ); // viewing transformation
    glUniformMatrix4fv( shader.uniform("projection"), 1, GL_FALSE, ::camera->projection.m );
    glUniform3f( shader.uniform("eye"), ::camera->pos[0], ::camera->pos[1], ::camera->pos[2] );
    glEnable(GL_PROGRAM_POINT_SIZE);

    using namespace Global;
    drawData = ::simulator->getObjectData();

    glBindBuffer(GL_ARRAY_BUFFER, posParticleBufferHandle);
    glBufferData(GL_ARRAY_BUFFER, drawData.vSize*3*sizeof(float), drawData.vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, colorParticleBufferHandle);
    glBufferData(GL_ARRAY_BUFFER, drawData.vSize*3*sizeof(float), drawData.color, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, alphaParticleBufferHandle);
    glBufferData(GL_ARRAY_BUFFER, drawData.vSize*sizeof(float), drawData.alpha, GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, normalParticleBufferHandle);
    glBufferData(GL_ARRAY_BUFFER, drawData.vSize*3*sizeof(float), drawData.normal, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBufferHandle);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, drawData.fSize*3*sizeof(int), drawData.faces, GL_STATIC_DRAW);

    glBindVertexArray(vaoParticleHandle);

    // glDrawArrays(GL_TRIANGLES, 0, drawData.vSize );
    glDrawElements(GL_TRIANGLES, drawData.fSize*3, GL_UNSIGNED_INT, 0);

    // glBindVertexArray(vaoObstacleHandle);
    // glDrawArrays(GL_TRIANGLES, 0, obstacleData.size);

    delete[] drawData.vertices;
    delete[] drawData.faces;
    delete[] drawData.color;
    delete[] drawData.alpha;
    delete[] drawData.normal;

    glutSwapBuffers();
    if (::simulator->getImagesDir()) {
        Screenshot().saveScreenshot(::mainWindow, ::simulator->getCurrentFrame(), ::simulator->getImagesDir());
    }
}

void reshape (int w, int h) {
    glViewport(0, 0, w, h);
    camera->resize(w, h);
    camera->setProjectionMatrix();
}

void keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case 27 : // ESCAPE
            exit(0);
        case ' ' :  // SPACE
            ::play = !::play;
            break;
    }
}

void special(int key, int x, int y) {

}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            ::mouseInfo.xprev = x;
            ::mouseInfo.yprev = y;
            ::mouseInfo.leftDragging = true;
        } else if (state == GLUT_UP) {
            ::mouseInfo.leftDragging = false;
        }
    } else if (button == GLUT_RIGHT_BUTTON) {
    } else if (button == GLUT_MIDDLE_BUTTON) {
        if (state == GLUT_DOWN) {
            ::mouseInfo.xprev = x;
            ::mouseInfo.yprev = y;
            ::mouseInfo.middleDraggin = true;
        } else if (state == GLUT_UP) {
            ::mouseInfo.middleDraggin = false;
        }
    } else if (button == 3 && state == GLUT_DOWN) {    // Wheel scroll up
        // ::camera->scale *= 1.2;
        ::camera->translateCamera(1);
        glutPostRedisplay();
    } else if (button == 4 && state == GLUT_DOWN) {     // Wheel scroll down
        // ::camera->scale /= 1.2;
        ::camera->translateCamera(-1);
        glutPostRedisplay();
    }
}

void motion(int x, int y) {
    if (::mouseInfo.leftDragging) {
        ::camera->rotateCamera(-(x - ::mouseInfo.xprev)*0.2*M_PI/180, -(y - ::mouseInfo.yprev)*0.2*M_PI/180);
        ::mouseInfo.xprev = x;
        ::mouseInfo.yprev = y;
        glutPostRedisplay();
    } else if (::mouseInfo.middleDraggin) {
        ::camera->pos[0] -= (x - ::mouseInfo.xprev)*0.02;
        ::camera->pos[1] += (y - ::mouseInfo.yprev)*0.02;
        ::camera->setViewMatrix();
        ::mouseInfo.xprev = x;
        ::mouseInfo.yprev = y;
        glutPostRedisplay();
    }
}
