#include "screenshot.hpp"
#include <GL/freeglut.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <png.h>

void Screenshot::saveScreenshot(const int window, const unsigned int frame, const char* imagesDir) {
	assert(imagesDir);
    glutSetWindow(window);
	int w = glutGet(GLUT_WINDOW_WIDTH);
	int h = glutGet(GLUT_WINDOW_HEIGHT);
	unsigned char *pixels = new unsigned char[w*h*3];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	for (int j = 0; j < h/2; j++)
        for (int i = 0; i < w; i++)
            for (int c = 0; c < 3; c++) {
            	int p0 = (i+w*j)*3+c;
            	int p1 = (i+w*(h-1-j))*3+c;
            	unsigned char temp = pixels[p0];
            	pixels[p0] = pixels[p1];
            	pixels[p1] = temp;
            }

	std::stringstream ss;
	ss << imagesDir << "/" << std::setw(6) << std::setfill('0') << frame << ".png";
	FILE* file = fopen(ss.str().c_str(), "wb");
	if (!file) {
		std::cerr << "#Error: could not create images file " << ss.str() << std::endl;
		std::cerr << "\tProbably the folder " << imagesDir << " doesn't exist." << std::endl;
		exit(-1);
	}
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL,
                                                  NULL, NULL);
    if (!png_ptr) {
        std::cerr << "Couldn't create a PNG write structure in file" << ss.str() << std::endl;
        fclose(file);
        return;
    }
	png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        std::cerr << "Couldn't create a PNG info structure in file " << ss.str() << std::endl;
        png_destroy_write_struct(&png_ptr, NULL);
        fclose(file);
        return;
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
        std::cout << "Had a problem writing " << ss.str() << std::endl;
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(file);
        return;
    }
    png_init_io(png_ptr, file);
    png_set_IHDR(png_ptr, info_ptr, w, h, 8,
                 false ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    int channels = 3;
    png_bytep* row_pointers = (png_bytep*) new unsigned char*[h];
    for (int y = 0; y < h; y++)
        row_pointers[y] = (png_bytep) &pixels[y*w*channels];
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    delete[] row_pointers;
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(file);
	
	delete[] pixels;
}