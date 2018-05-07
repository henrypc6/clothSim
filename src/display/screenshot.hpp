#ifndef DISPLAY_SCREENSHOT_HPP
#define DISPLAY_SCREENSHOT_HPP

class Screenshot
{
public:
	Screenshot() {}
	~Screenshot() {}

	void saveScreenshot(const int window, const unsigned int frame, const char* imagesDir);
	
};

#endif