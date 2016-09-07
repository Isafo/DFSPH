#pragma once

#include "GLFW/glfw3.h"
#include <vector>

struct glfwWindowInfo
{
	GLFWwindow* object;
	int width;
	int height;
};
struct glfwMonitorInfo
{
	GLFWmonitor* object;
};


class glfwContext
{
public:
	glfwContext();
	~glfwContext();

	void init(int _w, int _h, const char* _wName);

	void printGLInfo();

	void getCurrentWindow(GLFWwindow* &_window);

	GLFWwindow* currentWindow;

private:
	std::vector<glfwWindowInfo> window;
	
	std::vector<glfwMonitorInfo> monitor;
};

