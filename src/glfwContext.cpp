#include "glfwContext.h"
#include <iostream>


glfwContext::glfwContext()
{
}


glfwContext::~glfwContext()
{
	for (int i = 0; i < window.size(); i++)
	{
		glfwDestroyWindow(window[i].object);
	}

	glfwTerminate();
}

void glfwContext::init(int _w, int _h, const char* _wName)
{
	
	if (!glfwInit()) {
		std::cout << "ERROR: could not start GLFW3" << std::endl;
		return;
	}
	
	window.resize(1);
	monitor.resize(1);

	window[0].width = _w; window[0].height = _h;

	monitor[0].object = NULL;
	window[0].object = glfwCreateWindow(window[0].width, window[0].height, _wName, monitor[0].object, NULL);

	if (!window[0].object) {
		glfwTerminate();
		std::cout << "failed to initialize" << std::endl;
	}

	glfwMakeContextCurrent(window[0].object);
	currentWindow = window[0].object;

	glViewport(0, 0, window[0].width, window[0].height);

	std::cout << "glfw initialization complete" << std::endl;
}

void glfwContext::getCurrentWindow(GLFWwindow* &_window)
{
	_window = currentWindow;
}


void glfwContext::printGLInfo()
{

	// Print some info about the OpenGL context...
	int l_Major = glfwGetWindowAttrib(currentWindow, GLFW_CONTEXT_VERSION_MAJOR);
	int l_Minor = glfwGetWindowAttrib(currentWindow, GLFW_CONTEXT_VERSION_MINOR);
	int l_Profile = glfwGetWindowAttrib(currentWindow, GLFW_OPENGL_PROFILE);
	std::cout << "OpenGL: " << l_Major << "." << l_Minor << std::endl;
	
	if (l_Major >= 3) { // Profiles introduced in OpenGL 3.0...
		if (l_Profile == GLFW_OPENGL_COMPAT_PROFILE) 
			std::cout << "GLFW_OPENGL_COMPAT_PROFILE" << std::endl; 
		else 
			std::cout << "GLFW_OPENGL_CORE_PROFILE" << std::endl;
	}
	
	//std::cout << "Vendor: " << (char*)glGetString(GL_VENDOR) << std::endl;
	//std::cout << "Renderer: " <<  (char*)glGetString(GL_RENDERER) << std::endl;
	printf("Vendor: %s\n", (char*)glGetString(GL_VENDOR));
	printf("Renderer: %s\n", (char*)glGetString(GL_RENDERER));
}	
