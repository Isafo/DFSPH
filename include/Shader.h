#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class Shader {

public:
	GLuint programID;

	Shader();
	Shader(const char *vertexFilePath, const char *fragmentFilePath);
	~Shader();

	void createShader();
	void createShader(const char *vertexFilePath, const char *fragmentFilePath);
	void createShader(const char *vertexFilePath, const char *fragmentFilePath, const char* geometryFilePath);
	void createTransformShader(const char * vertexFilePath, const char * fragmentFilePath, const char * geometryFilePath);

private:
	std::string readFile(const char *filePath);

};