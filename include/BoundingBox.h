#pragma once
#include "GL/glew.h"
#include "glm/glm.hpp"

class BoundingBox
{
public:
	//! draws a cube, using lines, with the center in (pX, pY, pZ) and dimensions dX, dY, dZ
	BoundingBox(float pX, float pY, float pZ, float dX, float dY, float dZ);
	~BoundingBox();

	glm::vec3* getPosition() { return &m_position; }
	void setPosition(glm::vec3 pos) { m_position = pos; }

	void render();

private:
	void clean();
	
	glm::vec3 m_position;
	
	GLuint m_vao;				// Vertex array object, the main handle for geometry
	int m_nverts;				// Number of vertices in the vertex array
	int m_nlines;				// Number of lines in the index array (may be zero)
	GLuint m_vertexbuffer;	// Buffer ID to bind to GL_ARRAY_BUFFER
	GLuint m_indexbuffer;		// Buffer ID to bind to GL_ELEMENT_ARRAY_BUFFER
	GLfloat* p_vertexarray;	// Vertex array on interleaved format: x y z nx ny nz s t
	GLuint* p_indexarray; // Element index array

};
