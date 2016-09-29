#pragma once
#include "Gl/glew.h"
#include "glm/glm.hpp"

class Sphere {
public:

	// Creates a sphere  
	Sphere(float x, float y, float z, float _rad);
	~Sphere(void);

	Sphere() {
		m_vao = 0;
		m_vertexbuffer = 0;
		m_indexbuffer = 0;
		p_vertexarray = nullptr;
		p_indexarray = nullptr;
		m_nverts = 0;
		m_ntris = 0;
	};
	void setRadius(float r) { m_radius = r; }
	void createSphere(float m_radius, int m_segments);
	void clean();
	void render();

	float getRadius() const { return m_radius; }
	glm::vec3* getPosition() { return &m_position; }
	void setPosition(glm::vec3 pos) { m_position = pos; }

private:
	GLuint m_vao;          // Vertex array object, the main handle for geometry
	int m_nverts; // Number of vertices in the vertex array
	int m_ntris;  // Number of triangles in the index array (may be zero)
	GLuint m_vertexbuffer; // Buffer ID to bind to GL_ARRAY_BUFFER
	GLuint m_indexbuffer;  // Buffer ID to bind to GL_ELEMENT_ARRAY_BUFFER
	GLfloat* p_vertexarray; // Vertex array on interleaved format: x y z nx ny nz s t
	GLuint* p_indexarray;   // Element index array

	float m_radius;
	glm::vec3 m_position;

};