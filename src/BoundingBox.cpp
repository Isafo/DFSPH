#include "BoundingBox.h"


BoundingBox::BoundingBox(float pX, float pY, float pZ, float dX, float dY, float dZ) {
	m_position = glm::vec3(pX, pY, pZ);

	m_nverts = 8;
	m_nlines = 12;

	p_vertexarray = new GLfloat[m_nverts * 3];
	p_indexarray = new GLuint[m_nlines * 2];

	//	x-coordinate				y-coordinate					z-coordinate
	//front
	p_vertexarray[0] = dX / 2;	p_vertexarray[1] = -dY / 2;		p_vertexarray[2] = dZ / 2;		// 0
	p_vertexarray[3] = dX / 2;	p_vertexarray[4] = dY / 2;		p_vertexarray[5] = dZ / 2;		// 1
	p_vertexarray[6] = -dX / 2;	p_vertexarray[7] = dY / 2;		p_vertexarray[8] = dZ / 2;		// 2
	p_vertexarray[9] = -dX / 2;	p_vertexarray[10] = -dY / 2;	p_vertexarray[11] = dZ / 2;		// 3

	// back
	p_vertexarray[12] = dX / 2;		p_vertexarray[13] = -dY / 2;	p_vertexarray[14] = -dZ / 2;		// 4
	p_vertexarray[15] = dX / 2;		p_vertexarray[16] = dY / 2;		p_vertexarray[17] = -dZ / 2;		// 5
	p_vertexarray[18] = -dX / 2;	p_vertexarray[19] = dY / 2;		p_vertexarray[20] = -dZ / 2;		// 6
	p_vertexarray[21] = -dX / 2;	p_vertexarray[22] = -dY / 2;	p_vertexarray[23] = -dZ / 2;		// 7

	// front lines
	p_indexarray[0] = 0; p_indexarray[1] = 1;
	p_indexarray[2] = 1; p_indexarray[3] = 2;
	p_indexarray[4] = 2; p_indexarray[5] = 3;
	p_indexarray[6] = 3; p_indexarray[7] = 0;

	// right side
	p_indexarray[8] = 0; p_indexarray[9] = 4;
	p_indexarray[10] = 1; p_indexarray[11] = 5;

	// right side
	p_indexarray[12] = 2; p_indexarray[13] = 6;
	p_indexarray[14] = 3; p_indexarray[15] = 7;

	// back
	p_indexarray[16] = 4; p_indexarray[17] = 5;
	p_indexarray[18] = 5; p_indexarray[19] = 6;
	p_indexarray[20] = 6; p_indexarray[21] = 7;
	p_indexarray[22] = 7; p_indexarray[23] = 4;

	// Generate one vertex array object (VAO) and bind it
	glGenVertexArrays(1, &(m_vao));
	glBindVertexArray(m_vao);

	// Generate two buffer IDs
	glGenBuffers(1, &m_vertexbuffer);
	glGenBuffers(1, &m_indexbuffer);

	// Activate the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, m_vertexbuffer);
	// Present our vertex coordinates to OpenGL
	glBufferData(GL_ARRAY_BUFFER,
		3 * m_nverts * sizeof(GLfloat), p_vertexarray, GL_STATIC_DRAW);

	// Specify how many attribute arrays we have in our VAO
	glEnableVertexAttribArray(0); // Vertex coordinates

	// Specify how OpenGL should interpret the vertex buffer data:
	// Attributes 0, 1, 2 (must match the lines above and the layout in the shader)
	// Number of dimensions (3 means vec3 in the shader, 2 means vec2)
	// Type GL_FLOAT
	// Not normalized (GL_FALSE)
	// Stride 8 floats (interleaved array with 8 floats per vertex)
	// Array buffer offset 0, 3 or 6 floats (offset into first vertex)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
		3 * sizeof(GLfloat), (void*)0);

	// Activate the index buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_indexbuffer);
	// Present our vertex indices to OpenGL
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		2 * m_nlines * sizeof(GLuint), p_indexarray, GL_STATIC_DRAW);

	// Deactivate (unbind) the VAO and the buffers again.
	// Do NOT unbind the index buffer while the VAO is still bound.
	// The index buffer is an essential part of the VAO state.
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


BoundingBox::~BoundingBox() {
	delete[] p_vertexarray;
	delete[] p_indexarray;
	clean();
}

void BoundingBox::clean() {
	if (glIsVertexArray(m_vao)) {
		glDeleteVertexArrays(1, &m_vao);
	}
	m_vao = 0;

	if (glIsBuffer(m_vertexbuffer)) {
		glDeleteBuffers(1, &m_vertexbuffer);
	}
	m_vertexbuffer = 0;

	if (glIsBuffer(m_indexbuffer)) {
		glDeleteBuffers(1, &m_indexbuffer);
	}
	m_indexbuffer = 0;
}

void BoundingBox::render() {
	glBindVertexArray(m_vao);
	glDrawElements(GL_LINES, 2 * m_nlines, GL_UNSIGNED_INT, (void*)0);
	// (mode, vertex count, type, element array buffer offset)
	glBindVertexArray(0);
}