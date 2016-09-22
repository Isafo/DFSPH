#include "Sphere.h"

#define M_PI 3.14159265358979323846f

Sphere::Sphere(float x, float y, float z, float _rad)
{
	m_position = glm::vec3(x, y, z);
	m_radius = _rad;
	createSphere(_rad, 32);
}

void Sphere::clean() {

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


void Sphere::render()
{
	glBindVertexArray(m_vao);
	glDrawElements(GL_TRIANGLES, 3 * m_ntris, GL_UNSIGNED_INT, (void*)0);
	// (mode, vertex count, type, element array buffer offset)
	glBindVertexArray(0);
}

Sphere::~Sphere(void)
{
	clean();
}


void Sphere::createSphere(float radius, int segments) {
	int i, j, base, i0;
	float x, y, z, R;
	double theta, phi;
	int vsegs, hsegs;
	int stride = 8;

	// Delete any previous content in the TriangleSoup object
	clean();

	vsegs = segments;
	if (vsegs < 2) vsegs = 2;
	hsegs = vsegs * 2;
	m_nverts = 1 + (vsegs - 1) * (hsegs + 1) + 1; // top + middle + bottom
	m_ntris = hsegs + (vsegs - 2) * hsegs * 2 + hsegs; // top + middle + bottom
	p_vertexarray = new float[m_nverts * 8];
	p_indexarray = new GLuint[m_ntris * 3];

	// The vertex array: 3D xyz, 3D normal, 2D st (8 floats per vertex)
	// First vertex: top pole (+z is "up" in object local coords)
	p_vertexarray[0] = 0.0f;
	p_vertexarray[1] = 0.0f;
	p_vertexarray[2] = radius;
	p_vertexarray[3] = 0.0f;
	p_vertexarray[4] = 0.0f;
	p_vertexarray[5] = 1.0f;
	p_vertexarray[6] = 0.5f;
	p_vertexarray[7] = 1.0f;
	// Last vertex: bottom pole
	base = (m_nverts - 1)*stride;
	p_vertexarray[base] = 0.0f;
	p_vertexarray[base + 1] = 0.0f;
	p_vertexarray[base + 2] = -radius;
	p_vertexarray[base + 3] = 0.0f;
	p_vertexarray[base + 4] = 0.0f;
	p_vertexarray[base + 5] = -1.0f;
	p_vertexarray[base + 6] = 0.5f;
	p_vertexarray[base + 7] = 0.0f;
	// All other vertices:
	// vsegs-1 latitude rings of hsegs+1 vertices each
	// (duplicates at texture seam s=0 / s=1)
	for (j = 0; j<vsegs - 1; j++) { // vsegs-1 latitude rings of vertices
		theta = (double)(j + 1) / vsegs*M_PI;
		z = cos(theta);
		R = sin(theta);
		for (i = 0; i <= hsegs; i++) { // hsegs+1 vertices in each ring (duplicate for texcoords)
			phi = (double)i / hsegs*2.0*M_PI;
			x = R*cos(phi);
			y = R*sin(phi);
			base = (1 + j*(hsegs + 1) + i)*stride;
			p_vertexarray[base] = radius*x;
			p_vertexarray[base + 1] = radius*y;
			p_vertexarray[base + 2] = radius*z;
			p_vertexarray[base + 3] = x;
			p_vertexarray[base + 4] = y;
			p_vertexarray[base + 5] = z;
			p_vertexarray[base + 6] = (float)i / hsegs;
			p_vertexarray[base + 7] = 1.0f - (float)(j + 1) / vsegs;
		}
	}

	// The index array: triplets of integers, one for each triangle
	// Top cap
	for (i = 0; i<hsegs; i++) {
		p_indexarray[3 * i] = 0;
		p_indexarray[3 * i + 1] = 1 + i;
		p_indexarray[3 * i + 2] = 2 + i;
	}
	// Middle part (possibly empty if vsegs=2)
	for (j = 0; j<vsegs - 2; j++) {
		for (i = 0; i<hsegs; i++) {
			base = 3 * (hsegs + 2 * (j*hsegs + i));
			i0 = 1 + j*(hsegs + 1) + i;
			p_indexarray[base] = i0;
			p_indexarray[base + 1] = i0 + hsegs + 1;
			p_indexarray[base + 2] = i0 + 1;
			p_indexarray[base + 3] = i0 + 1;
			p_indexarray[base + 4] = i0 + hsegs + 1;
			p_indexarray[base + 5] = i0 + hsegs + 2;
		}
	}
	// Bottom cap
	base = 3 * (hsegs + 2 * (vsegs - 2)*hsegs);
	for (i = 0; i<hsegs; i++) {
		p_indexarray[base + 3 * i] = m_nverts - 1;
		p_indexarray[base + 3 * i + 1] = m_nverts - 2 - i;
		p_indexarray[base + 3 * i + 2] = m_nverts - 3 - i;
	}

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
		8 * m_nverts * sizeof(GLfloat), p_vertexarray, GL_STATIC_DRAW);
	// Specify how many attribute arrays we have in our VAO
	glEnableVertexAttribArray(0); // Vertex coordinates
	glEnableVertexAttribArray(1); // Normals
	glEnableVertexAttribArray(2); // Texture coordinates
	// Specify how OpenGL should interpret the vertex buffer data:
	// Attributes 0, 1, 2 (must match the lines above and the layout in the shader)
	// Number of dimensions (3 means vec3 in the shader, 2 means vec2)
	// Type GL_FLOAT
	// Not normalized (GL_FALSE)
	// Stride 8 (interleaved array with 8 floats per vertex)
	// Array buffer offset 0, 3, 6 (offset into first vertex)
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
		8 * sizeof(GLfloat), (void*)0); // xyz coordinates
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE,
		8 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat))); // normals
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE,
		8 * sizeof(GLfloat), (void*)(6 * sizeof(GLfloat))); // texcoords

	// Activate the index buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_indexbuffer);
	// Present our vertex indices to OpenGL
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		3 * m_ntris*sizeof(GLuint), p_indexarray, GL_STATIC_DRAW);
	
	// Deactivate (unbind) the VAO and the buffers again.
	// Do NOT unbind the buffers while the VAO is still bound.
	// The index buffer is an essential part of the VAO state.
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

};