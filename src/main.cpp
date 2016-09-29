#include "GL/glew.h"

#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
#include <glm/gtc/type_ptr.hpp> //glm::make:mat4
#include <glm/gtx/rotate_vector.hpp> // rotate vector

#include "glfwContext.h"
#include "Shader.h"
#include "MatrixStack.h"
#include "Camera.h"
#include "Sphere.h"
#include "BoundingBox.h"
#include "SPH.h"
#include <iostream>

void inputHandler(GLFWwindow* _window, double _dT);
void cameraHandler(GLFWwindow* _window, double _dT, Camera* _cam);
void GLcalls();



int main() {

	glfwContext glfw;
	GLFWwindow* currentWindow = nullptr;

	glfw.init(1920, 1080, "Waves4Life");
	glfw.getCurrentWindow(currentWindow);
	glfwSetCursorPos(currentWindow, 960, 540);

	//start GLEW extension handler
	glewExperimental = GL_TRUE;
	GLenum l_GlewResult = glewInit();
	if (l_GlewResult != GLEW_OK)
		std::cout << "glewInit() error." << std::endl;

	// Print some info about the OpenGL context...
	glfw.printGLInfo();

	Shader sceneLight;
	sceneLight.createShader("shaders/scene.vert", "shaders/scene.frag");

	GLint locationP = glGetUniformLocation(sceneLight.programID, "P"); //perspective matrix
	GLint locationMV = glGetUniformLocation(sceneLight.programID, "MV"); //modelview matrix
	GLint locationLP = glGetUniformLocation(sceneLight.programID, "LP"); // light position
	GLint locationTex = glGetUniformLocation(sceneLight.programID, "tex"); //texcoords

	MatrixStack MVstack; MVstack.init();

	//BoundingBox bbox(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
	BoundingBox bbox(0.f, 0.f, 0.f, 1.f, 1.f, 1.f);

	SPH s;
	s.init_positions(bbox.getPosition(), 5, 5);

	Sphere sphere(0.0f, 0.0f, 0.0f, s.get_particle_radius());

	Camera mCamera;
	// mCamera.setPosition(&glm::vec3(0f, 0.f, 0.f));
	mCamera.setPosition(&glm::vec3(0.f, 0.f, 1.5f));
	mCamera.update();

	double lastTime = glfwGetTime() - 0.001f;
	double dT = 0.0;
	while (!glfwWindowShouldClose(currentWindow))
	{
		dT = glfwGetTime() - lastTime;
		lastTime = glfwGetTime();

		//glfw input handler
		inputHandler(currentWindow, dT);
		mCamera.fpsCamera(currentWindow, dT);

		GLcalls();

		glUseProgram(sceneLight.programID);

		s.update(dT / 10);

		MVstack.push();//Camera transforms --<
		glUniformMatrix4fv(locationP, 1, GL_FALSE, mCamera.getPerspective());
		MVstack.multiply(mCamera.getTransformM());

		glm::vec3 particlePos;
		for (int i = 0; i < D_NR_OF_PARTICLES; ++i) {
			MVstack.push();
			particlePos = glm::vec3(
				s.get_particle_positions()->x[i],
				s.get_particle_positions()->y[i],
				s.get_particle_positions()->z[i]
			);
			MVstack.translate(&particlePos);
			//MVstack.translate(particle.getPosition());
			glUniformMatrix4fv(locationMV, 1, GL_FALSE, MVstack.getCurrentMatrix());
			sphere.render();
			MVstack.pop();
		}

		MVstack.push();
		MVstack.translate(bbox.getPosition());
		glUniformMatrix4fv(locationMV, 1, GL_FALSE, MVstack.getCurrentMatrix());
		bbox.render();
		MVstack.pop();
		MVstack.pop(); //Camera transforms >--

		glfwSwapBuffers(currentWindow);
		glfwPollEvents();
	}

	return 0;
}

void inputHandler(GLFWwindow* _window, double _dT)
{
	if (glfwGetKey(_window, GLFW_KEY_ESCAPE)) {
		glfwSetWindowShouldClose(_window, GL_TRUE);
	}

}

void GLcalls()
{
	glClearColor(0.01f, 0.01f, 0.01f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	//glDisable(GL_TEXTURE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}