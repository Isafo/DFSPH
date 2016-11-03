#include "GL/glew.h"
#include <imgui.h>
#include "imgui_impl_glfw.h"
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
	glfwSetInputMode(currentWindow, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
	glfwSetCursorPos(currentWindow, 960, 540);

	// Setup ImGui binding
	ImGui_ImplGlfw_Init(currentWindow, true);

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

	SPH s(&glm::vec3(-0.15, 0.0, 0.0));
	
	Sphere sphere(0.0f, 0.0f, 0.0f, s.get_particle_radius());

	Camera mCamera;
	mCamera.setPosition(&glm::vec3(0.f, 0.f, 1.0f));
	mCamera.update();

	bool fpsResetBool = false;
	
	double lastTime = glfwGetTime() - 0.001f;
	double dT = 0.0;
	while (!glfwWindowShouldClose(currentWindow))
	{
		glfwPollEvents();
		
		dT = glfwGetTime() - lastTime;
		lastTime = glfwGetTime();

		//glfw input handler
		inputHandler(currentWindow, dT);

		if (glfwGetKey(currentWindow, GLFW_KEY_LEFT_CONTROL)) 
		{
			if (!fpsResetBool)
			{
				fpsResetBool = true;
				glfwSetCursorPos(currentWindow, 960, 540);
			}
			
			mCamera.fpsCamera(currentWindow, dT);
		}
		else
		{
			fpsResetBool = false;
		}
		
		GLcalls();

		glUseProgram(sceneLight.programID);

		s.update(dT);

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

		glUseProgram(0);

		ImGui_ImplGlfw_NewFrame();


		{
			static int dParticle = 0;
			ImGui::Text("Hello, world!");
			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

			ImGui::SliderInt("particle: ", &dParticle, 0, 10);

			Float3s pos = s.get_pos_i(dParticle);
			Float3s vel = s.get_vel_i(dParticle);
			Float3s F_adv = s.get_F_adv_i(dParticle);
			float p = s.get_p_i(dParticle);
			float dens = s.get_dens_i(dParticle);

			ImGui::Text("pos: %.1f %.1f %.1f", pos.x, pos.y, pos.z);
			ImGui::Text("vel: %.1f %.1f %.1f", vel.x, vel.y, vel.z);
			ImGui::Text("F_adv: %.1f %.1f %.1f", F_adv.x, F_adv.y, F_adv.z);
			ImGui::Text("p: %.1f", p);
			ImGui::Text("dens: %.1f", dens);

		}

		// Rendering imgui
		int display_w, display_h;
		glfwGetFramebufferSize(currentWindow, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		ImGui::Render();
		glfwSwapBuffers(currentWindow);
	}

	ImGui_ImplGlfw_Shutdown();

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