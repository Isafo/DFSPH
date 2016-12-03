#include "GL/glew.h"
#include <imgui\imgui.h>
#include "imgui\imgui_impl_glfw.h"
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

void start_new_simulation(SPH* sph, int n_particles, int x, int y, int z);

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
	GLint locationColor = glGetUniformLocation(sceneLight.programID, "Color");

	MatrixStack MVstack; MVstack.init();

	//BoundingBox bbox(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
	BoundingBox bbox(0.f, 0.f, 0.f, 1.f, 1.f, 1.f);

	// ugly hax with sending in null
	SPH* current_simulation = new SPH(4000, -0.15f, 0.f, 0.f);
	Sphere sphere(0.0f, 0.0f, 0.0f, current_simulation->get_particle_radius());

	// for testing.
	//delete current_simulation;
	//current_simulation = nullptr;

	Camera mCamera;
	mCamera.setPosition(&glm::vec3(0.f, 0.f, 1.0f));
	mCamera.update();

	bool fpsResetBool = false;

	double lastTime = glfwGetTime() - 0.001;
	double dT = 0.0;

	// GUI variables
	int dParticle = 0;
	int n_particles = current_simulation->get_nr_of_particles();

	while (!glfwWindowShouldClose(currentWindow))
	{
		// Loop for each frame...

		glfwPollEvents();
		if (dT > 1.0 / 30.0) {
			if (current_simulation) {
				current_simulation->update(dT);
			}

			lastTime = glfwGetTime();
		}
		else {
			dT = glfwGetTime() - lastTime;
		}

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

		MVstack.push();//Camera transforms --<
		glUniformMatrix4fv(locationP, 1, GL_FALSE, mCamera.getPerspective());
		MVstack.multiply(mCamera.getTransformM());

		// if there is a simulation running
		glm::vec3 particlePos;
		Float3* particle_pos;
		if (current_simulation != nullptr) {

			particle_pos = current_simulation->get_particle_positions();
			for (int i = 0; i < current_simulation->get_nr_of_particles(); ++i) {
				MVstack.push();
				particlePos = glm::vec3(
					particle_pos->x[i],
					particle_pos->y[i],
					particle_pos->z[i]
					);

				if (i == dParticle)
				{
					float color[] = { 1.0f, 0.3f, 0.0f };
					glUniform3fv(locationColor, 1, &color[0]);
				}
				else
				{
					float color[] = { 0.0f, 0.3f, 1.0f };
					glUniform3fv(locationColor, 1, &color[0]);
				}
				MVstack.translate(&particlePos);
				//MVstack.translate(particle.getPosition());
				glUniformMatrix4fv(locationMV, 1, GL_FALSE, MVstack.getCurrentMatrix());
				sphere.render();
				MVstack.pop();
			}
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
			ImGui::Text("Simulation properties:");
			ImGui::SliderInt("Number of particles: ", &n_particles, 0, 10000);

			if (ImGui::Button("Start")) {
				start_new_simulation(current_simulation, n_particles, 0, 0, 0);
			}

			ImGui::Separator();
			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
			ImGui::Text("Simulation average %.3f ms/frame", dT);

			ImGui::SliderInt("particle: ", &dParticle, 0, current_simulation->get_nr_of_particles() - 1);

			if (current_simulation) {
				Float3s pos = current_simulation->get_pos_i(dParticle);
				Float3s vel = current_simulation->get_vel_i(dParticle);
				Float3s pred_vel = current_simulation->get_predvel_i(dParticle);
				Float3s F_adv = current_simulation->get_F_adv_i(dParticle);
				float dens = current_simulation->get_dens_i(dParticle);

				ImGui::Text("pos: %.4f %.4f %.4f", pos.x, pos.y, pos.z);
				ImGui::Text("vel: %.4f %.4f %.4f", vel.x, vel.y, vel.z);
				ImGui::Text("pred. vel: %.4f %.4f %.4f", pred_vel.x, pred_vel.y, pred_vel.z);
				ImGui::Text("F_adv: %.4f %.4f %.4f", F_adv.x, F_adv.y, F_adv.z);
				ImGui::Text("dens: %.4f", dens);
			}
		}

		// Rendering imgui
		int display_w, display_h;
		glfwGetFramebufferSize(currentWindow, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		ImGui::Render();
		glfwSwapBuffers(currentWindow);

		// check if user respawns particle system
		if (glfwGetKey(currentWindow, GLFW_KEY_R)) {
			start_new_simulation(current_simulation, n_particles, -0.15f, 0.0f, 0.0f);
		}

	}

	ImGui_ImplGlfw_Shutdown();

	if (current_simulation)
		delete current_simulation;

	return 0;
}


void inputHandler(GLFWwindow* _window, double _dT)
{
	if (glfwGetKey(_window, GLFW_KEY_ESCAPE)) {
		glfwSetWindowShouldClose(_window, GL_TRUE);
	}
}


void start_new_simulation(SPH* sph, int n_particles, int x = 0, int y = 0, int z = 0)
{
	delete sph;

	sph = new SPH(n_particles, x, y, z);
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