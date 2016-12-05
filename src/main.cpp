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
	GLint locationColor = glGetUniformLocation(sceneLight.programID, "Color");

	MatrixStack MVstack; MVstack.init();

	BoundingBox bbox(0.f, 0.f, 0.f, 1.f, 1.f, 1.f);

	SPH sph(5000);
	sph.init();

	Sphere sphere(0.0f, 0.0f, 0.0f, sph.get_particle_radius());
	Sphere static_sphere(0.0f, -0.5f, 0.0f, 0.25f);
	
	Camera mCamera;
	mCamera.setPosition(&glm::vec3(0.f, 0.f, 1.0f));
	mCamera.update();

	// Time related variables
	bool is_running = false;
	bool is_paused = false;

	bool fpsResetBool = false;

	double lastTime = glfwGetTime() - 0.001;
	double dT = 0.0;

	// GUI variables
	int dParticle = 0;
	int n_particles = 4000;

	int Miter_v = 200, Miter = 200;

	float dens_error = 0.1f, div_error = 0.01f;
	float max_error = 0.2f , min_error = 0.01f;
	float time_factor = 0.4f;
	static bool addImplicitSphere = false;

	float original_sphereRad = 1.0f;
	float sphereX = 0.f, sphereY = -0.25f, sphereZ = 0.0f, sphereRad = 0.16f;


	while (!glfwWindowShouldClose(currentWindow))
	{
		glfwPollEvents();

		ImGui_ImplGlfw_NewFrame();
		{

		
			ImGui::Text("Simulation properties");
			ImGui::InputInt("Number of particles", &n_particles, 100, 1);
			n_particles = glm::clamp(n_particles, 100, 8000);

			ImGui::InputInt("Max iter. (Divergence solv.)", &Miter_v, 10, 1);
			ImGui::InputInt("Max iter. (Density solv.)", &Miter, 10, 1);
			Miter_v = glm::clamp(Miter_v, 10, 200);
			Miter = glm::clamp(Miter, 10, 200);
			
			ImGui::InputFloat("Error % (Divergence)", &div_error, 0.01f, 1, 2);
			ImGui::InputFloat("Error % (Density)", &dens_error, 0.01f, 1, 2);
			div_error = glm::clamp(div_error, min_error, max_error);
			dens_error = glm::clamp(dens_error, min_error, max_error);

			ImGui::InputFloat("CFL - Factor", &time_factor, 0.01f, 1, 2);
			time_factor = glm::clamp(time_factor, 0.4f, 0.6f);
			//ImGui::Spacing();

			ImGui::Checkbox("checkbox", &addImplicitSphere);
			if (addImplicitSphere) {
				ImGui::InputFloat("input x", &sphereX);
				ImGui::InputFloat("input y", &sphereY);
				ImGui::InputFloat("input z", &sphereZ);
				ImGui::InputFloat("Radius", &sphereRad);
				sph.setStaticSphere(sphereX, sphereY, sphereZ, sphereRad);
				static_sphere.setRadius(sphereRad);
				static_sphere.setPosition(glm::vec3(sphereX, sphereY, sphereZ));
				
			}

			if (ImGui::Button("Start")) {
				sph.set_n_particles(n_particles);
				sph.set_max_dens_iter(Miter_v);
				sph.set_max_div_iter(Miter);
				sph.set_divergence_error(div_error);
				sph.set_density_error(dens_error);
				sph.set_timestep(time_factor);
				
				sph.reset();
				sph.init_positions(0, 0, 0, 20, 20);
				is_running = true;
			}ImGui::SameLine();

			if (ImGui::Button("Pause")) {
				is_paused = !is_paused;
			}

			if (is_running) {
				ImGui::Separator();
				ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
				ImGui::Text("Simulation average %.3f ms/frame", dT);

				ImGui::SliderInt("particle: ", &dParticle, 0, n_particles - 1);

				Float3s pos = sph.get_pos_i(dParticle);
				Float3s vel = sph.get_vel_i(dParticle);
				Float3s pred_vel = sph.get_predvel_i(dParticle);
				//Float3s F_adv = sph.get_F_adv_i(dParticle);
				float dens = sph.get_dens_i(dParticle);

				ImGui::Text("pos: %.4f %.4f %.4f", pos.x, pos.y, pos.z);
				ImGui::Text("vel: %.4f %.4f %.4f", vel.x, vel.y, vel.z);
				ImGui::Text("pred. vel: %.4f %.4f %.4f", pred_vel.x, pred_vel.y, pred_vel.z);
				//ImGui::Text("F_adv: %.4f %.4f %.4f", F_adv.x, F_adv.y, F_adv.z);
				ImGui::Text("dens: %.4f", dens);
				ImGui::Text("Time Step: %.4f", sph.get_timestep());
			}
		
		}

		if (dT > 1.0 / 30.0) {
			if (is_running && !is_paused) {
				sph.update();
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
		if (is_running) {

			particle_pos = sph.get_particle_positions();
			for (int i = 0; i < sph.get_n_particles(); ++i) {
				MVstack.push();
				particlePos = glm::vec3(
					particle_pos->x[i],
					particle_pos->y[i],
					particle_pos->z[i]
					);

				if (i == dParticle)
				{
					float color[] = { 1.0f, 0.3f, 0.0f, 1.0f };
					glUniform4fv(locationColor, 1, &color[0]);
				}
				else
				{
					float color[] = { 0.0f, 0.3f, 1.0f, 1.0f };
					glUniform4fv(locationColor, 1, &color[0]);
				}
				MVstack.translate(&particlePos);
				//MVstack.translate(particle.getPosition());
				glUniformMatrix4fv(locationMV, 1, GL_FALSE, MVstack.getCurrentMatrix());
				sphere.render();
				MVstack.pop();
			}
		}

		if (addImplicitSphere) 
		{
			MVstack.push();
			float color[] = { 0.0f, 1.0f, 0.5f };
			glUniform3fv(locationColor, 1, &color[0]);
			MVstack.translate(static_sphere.getPosition());
			MVstack.scale(static_sphere.getRadius()/(original_sphereRad+0.01f));
			glUniformMatrix4fv(locationMV, 1, GL_FALSE, MVstack.getCurrentMatrix());
			static_sphere.render();
			MVstack.pop();
		}

		MVstack.push();
		MVstack.translate(bbox.getPosition());
		glUniformMatrix4fv(locationMV, 1, GL_FALSE, MVstack.getCurrentMatrix());
		bbox.render();
		MVstack.pop();
		MVstack.pop(); //Camera transforms >--

		glUseProgram(0);

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
	glDisable(GL_TEXTURE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}