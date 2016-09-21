#include "Camera.h"
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
#include <glm/gtc/type_ptr.hpp> //glm::make:mat4


Camera::Camera()
{
	transform = glm::mat4(1.0f);
	perspective = glm::perspective(45.0f, 16.0f / 9.0f, 0.1f, 100.f);

	direction = glm::vec3(0.0f, 0.0f, -1.0f);
	upDirection = glm::vec3(0.0f, 1.0f, 0.0f);
	rightDirection = glm::vec3(1.0f, 0.0f, 0.0f);

	position = glm::vec3(0.0f, 0.0f, 0.0f);

	pitch = 0.0f;
	yaw = -3.14f;
	//yaw = 0.0f;
}


Camera::~Camera()
{

}

void Camera::getPosition(glm::vec3& _Position)
{
	//_Position[0] = transform[3][0];
	//_Position[1] = transform[3][1];
	//_Position[2] = transform[3][2];
	_Position = position;
}

float* Camera::getPositionF()
{
	//return &transform[3][0];
	return &position[0];
}

glm::vec3* Camera::getDirection()
{
	return &direction;
}

glm::vec3* Camera::getUpDirection()
{
	return &upDirection;
}

float* Camera::getTransformF()
{
	return &transform[0][0];
}

glm::mat4* Camera::getTransformM()
{
	return &transform;
}

float* Camera::getPerspective()
{
	return &perspective[0][0];
}

void Camera::translate(glm::vec3* _Translation)
{
	//transform[3][0] += (*_Translation)[0];
	//transform[3][1] += (*_Translation)[1];
	//transform[3][2] += (*_Translation)[2];
	position += *_Translation;
}

void Camera::setPosition(glm::vec3* _Position)
{
	/*transform[3][0] = (*_Position)[0];
	transform[3][1] = (*_Position)[1];
	transform[3][2] = (*_Position)[2];*/
	position = *_Position;
}

void Camera::setTransform(glm::mat4* _Transform)
{
	transform = *_Transform;
}

void Camera::setPerspective(glm::mat4* _Perspective)
{
	perspective = *_Perspective;
}

void Camera::update()
{
	direction = glm::vec3(cos(pitch)*sin(yaw),
						  sin(pitch),
						  cos(pitch)*cos(yaw));

	rightDirection = glm::vec3(sin(yaw - 3.14f / 2.0f),
							   0,
							   cos(yaw - 3.14f / 2.0f));
	upDirection = glm::cross(rightDirection, direction);

	
	glm::vec3 test = glm::vec3(transform[3]);
	transform = glm::lookAt(position, position + direction, upDirection);
}

void Camera::fpsCamera(GLFWwindow* _window, double _dT)
{
	double X, Y, dX, dY;
	glfwGetCursorPos(_window, &X, &Y);

	yaw -= (X - 960.0) / 1920.0;
	pitch -= (Y - 540.0) / 1080.0;
	pitch = std::fmax(std::fmin(pitch, 1.57079f), -1.57079f);

	this->update();

	glm::vec3 translation;
	float movementSpeed = 0.0f;

	if (glfwGetKey(_window, GLFW_KEY_LEFT_SHIFT)){
		movementSpeed = 10.0f;
	}
	else{
		movementSpeed = 1.0f;
	}
	if (glfwGetKey(_window, GLFW_KEY_W)) {
		translation = direction*movementSpeed*(float)_dT;
		this->translate(&translation);
	}
	if (glfwGetKey(_window, GLFW_KEY_S)){
		translation = direction*-movementSpeed*(float)_dT;
		this->translate(&translation);
	}
	if (glfwGetKey(_window, GLFW_KEY_A)) {
		translation = rightDirection*-movementSpeed*(float)_dT;
		this->translate(&translation);
	}
	if (glfwGetKey(_window, GLFW_KEY_D)) {
		translation = rightDirection*movementSpeed*(float)_dT;
		this->translate(&translation);
	}

	glfwSetCursorPos(_window, 960, 540);

}