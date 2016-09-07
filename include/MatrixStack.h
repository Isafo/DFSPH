#pragma once
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GL/glew.h"
#include <glm/mat4x4.hpp> // glm::mat4

// A linked list implementing a matrix stack

#ifndef MATRIXSTACK_H
#define MATRIXSTACK_H
struct Matrix {
  glm::mat4 m;
  Matrix *previous;
};
#endif
class MatrixStack {

public:
    Matrix *currentMatrix;

    /* Constructor: Create one matrix and set it to identity */
    MatrixStack();

    /* Destructor: Delete all elements from the stack */
    ~MatrixStack();

    // Get the pointer to the topmost (current) matrix array
    float* getCurrentMatrix();
	glm::mat4* getCurrentMatrixM();

    // Set the topmost (current) matrix to the identity matrix
    void init();

    // Multiply the topmost (current) matrix with a rotation around X
    void rotX(float _a);

    // Multiply the topmost (current) matrix with a rotation around Y
    void rotY(float _a);

    // Multiply the topmost (current) matrix with a rotation around Z
    void rotZ(float _a);

	void rotAxis(glm::vec3* _Axis, float _a);

    // Multiply the topmost (current) matrix with a uniform scaling
    void scale(float _s);

    // Multiply the topmost (current) matrix with a translation
    void translate(glm::vec3* _Position);

	void multiply(glm::mat4* _M);
	//void multiply(float **M);

    // Add a new level on the stack, making a copy of the topmost matrix
    void push();

    // Remove one element from the stack, exposing the element below
    void pop();

    // Remove all elements except the last one from the stack
    void flush();

    // Count the number of elements on the stack
    int depth();

    // Print the entire contents of the stack (for debugging purposes)
    void print();

private:
    // These functions should not be used outside of the class.
	void matrixMult(glm::mat4* M1, glm::mat4* M2, glm::mat4* Mout);
	void matrixPrint(glm::mat4* M);
};
