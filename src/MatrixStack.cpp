#include "MatrixStack.h"
#include <math.h>
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
#include <glm/gtc/type_ptr.hpp> //glm::make:mat4
#include <iostream>

/* Constructor: Create one matrix and set it to identity */
MatrixStack::MatrixStack() {
    currentMatrix = new Matrix;
    init();
};

/* Destructor: Delete all elements from the stack */
MatrixStack::~MatrixStack() {
    flush();
};

// Get the pointer to the topmost (current) matrix array.
float* MatrixStack::getCurrentMatrix() {
    if(currentMatrix)
        return &currentMatrix->m[0][0];
    else
        return nullptr;
};

glm::mat4* MatrixStack::getCurrentMatrixM() {
	if (currentMatrix)
		return &currentMatrix->m;
	else
		return nullptr;
};

// Set the topmost (current) matrix to the identity matrix
// and make it the only element on the stack.
// This should be used with caution.
void MatrixStack::init() {
    currentMatrix->previous = NULL;
	currentMatrix->m = glm::mat4(1);
	
};

// Multiply the topmost (current) matrix with a rotation around X
void MatrixStack::rotX(float _a){
    glm::mat4 rx;
	rx[0][0] = 1.0f; rx[0][1] = 0.0f; rx[0][2] = 0.0f; rx[0][3] = 0.0f;
	rx[1][0] = 0.0f; rx[1][1] = cosf(_a); rx[1][2] = sinf(_a); rx[1][3] = 0.0f;
	rx[2][0] = 0.0f; rx[2][1] = -sinf(_a); rx[2][2] = cosf(_a); rx[2][3] = 0.0f;
	rx[3][0] = 0.0f; rx[3][1] = 0.0f; rx[3][2] = 0.0f; rx[3][3] = 1.0f;
    matrixMult(&currentMatrix->m, &rx, &currentMatrix->m);
};

// Multiply the topmost (current) matrix with a rotation around Y
void MatrixStack::rotY(float _a){
	glm::mat4 ry;
	ry[0][0] = cosf(_a); ry[0][1] = 0.0f; ry[0][2] = -sinf(_a); ry[0][3] = 0.0f;
	ry[1][0]= 0.0f; ry[1][1]=1.0f; ry[1][2]=0.0f; ry[1][3]= 0.0f;
	ry[2][0] = sinf(_a); ry[2][1] = 0.0f; ry[2][2] = cosf(_a); ry[2][3] = 0.0f;
	ry[3][0]=0.0f; ry[3][1]=0.0f; ry[3][2]=0.0f; ry[3][3]=1.0f;
	matrixMult(&currentMatrix->m, &ry, &currentMatrix->m);
};

// Multiply the topmost (current) matrix with a rotation around Z
void MatrixStack::rotZ(float _a) {
    glm::mat4 rz;
	rz[0][0] = cosf(_a); rz[0][1] = sinf(_a); rz[0][2] = 0.0f; rz[0][3] = 0.0f;
	rz[1][0] = -sinf(_a); rz[1][1] = cosf(_a); rz[1][2] = 0.0f; rz[1][3] = 0.0f;
	rz[2][0]= 0.0f; rz[2][1]=0.0f; rz[2][2]=1.0f; rz[2][3]=0.0f;
	rz[3][0]=0.0f; rz[3][1]=0.0f; rz[3][2]=0.0f; rz[3][3]=1.0f;
    matrixMult(&currentMatrix->m, &rz, &currentMatrix->m);
};

void MatrixStack::rotAxis(glm::vec3* _Axis, float _a) {

	glm::mat4 ra;
	glm::rotate(ra, _a, *_Axis);

	matrixMult(&currentMatrix->m, &ra, &currentMatrix->m);
}

// Multiply the topmost (current) matrix with a uniform scaling
void MatrixStack::scale(float _s) {
    glm::mat4 sm;
	sm = glm::scale(sm, glm::vec3(_s, _s, _s));
    matrixMult(&currentMatrix->m, &sm, &currentMatrix->m);
};

// Multiply the topmost (current) matrix with a translation
void MatrixStack::translate(glm::vec3* _Position) {
    glm::mat4 tm;
	//glm::vec3 test = *_Position;
	tm = glm::translate(tm, *_Position);
    matrixMult(&currentMatrix->m, &tm, &currentMatrix->m);
};

void MatrixStack::multiply(glm::mat4* _M) {

	matrixMult(&currentMatrix->m, _M, &currentMatrix->m);
};

// Add a new level on the stack, by making a copy of the current matrix and
// adding it to the top of the linked list, making it the new current matrix.
void MatrixStack::push() {
    Matrix *new_mtx;
    int i;

    new_mtx = new Matrix;

    
    new_mtx->m = currentMatrix->m;
    
    new_mtx->previous = currentMatrix;
    currentMatrix = new_mtx;
};

// Remove one element from the stack, exposing the element below.
void MatrixStack::pop() {
    Matrix *old_mtx;

    old_mtx = currentMatrix;

    if(old_mtx == NULL) {
        fprintf(stderr, "Stack underflow: Pop attemted on empty stack.\n");
        return;
    }
    currentMatrix = currentMatrix->previous;
    delete old_mtx;
};

// Remove all elements from the stack
void MatrixStack::flush() {
    while(currentMatrix != NULL) {
        pop();
    }
}

// Count the number of elements on the stack
int MatrixStack::depth() {
    Matrix *tmp_mtx;
    int i=0;

    tmp_mtx = currentMatrix;
    while(tmp_mtx != NULL) {
        i++;
        tmp_mtx = tmp_mtx->previous;
    }
    return i;
}

// Print the entire contents of a matrix stack (for debugging purposes)
void MatrixStack::print() {
    Matrix *tmp_mtx;

    tmp_mtx = currentMatrix;
    printf("Stack contents:\n");
    while(tmp_mtx != NULL) {
        matrixPrint(&tmp_mtx->m);
        tmp_mtx = tmp_mtx->previous;
        printf("\n");
    }
}


// private:
// Print the contents of a matrix M to stdout.
void MatrixStack::matrixPrint(glm::mat4* M) {
	std::cout << (*M)[0][0] << " " << (*M)[0][1] << " " << (*M)[0][2] << " " << (*M)[0][3] << " "
		<< (*M)[1][0] << " " << (*M)[1][1] << " " << (*M)[1][2] << " " << (*M)[1][3] << " "
		<< (*M)[2][0] << " " << (*M)[2][1] << " " << (*M)[2][2] << " " << (*M)[2][3] << " "
		<< (*M)[3][0] << " " << (*M)[3][1] << " " << (*M)[3][2] << " " << (*M)[3][3] << std::endl;
}

// private:
// Multiply two matrices.
void MatrixStack::matrixMult(glm::mat4* M1, glm::mat4* M2, glm::mat4* Mout) {
    // Compute result in a local variable to avoid conflicts
    // with overwriting if Mout is the same variable as either
    // M1 or M2.
	*Mout = (*M1) * (*M2);

}
