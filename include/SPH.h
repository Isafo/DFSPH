#pragma once

#include "glfwContext.h"
#include "Shader.h"
#include "MatrixStack.h"

class SPH
{
public:
	SPH();
	~SPH();
	void draw_linebox();

	//before loop
	void find_neighborhoods();
	void calculate_densities();
	void calculate_factors();
	
	//in loop
	void non_pressure_forces();
	void calculate_time_step();
	void predict_velocities();
	void correct_density_error();
	void correct_strain_rate_error();
	void update_positions();
	void update_neighbors();
	//
	//void find_neighborhoods();
	//
	void correct_divergence_error();
	void update_velocities();

private:
	struct Float3	
	{
		float* x;
		float* y;
		float* z;
	};

	struct Particles
	{
		Float3 pos;
		Float3 vel;

		float* rad;
		float* p;
		float* dens;
		float* mass;
		float* alpha;
	};
	Particles particles;


};

