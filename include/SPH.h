#pragma once
#include "glm\glm.hpp"
#include "glm\glm.hpp"

#define D_NR_OF_PARTICLES 100
#define D_MAX_NR_OF_NEIGHBORS 100

struct Float3
{
	float* x;
	float* y;
	float* z;
};

struct Neighbor_Data
{
	int neighbor[D_MAX_NR_OF_NEIGHBORS];
	float g_value[D_MAX_NR_OF_NEIGHBORS];
	unsigned int n;
};


class SPH
{
public:
	SPH();
	SPH(int size);

	~SPH();
	//render
	void update(float dT);

	//before loop
	void find_neighborhoods();
	void calculate_densities();
	//void calculate_factors();
	void init_positions(glm::vec3* start_pos, int rows = 3, int cols = 3);

	unsigned int get_nr_of_particles() const { return m_nr_of_particles; }
	float get_particle_radius() const { return m_particles.rad; }
	Float3* get_particle_positions() { return &m_particles.pos; }

private:	

	//in loop
	void pressure_forces();
	void non_pressure_forces();
	void calculate_time_step();
	void predict_velocities(float dT);
	void correct_density_error(float* alpha);
	void correct_strain_rate_error();
	void update_positions(float dT);
	void update_function_g();
	void correct_divergence_error(float* alpha);
	void update_velocities(float dT);
	//void calculate_kvi();
	
	
	struct Particles
	{
		Float3 pos;
		Float3 vel;
		Float3 F_adv;

		float rad;
		float* p;
		float* dens;
		float* mass;
	};
	float m_delta_t;
	Particles m_particles;
	Neighbor_Data *m_neighbor_data;
	unsigned int m_nr_of_particles;

	const float C_REST_DENS = 0.1f;
};


inline void calculate_factors(float* mass, float* dens, float nr_particles, Neighbor_Data* neighbor_data, float* alpha);
inline void calculate_kvi(float* alpha, Float3* vel, float* mass, int nr_particles, float delta_t, float* k_v_i);