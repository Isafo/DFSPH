#pragma once
#include "glm/glm.hpp"
#include "iostream"
#define D_NR_OF_PARTICLES 100
#define D_MAX_NR_OF_NEIGHBORS 100

struct Float3
{
	float* x;
	float* y;
	float* z;
};
struct Float3s
{
	float x;
	float y;
	float z;
};

struct Neighbor_Data
{
	int neighbor[D_MAX_NR_OF_NEIGHBORS];
	unsigned int n;
};


class SPH
{
public:
	SPH();

	~SPH();
	//render
	void update(float dT);

	//before loop
	void find_neighborhoods();
	void calculate_densities();
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
	void correct_density_error(float* alpha,float dT, float* g_values, Float3s* f_tot, Float3s* k_v_i);
	void correct_strain_rate_error();
	void update_positions(float dT);
	void correct_divergence_error(float* alpha, Float3* k_v_i);
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
	const float C_NEIGHBOR_RAD = 0.3f;
};

inline void calculate_factors(float* mass, Float3* pos, float* dens, float* g_value, float nr_particles, Neighbor_Data* neighbor_data, float* alpha);
inline void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data, const float NEIGHBOR_RAD);
inline void calculate_kvi(float* alpha, Float3* vel, Float3* pos, float* mass, int nr_particles, float delta_t, Float3s *k_v_i, Neighbor_Data* neighbor_data, float* g_value);
inline void update_function_g(Float3* pos, Neighbor_Data* neighbor_data, float* g, const float NEIGHBOR_RADIUS);
inline void calculate_pressure_force(Float3s* f_tot, Float3s* k_v_i, Float3* pos, float* mass, float* g_val, Neighbor_Data* neighbor_data, float* dens);
inline void calculate_predicted_pressure(Float3s * predicted_pressure, Float3s* f_p, float* mass, float_t*dens, float* g_val, float delta_t, Neighbor_Data* n_data, Float3 * pos, const float rest_Dens);