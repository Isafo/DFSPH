#pragma once
#include "glm/glm.hpp"

#define D_NR_OF_PARTICLES 500
#define D_MAX_NR_OF_NEIGHBORS 500

// A struct containing three arrays (SoA)
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

// Containing information about the index of each neighbor to a particle 
// and the number of neighbors the particle has
struct Neighbor_Data
{
	int neighbor[D_MAX_NR_OF_NEIGHBORS];
	unsigned int n;
};

class SPH
{
public:
	SPH(glm::vec3* start_pos);
	~SPH();
	// performs the simulation steps and updates the particles
	void update(float dT);
	
	// initializes the particles in a given grid formation
	void init_positions(glm::vec3* start_pos, int rows = 3, int cols = 3) const;

	static unsigned int get_nr_of_particles() { return D_NR_OF_PARTICLES; }
	float get_particle_radius() const { return m_particles.rad; }
	Float3* get_particle_positions() { return &m_particles.pos; }

private:

	// Finds the neighbors of a particle within the given radius D_NEIGBBOR_RAD
	void find_neighborhoods() const;

	void pressure_forces() const;
	// Calculates the non-pressure forces: Gravity, surface-tension and vicosity
	void non_pressure_forces() const;
	// Calculates a stable time-step
	void calculate_time_step(float dT);
	// Calculates an unstable predicted velocity
	void predict_velocities();
	void correct_density_error(float* alpha,float dT, float* g_values, Float3s* f_tot, float* k_v_i);
	void correct_strain_rate_error();
	void update_positions(float dT) const;
	void correct_divergence_error(float* k_v_i, float* scalar_values, float* alpha);

	void update_velocities();

	struct Particles
	{
		Float3 pos;
		Float3 vel;
		Float3 F_adv;

		float* p;
		float* dens;
		float mass;
		float rad;
	};

	float m_delta_t;
	Particles m_particles;
	Neighbor_Data *m_neighbor_data;

	const float C_REST_DENS{ 2.861f };
};

// calculates the density and the alpha particle factors
inline void update_density_and_factors(float mass, Float3* pos, float* dens, float* scalar_values,
										Neighbor_Data* neighbor_data, float* alpha, float* kernel_values);

inline void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data);

inline void calculate_pressure_force(Float3s* f_tot, float* k_v_i, Float3* pos, float mass, float* scalar_values, Neighbor_Data* neighbor_data, float* dens);
inline void calculate_predicted_pressure(Float3s* predicted_pressure, Float3s* f_p, float mass, float_t*dens, float* scalar_values, float delta_t, Neighbor_Data* n_data, Float3 * pos, const float rest_dens);

// calculates the k^v_i variable for all particles
inline float calculate_kv(float* alpha, Float3* vel, Float3* pos, float* dens, float delta_t, float *k_v_i, Neighbor_Data* neighbor_data, float* scalar_values);
// updates the scalar values g(q) for all particles
inline void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values);