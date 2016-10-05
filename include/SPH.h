#pragma once
#include "glm/glm.hpp"

#define D_NR_OF_PARTICLES 2
#define D_MAX_NR_OF_NEIGHBORS 100

// A struct containing three arrays (SoA)
struct Float3
{
	float* x;
	float* y;
	float* z;
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
	SPH();

	void update_velocities(float dT) const;

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
	void predict_velocities(float dT);
	// Correct the density error in the predicted velocity
	void correct_density_error(float* alpha);

	void correct_strain_rate_error();
	// Update the particle positions
	void update_positions(float dT) const;
	// correct the divergence error in the predicted velocity
	void correct_divergence_error(float* alpha);
	// Update the velocity with the stable corrected predicted velocity
	void update_velocities(float dT);

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

	const float C_REST_DENS{ 0.1f };
	const float C_NEIGHBOR_RAD{ 0.3f };
};

// calculates the density and the alpha particle factors
inline void update_density_and_factors(float* mass, Float3* pos, float* dens, float* g_value, float nr_particles, 
										Neighbor_Data* neighbor_data, float* alpha, float* kernel_values);

inline void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data, const float NEIGHBOR_RAD);
// calculates the k^v_i variable for all particles
inline void calculate_kvi(float* alpha, Float3* vel, float* mass, int nr_particles, float delta_t, float* k_v_i);
// updates the scalar values g(q) for all particles
inline void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* g, const float NEIGHBOR_RADIUS);
