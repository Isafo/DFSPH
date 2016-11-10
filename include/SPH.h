#pragma once
#include "glm/glm.hpp"

#define D_NR_OF_PARTICLES 400
#define D_MAX_NR_OF_NEIGHBORS 400

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

	// for debug
	Float3s get_pos_i(int i) const
	{
		Float3s i_pos;
		i_pos.x = m_particles.pos.x[i];
		i_pos.y = m_particles.pos.y[i];
		i_pos.z = m_particles.pos.z[i];

		return i_pos;
	}

	Float3s get_vel_i(int i) const
	{
		Float3s i_vel;
		i_vel.x = m_particles.vel.x[i];
		i_vel.y = m_particles.vel.y[i];
		i_vel.z = m_particles.vel.z[i];

		return i_vel;
	}

	Float3s get_predvel_i(int i) const
	{
		Float3s i_vel;
		i_vel.x = m_particles.pred_vel.x[i];
		i_vel.y = m_particles.pred_vel.y[i];
		i_vel.z = m_particles.pred_vel.z[i];

		return i_vel;
	}


	Float3s get_F_adv_i(int i) const
	{
		Float3s i_f;
		i_f.x = m_particles.F_adv.x[i];
		i_f.y = m_particles.F_adv.y[i];
		i_f.z = m_particles.F_adv.z[i];

		return i_f;
	}

	float get_p_i(int i) const { return m_particles.p[i]; }
	float get_dens_i(int i) const { return m_particles.dens[i]; }

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
	void correct_density_error(float* alpha, float dT, float* g_values, Float3s* f_tot, float* k_v_i);
	void correct_strain_rate_error();
	void update_positions(float dT) const;
	void correct_divergence_error(float* k_v_i, float* scalar_values, float* alpha);

	void update_velocities();

	struct Particles
	{
		Float3 pos;
		Float3 vel;
		Float3 pred_vel;
		Float3 F_adv;

		float* p;
		float* dens;
		float mass;
		float rad;
	};

	int iter = 0;
	float m_delta_t;
	Particles m_particles;
	Neighbor_Data *m_neighbor_data;

	// TODO dont forfill the condition (p - p0 = 0) right now
	const float C_REST_DENS{ 95.f };
};

// calculates the density and the alpha particle factors
void update_density_and_factors(float mass, Float3* pos, float* dens, float* scalar_values,
	Neighbor_Data* neighbor_data, float* alpha, float* kernel_values);

void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data);

void calculate_pressure_force(Float3s* f_tot, float* k_v_i, Float3* pos, float mass, float* scalar_values, Neighbor_Data* neighbor_data, float* dens);
void calculate_predicted_pressure(Float3s* predicted_pressure, Float3s* pred_vel, float mass, float_t*dens, float* scalar_values, float delta_t, Neighbor_Data* n_data, Float3 * pos);

// calculates the k^v_i variable for all particles
float calculate_stiffness(float* alpha, Float3* vel, Float3* pred_vel, Float3* pos, float* dens, float delta_t, float *k_v_i, Neighbor_Data* neighbor_data, float* scalar_values,float mass);
// updates the scalar values g(q) for all particles
void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values);