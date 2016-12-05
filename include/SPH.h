#pragma once

#define D_MAX_NR_OF_NEIGHBORS 50

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

	SPH(int n_particles);

	void init();
	void reset();
	// Free the memory
	~SPH();

	// performs the simulation steps and updates the particles
	void update(float dT);

	// initializes the particles in a given grid formation
	void init_positions(int x = 0, int y = 0, int z = 0, int rows = 3, int cols = 3) const;

	int get_nr_of_particles() { return C_N_PARTICLES; }

	float get_particle_radius() const { return m_rad; }

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

	float get_dens_i(int i) const { return m_particles.dens[i]; }
	float get_timestep() const { return m_delta_t; }

	int current_n_particles;
	int Viter_max{ 100 };
	int iter_max{ 100 };

	float divergence_error{ 0.10f };
	float density_error{ 0.01f };
	float time_factor{ 0.5f };


private:

	// Finds the neighbors of a particle within the given radius D_NEIGBBOR_RAD
	void find_neighborhoods() const;

	// Calculates the non-pressure forces: Gravity, surface-tension and vicosity
	void non_pressure_forces() const;

	// Calculates a stable time-step
	void calculate_time_step();

	// Calculates an unstable predicted velocity
	void predict_velocities();

	void correct_density_error();

	void update_positions() const;

	void correct_divergence_error();

	void update_velocities();

	// former inline functions
	void calculate_derived_density_pred_dens(Neighbor_Data* neighbor_data);

	// calculates the density and the alpha particle factors
	void update_density_and_factors(Neighbor_Data* neighbor_data);

	struct Particles
	{
		Float3 pos;
		Float3 vel;
		Float3 pred_vel;
		Float3 F_adv;

		float* p;
		float* dens;
	};

	struct sphereConstaint
	{
		Float3s center;
		Float3s normal;
		float radius_2;
	};

	float m_delta_t;
	float m_mass;
	float m_rad;

	Particles m_particles;
	Neighbor_Data *m_neighbor_data;

	const float C_REST_DENS{ 1000.f };
	const int C_N_PARTICLES;

	float* m_alpha;
	float* m_dens_derive;
	float* m_pred_dens;
	float* m_scalar_values;
	float* m_kernel_values;

	float m_dens_derive_avg;
	float m_pred_dens_avg;
};

inline void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data, const int N_PARTICLES);

// updates the scalar values g(q) for all particles
inline void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values, const int N_PARTICLES);