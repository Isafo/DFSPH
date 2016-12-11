#pragma once

#define D_MAX_NR_OF_PARTICLES 20000
#define D_MAX_NR_OF_NEIGHBORS 30

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

struct sphereConstaint
{
	Float3s center;
	float radius_2;
};

class SPH
{
public:

	SPH();

	~SPH();

	void reset();

	// performs the simulation steps and updates the particles
	void update();

	// initializes the particles in a given grid formation
	void init_positions(Float3s pos, int rows = 3, int cols = 3) const;

	float get_particle_radius() const { return m_rad; }

	Float3* get_particle_positions() { return &m_particles.pos; }

	float get_dens_i(int i) const { return m_particles.dens[i]; }
	float get_timestep() const { return m_delta_t; }
	int	get_nr_of_particles() const { return m_simulated_particles; }
	float get_gravity() const { return m_gravity; }

	void set_timestep(float timestep) { m_delta_t = timestep; }
	void set_nr_of_particles(int n_particles) { m_simulated_particles = n_particles; }
	void set_max_dens_iter(int iter) { m_Viter_max = iter; }
	void set_max_div_iter(int iter) { m_iter_max = iter; }
	void set_divergence_error(float error) { m_divergence_error = error; }
	void set_density_error(float error) { m_density_error = error; }
	void set_gravity(float gravity) { m_gravity = gravity; }

	void set_wind(Float3s wind)
	{
		m_wind.x = wind.x;
		m_wind.y = wind.y;
		m_wind.z = wind.z;
	}

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

	Float3* get_vel()
	{
		return &m_particles.vel;
	}

	void setStaticSphere(float pos_x, float pos_y, float pos_z, float radius)
	{
		sc.center.x = pos_x;
		sc.center.y = pos_y;
		sc.center.z = pos_z;
		sc.radius_2 = radius*radius;
	}

private:
	sphereConstaint sc;

	// Finds the neighbors of a particle within the given radius D_NEIGBBOR_RAD
	void find_neighborhoods() const;

	// Gravity and wind
	void non_pressure_forces();

	void calculate_time_step();

	void predict_velocities();

	void correct_density_error();

	void update_positions() const;

	void correct_divergence_error();

	void update_velocities();

	void calculate_derived_density_pred_dens(Neighbor_Data* neighbor_data);

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


	float m_delta_t;
	float m_mass;
	float m_rad;

	Particles m_particles;
	Neighbor_Data *m_neighbor_data;

	const float C_REST_DENS{ 1000.f };

	float* m_alpha;
	float* m_dens_derive;
	float* m_pred_dens;
	float* m_scalar_values;
	float* m_kernel_values;

	float m_dens_derive_avg;
	float m_pred_dens_avg;

	int m_simulated_particles;
	int m_Viter_max{ 100 };
	int m_iter_max{ 100 };

	float m_divergence_error{ 0.10f };
	float m_density_error{ 0.01f };
	float m_time_factor{ 0.5f };

	// effects
	Float3s m_wind{ 0.0f };
	float m_gravity{ 9.82f };
};

inline void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data, const int N_PARTICLES);

// updates the scalar values g(q) for all particles
inline void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values, const int N_PARTICLES);