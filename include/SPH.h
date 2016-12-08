#pragma once

#define D_MAX_NR_OF_NEIGHBORS 100

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

	SPH(int n_particles);

	~SPH();

	void init();

	void reset();

	// performs the simulation steps and updates the particles
	void update(float windX, float windY, float windZ);

	// initializes the particles in a given grid formation
	void init_positions(float x = 0.0f, float y = 0.0f, float z = 0.0f, int rows = 3, int cols = 3) const;

	int get_nr_of_particles() const { return get_n_particles(); }

	float get_particle_radius() const { return m_rad; }

	Float3* get_particle_positions() { return &m_particles.pos; }

	float get_dens_i(int i) const { return m_particles.dens[i]; }
	float get_timestep() const { return m_delta_t; }
	int	get_n_particles() const { return current_n_particles; }

	void set_timestep(float timestep) { m_delta_t = timestep; }
	void set_n_particles(int n_particles) { current_n_particles = n_particles; }
	void set_max_dens_iter(int iter) { Viter_max = iter; }
	void set_max_div_iter(int iter) { iter_max = iter; }
	void set_divergence_error(float error) { divergence_error = error; }
	void set_density_error(float error) { density_error = error; }

	void set_wind(float w_x, float w_y, float w_z)
	{
		m_wind.x = w_x;
		m_wind.y = w_y;
		m_wind.z = w_z;
	}
	void set_gravity(float gravity) { m_gravity = gravity; }

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

	void setStaticSphere(float x, float y, float z, float radius)
	{
		sc.center.x = x;
		sc.center.y = y;
		sc.center.z = z;
		sc.radius_2 = radius*radius;
	}

private:
	sphereConstaint sc;

	// Finds the neighbors of a particle within the given radius D_NEIGBBOR_RAD
	void find_neighborhoods() const;

	// Gravity and wind
	void non_pressure_forces() const;

	void calculate_time_step();

	void predict_velocities(float windX = 0.0f, float windY = 0.0f, float windZ = 0.0f);

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
	const int C_N_PARTICLES;

	float* m_alpha;
	float* m_dens_derive;
	float* m_pred_dens;
	float* m_scalar_values;
	float* m_kernel_values;

	float m_dens_derive_avg;
	float m_pred_dens_avg;

	int current_n_particles;
	int Viter_max{ 100 };
	int iter_max{ 100 };

	float divergence_error{ 0.10f };
	float density_error{ 0.01f };
	float time_factor{ 0.5f };

	// effects
	Float3s m_wind;
	float m_gravity;

};

inline void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data, const int N_PARTICLES);

// updates the scalar values g(q) for all particles
inline void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values, const int N_PARTICLES);