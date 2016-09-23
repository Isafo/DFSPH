#pragma once

struct Float3
{
	float* x;
	float* y;
	float* z;
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
	void calculate_factors();

	

	unsigned int get_nr_of_particles() const { return m_nr_of_particles; }
	Float3* get_particle_positions() { return &m_particles.pos; }
	
private:
	
	//in loop
	void non_pressure_forces();
	void calculate_time_step();
	void predict_velocities();
	void correct_density_error();
	void correct_strain_rate_error();
	void update_positions(float dT);
	void update_neighbors();
	//
	//void find_neighborhoods();
	//
	void correct_divergence_error();
	void update_velocities(float dT);

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
	Particles m_particles;
	unsigned int m_nr_of_particles;

	const float REST_DENS = 0.1f;
};

