#include "SPH.h"

#define gravity -9.82f
#define neighbor_rad 0.3f;

SPH::SPH()
{
}

SPH::SPH(int size) {
	m_nr_of_particles = size;

	m_particles.alpha = new float[m_nr_of_particles];
	m_particles.dens = new float[m_nr_of_particles];
	m_particles.mass = new float[m_nr_of_particles];
	m_particles.p = new float[m_nr_of_particles];
	m_particles.pos.x = new float[m_nr_of_particles];
	m_particles.pos.y = new float[m_nr_of_particles];
	m_particles.pos.z = new float[m_nr_of_particles];
	m_particles.vel.x = new float[m_nr_of_particles];
	m_particles.vel.y = new float[m_nr_of_particles];
	m_particles.vel.z = new float[m_nr_of_particles];
	m_neighboors = new int*[m_nr_of_particles];
	
	/* The radius should be read from a
	settings class with static values
	instead of being defined here. */
	m_particles.rad = 0.01f;

	for (int i = 0; i < m_nr_of_particles; ++i) {
		m_particles.alpha[i] = 0.f;
		m_particles.dens[i] = 0.f;
		m_particles.mass[i] = 0.f;
		m_particles.p[i] = 0.f;
		m_particles.pos.x[i] = 0.f;
		m_particles.pos.y[i] = 0.f;
		m_particles.pos.z[i] = 0.f;
		m_particles.vel.x[i] = 0.f;
		m_particles.vel.y[i] = 0.f;
		m_particles.vel.z[i] = 0.f;
		m_neighboors[i] = new int[100];
	}
}

//v = v0+a*dt
//s = v*dt
void SPH::update(float dT) {
	update_positions(dT);
	update_velocities(dT);
}

void SPH::find_neighborhoods()
{

}
void SPH::calculate_densities()
{
	// TODO: LATER UPDATE THIS TO ONLY LOOP OVER NEIGHBOURS
	/*for (int i = 0; i < m_nr_of_particles; ++i)
	{*/
	m_particles.dens[0] = m_particles.mass[1];
	m_particles.dens[1] = m_particles.mass[0];
	//}
}

void SPH::calculate_factors()
{

}

void SPH::init_positions(glm::vec3* start_pos, int rows, int cols)
{
	float dist_between = 2.f * m_particles.rad;
	float padding_factor = 1.1f;

	float x, y, z;
	int ind;

	for (int i = 0; i < m_nr_of_particles; ++i)
	{
		z = i / (rows*cols);
		ind = i - z*rows*cols;
		y = ind / rows;
		x = ind % rows;

		m_particles.pos.x[i] = start_pos->x + x * dist_between*padding_factor;
		m_particles.pos.y[i] = start_pos->y + y * dist_between*padding_factor;
		m_particles.pos.z[i] = start_pos->z + z * dist_between*padding_factor;
	}
}

void SPH::non_pressure_forces()
{
	// calculate the particle pressure
	for (int i = 0; i < m_nr_of_particles; ++i)
	{
		m_particles.p[i] = m_particles.dens[i] - REST_DENS;
	}
}

void SPH::calculate_time_step() {}
void SPH::predict_velocities() {}
void SPH::correct_density_error() {}
void SPH::correct_strain_rate_error() {}
void SPH::update_positions(float dT)
{
	for (unsigned int i = 0; i < m_nr_of_particles; ++i)
	{
		m_particles.pos.x[i] += m_particles.vel.x[i] * dT;
		m_particles.pos.y[i] += m_particles.vel.y[i] * dT;
		m_particles.pos.z[i] += m_particles.vel.z[i] * dT;
	}
}
//TODO, update with lookuptable. instead of bruteforcing it. 
void SPH::update_neighbors() 
{
	float vector_i_n[3];
	float dist_i_n_2;
	const float neigborhod_rad_2 = neighbor_rad;
	int count = 0;
	for (int i = 0; i < m_nr_of_particles; i++) 
	{
		for(int n = 0; n < m_nr_of_particles; n++)
		{
			if (i != n )
			{
				vector_i_n[0] = m_particles.pos.x[n] - m_particles.pos.x[i];
				vector_i_n[1] = m_particles.pos.y[n] - m_particles.pos.y[i];
				vector_i_n[2] = m_particles.pos.z[n] - m_particles.pos.z[i];
				dist_i_n_2 = vector_i_n[0] * vector_i_n[0] + vector_i_n[1] * vector_i_n[1] + vector_i_n[0] * vector_i_n[0];
				
				if (dist_i_n_2 < neigborhod_rad_2*neigborhod_rad_2)
				{
					m_neighboors[i][count] = n;
					++count;
				}
			}
			count = 0;
		}
	}
}
void SPH::correct_divergence_error() {}
void SPH::update_velocities(float dT)
{
	for (unsigned int i = 0; i < m_nr_of_particles; ++i)
	{
		//m_particles.vel.x[i] += 0.f * dT;
		m_particles.vel.y[i] += gravity * dT;
		//m_particles.pos.z[i] += 0.f * dT;
	}
}

SPH::~SPH()
{
	delete[] m_particles.vel.x;
	delete[] m_particles.vel.y;
	delete[] m_particles.vel.z;
	delete[] m_particles.pos.x;
	delete[] m_particles.pos.y;
	delete[] m_particles.pos.z;
	delete[] m_particles.alpha;
	delete[] m_particles.dens;
	delete[] m_particles.mass;
	delete[] m_particles.p;
}
