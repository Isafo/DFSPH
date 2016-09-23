#include "SPH.h"

#define gravity -9.82f

SPH::SPH()
{
}

SPH::SPH(int size) {
	m_nr_of_particles = size;
	m_particles.alpha = new float[size];
	m_particles.dens = new float[size];
	m_particles.mass = new float[size];
	m_particles.p = new float[size];
	m_particles.rad = new float[size];
	m_particles.pos.x = new float[size];
	m_particles.pos.y = new float[size];
	m_particles.pos.z = new float[size];
	m_particles.vel.x = new float[size];
	m_particles.vel.y = new float[size];
	m_particles.vel.z = new float[size];

	for (int i = 0; i < m_nr_of_particles; ++i) {
		m_particles.alpha[i] = 0.f;
		m_particles.dens[i] = 0.f;
		m_particles.mass[i] = 0.f;
		m_particles.p[i] = 0.f;
		m_particles.rad[i] = 0.f;
		m_particles.pos.x[i] = 0.f;
		m_particles.pos.y[i] = 0.f;
		m_particles.pos.z[i] = 0.f;
		m_particles.vel.x[i] = 0.f;
		m_particles.vel.y[i] = 0.f;
		m_particles.vel.z[i] = 0.f;
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
	for(unsigned int i = 0; i < m_nr_of_particles; ++i)
	{
		m_particles.pos.x[i] += m_particles.vel.x[i] * dT;
		m_particles.pos.y[i] += m_particles.vel.y[i] * dT;
		m_particles.pos.z[i] += m_particles.vel.z[i] * dT;
	}
}
void SPH::update_neighbors() {}
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
	delete[] m_particles.rad;
}
