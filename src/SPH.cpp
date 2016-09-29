#include "SPH.h"

#include <math.h>

#define D_GRAVITY -9.82f
#define D_NEIGHBOR_RAD 0.3f;
#define D_PI 3.1415926559f;
#define D_EPSILON 0.000000000000001f;

SPH::SPH() 
{
	//m_particles.alpha = new float[D_NR_OF_PARTICLES];
	m_particles.dens = new float[D_NR_OF_PARTICLES];
	m_particles.mass = new float[D_NR_OF_PARTICLES];
	m_particles.p = new float[D_NR_OF_PARTICLES];
	m_particles.pos.x = new float[D_NR_OF_PARTICLES];
	m_particles.pos.y = new float[D_NR_OF_PARTICLES];
	m_particles.pos.z = new float[D_NR_OF_PARTICLES];
	m_particles.vel.x = new float[D_NR_OF_PARTICLES];
	m_particles.vel.y = new float[D_NR_OF_PARTICLES];
	m_particles.vel.z = new float[D_NR_OF_PARTICLES];
	m_particles.F_adv.x = new float[D_NR_OF_PARTICLES];
	m_particles.F_adv.y = new float[D_NR_OF_PARTICLES];
	m_particles.F_adv.z = new float[D_NR_OF_PARTICLES];
	m_neighbor_data = new Neighbor_Data[D_NR_OF_PARTICLES];

	/* The radius should be read from a
	settings class with static values
	instead of being defined here. */
	m_particles.rad = 0.01f;
	  
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i) {
		//m_particles.alpha[i] = 0.f;
		m_particles.dens[i] = 0.f;
		m_particles.mass[i] = 0.f;
		m_particles.p[i] = 0.f;
		m_particles.pos.x[i] = 0.f;
		m_particles.pos.y[i] = 0.f;
		m_particles.pos.z[i] = 0.f;
		m_particles.vel.x[i] = 0.f;
		m_particles.vel.y[i] = 0.f;
		m_particles.vel.z[i] = 0.f;
	}
}

void SPH::update(float dT) 
{
	float alpha[D_NR_OF_PARTICLES];
	static float kernel_values[D_NR_OF_PARTICLES*D_MAX_NR_OF_NEIGHBORS];

	find_neighborhoods();

	find_neighborhoods();

	update_kernel_values(kernel_values, &m_particles.pos, m_neighbor_data);

	calculate_densities();

	calculate_factors(m_particles.mass, &m_particles.pos, m_particles.dens, D_NR_OF_PARTICLES, m_neighbor_data, alpha);
	
	non_pressure_forces();
	
	update_positions(dT);
	
	update_velocities(dT);
}



void SPH::find_neighborhoods()
{
	float vector_i_n[3];
	float dist_i_n_2;
	const float neigborhod_rad_2 = D_NEIGHBOR_RAD;
	int count = 0;
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		for (int n = 0; n < D_NR_OF_PARTICLES; ++n)
		{
			if (i != n)
			{
				vector_i_n[0] = m_particles.pos.x[n] - m_particles.pos.x[i];
				vector_i_n[1] = m_particles.pos.y[n] - m_particles.pos.y[i];
				vector_i_n[2] = m_particles.pos.z[n] - m_particles.pos.z[i];
				dist_i_n_2 = vector_i_n[0] * vector_i_n[0] + vector_i_n[1] * vector_i_n[1] + vector_i_n[0] * vector_i_n[0];

				if (dist_i_n_2 < neigborhod_rad_2*neigborhod_rad_2)
				{
					m_neighbor_data[i].neighbor[count] = n;
					++count;
				}
			}
		}
		//save nr of neighbor to first position 
		m_neighbor_data[i].n = count;
		count = 0;
	}
}

// TODO: Add kernal function
void SPH::calculate_densities()
{
	int n_neighbors;
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{	
		n_neighbors = m_neighbor_data[i].n;
		for (int n = 0; n < n_neighbors ; ++n)
		{
			m_particles.dens[i] += m_particles.mass[m_neighbor_data[i].neighbor[n]]; // * kernel function
		}
	}
}

// TODO: Add gradient kernal function
inline void calculate_factors(float* mass, Float3* pos, float* dens, float nr_particles, Neighbor_Data* neighbor_data, float* alpha)
{
	int nr_neighbors;
	float abs_sum_denom{ 0 };
	float temp;

	float sum_abs_denom = 0;
	float denom;

	float particle_mass;
	glm::vec3 particle_pos;
	glm::vec3 kernel_gradient;

	unsigned int neighbor_index;

	for (int i = 0; i < nr_particles; i++)
	{
		nr_neighbors = neighbor_data[i].n;
		for (int neighbor = 0; neighbor < nr_neighbors; ++neighbor)
		{
			neighbor_index = neighbor_data[i].neighbor[neighbor];

			particle_mass = mass[neighbor_index];

			particle_pos = glm::vec3(pos->x[neighbor_index],
									pos->y[neighbor_index],
									pos->z[neighbor_index]);

			kernel_gradient = particle_pos * neighbor_data[i].g_value[neighbor_index];

			sum_abs_denom = glm::length(particle_mass*kernel_gradient);
			abs_sum_denom += glm::length(particle_mass*kernel_gradient)*glm::length(particle_mass*kernel_gradient);
		}
		alpha[i] = glm::length(sum_abs_denom)*glm::length(sum_abs_denom) + abs_sum_denom;
	}
}

void SPH::init_positions(glm::vec3* start_pos, int rows, int cols)
{
	float dist_between = 2.f * m_particles.rad;
	float padding_factor = 1.1f;

	float x, y, z;
	int ind;

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
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
	for (int i = 0; i < D_NR_OF_PARTICLES; i++) 
	{
		m_particles.F_adv.x[i] = 0.f;
		m_particles.F_adv.y[i] = D_GRAVITY;
		m_particles.F_adv.z[i] = 0.f;
	}
}
void SPH::pressure_forces()
{
	// calculate the particle pressure
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.p[i] = m_particles.dens[i] - C_REST_DENS;
	}
}
void SPH::calculate_time_step() 
{
	float v_max_2 = 0;
	float x, y, z;
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		x = m_particles.pos.x[i] * m_particles.pos.x[i];
		y = m_particles.pos.y[i] * m_particles.pos.y[i];
		z = m_particles.pos.z[i] * m_particles.pos.z[i];
		if (v_max_2 < (x + y + z))v_max_2 = (x + y + z);
	}
	m_delta_t = 0.4f*(m_particles.rad * 2) / sqrtf(v_max_2) - D_EPSILON;
}

void SPH::predict_velocities(float dT) 
{
	for (unsigned int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.vel.x[i] += m_particles.F_adv.x[i] * dT;
		m_particles.vel.y[i] += m_particles.F_adv.y[i] * dT;
		m_particles.pos.z[i] += m_particles.F_adv.z[i] * dT;
	}
}

//TODO: add kernel function
void SPH::correct_density_error(float* alpha)
{

}

void SPH::correct_strain_rate_error() {}

void SPH::update_positions(float dT)
{
	for (unsigned int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.pos.x[i] += m_particles.vel.x[i] * dT;
		m_particles.pos.y[i] += m_particles.vel.y[i] * dT;
		m_particles.pos.z[i] += m_particles.vel.z[i] * dT;
	}
}

void SPH::update_function_g()
{
	float kernel_derive;
	//Loop through all particles
	for (int particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		//Loop through particle neibors
		for (int neighbor = 0; neighbor < m_neighbor_data[particle].n; ++neighbor)
		{
			auto neighbor_ind = m_neighbor_data[particle].neighbor[neighbor];

			auto dx = m_particles.pos.x[particle] - m_particles.pos.x[neighbor_ind];
			auto dy = m_particles.pos.y[particle] - m_particles.pos.y[neighbor_ind];
			auto dz = m_particles.pos.z[particle] - m_particles.pos.z[neighbor_ind];

			auto dist = std::sqrt(dx*dx + dy*dy + dz*dz);
			auto q = dist / m_particles.rad;

			//Compute the derivitive of the kernel function
			
			if (q >= 0 || q <= 0.5f)
			{
				kernel_derive = (-12.f*q + 18.f*q*q);
			}
			else if (q > 0.5f || q <= 1.0f)
			{
				kernel_derive = -6.0f*(1 - q)*(1 - q);
			}
			else
			{
				kernel_derive = 0;
			}

			m_neighbor_data[particle].g_value[neighbor] = kernel_derive / (m_particles.rad*dist);
		}
	}
}

inline void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data)
{
	unsigned int ind;
	float x, y, z, q;
	float kernel_val;
	for (int particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		float particle_pos_x = pos->x[particle];
		float particle_pos_y = pos->y[particle];
		float particle_pos_z = pos->z[particle];
		for (int neighbor = 0; neighbor < neighbor_data[particle].n; ++neighbor)
		{
			// compute q
			ind = neighbor_data[particle].neighbor[neighbor];
			x = pos->x[ind] - particle_pos_x;
			y = pos->y[ind] - particle_pos_y;
			z = pos->z[ind] - particle_pos_z;

			float len = x*x + y*y + z*z;
			q = sqrt(x*x + y*y + z*z) / D_NEIGHBOR_RAD;

			kernel_val = 8.0f / D_PI;

			if (q >= 0 || q <= 0.5f)
			{
				kernel_val *= (1 - 6.f*q*q + 6.f*q*q*q);
			}
			else if (q > 0.5f || q <= 1.0f)
			{
				kernel_val *= 2.0f*(1 - q)*(1 - q)*(1 - q);
			}
			else
			{
				kernel_val = 0;
			}

			kernel_values[particle*D_NR_OF_PARTICLES + ind] = kernel_val;
		}
	}
}


//TODO add kernelgradient
void SPH::correct_divergence_error(float* alpha)
{
	float div_i;
	float dens_i;
	float sum;
	int neighbor_index;
	float k_v_i[D_NR_OF_PARTICLES];

	calculate_kvi(alpha, &m_particles.vel, m_particles.mass, D_NR_OF_PARTICLES, m_delta_t, k_v_i);
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		dens_i = m_particles.dens[i];
		div_i = k_v_i[i] / dens_i;

		for (int j = 0; j < m_neighbor_data[i].n; ++j)
		{
			neighbor_index = m_neighbor_data[i].neighbor[j];
			sum += m_particles.mass[neighbor_index] * (div_i + k_v_i[neighbor_index] / m_particles.dens[neighbor_index]); //here
		}
		sum = m_delta_t * sum;
		m_particles.vel.x[i] -= sum;
		m_particles.vel.y[i] -= sum;
		m_particles.vel.z[i] -= sum;
	}
}

inline void calculate_kvi(float* alpha, Float3* vel, float* mass, int nr_particles, float delta_t, float* k_v_i)
{
	float d_dens = 0.0f;
	float x, y, z;
	for (int i = 0; i < nr_particles; ++i)
	{
		for (int n = 0; n < nr_particles; ++n)
		{
			x = (vel->x[n] - vel->x[i]); // here 
			y = (vel->y[n] - vel->y[i]); // here
			z = (vel->z[n] - vel->z[i]); // and here
			d_dens += mass[n] * (x + y + z);
		}

		k_v_i[i] = (1.f / delta_t)* d_dens *  alpha[i];
	}
}

//TODO: remake this function using the predicted velocity
void SPH::update_velocities(float dT)
{
	for (unsigned int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.vel.x[i] += m_particles.F_adv.x[i] * dT;
		m_particles.vel.y[i] += m_particles.F_adv.y[i] * dT;
		m_particles.pos.z[i] += m_particles.F_adv.z[i] * dT;
	}
}

SPH::~SPH()
{
	delete[] m_particles.F_adv.x;
	delete[] m_particles.F_adv.y;
	delete[] m_particles.F_adv.z;
	delete[] m_particles.vel.x;
	delete[] m_particles.vel.y;
	delete[] m_particles.vel.z;
	delete[] m_particles.pos.x;
	delete[] m_particles.pos.y;
	delete[] m_particles.pos.z;
	delete[] m_particles.dens;
	delete[] m_particles.mass;
	delete[] m_particles.p;

	delete[] m_neighbor_data;
}
