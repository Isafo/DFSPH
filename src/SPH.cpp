#include "SPH.h"

#define D_GRAVITY -9.82f
#define D_NEIGHBOR_RAD 0.3f;
#define D_PI 3.1415926559f;
#define D_EPSILON 0.000000000000001f;

SPH::SPH()
{
}

SPH::SPH(int size) {
	m_nr_of_particles = size;

	m_particles.alpha = new float[m_nr_of_particles];
	m_particles.dens = new float[m_nr_of_particles];
	m_particles.mass = new float[m_nr_of_particles];
	m_particles.p = new float[m_nr_of_particles];
	m_particles.k_v_i = new float[m_nr_of_particles];
	m_particles.pos.x = new float[m_nr_of_particles];
	m_particles.pos.y = new float[m_nr_of_particles];
	m_particles.pos.z = new float[m_nr_of_particles];
	m_particles.vel.x = new float[m_nr_of_particles];
	m_particles.vel.y = new float[m_nr_of_particles];
	m_particles.vel.z = new float[m_nr_of_particles];
	m_particles.F_adv.x = new float[m_nr_of_particles];
	m_particles.F_adv.y = new float[m_nr_of_particles];
	m_particles.F_adv.z = new float[m_nr_of_particles];
	m_neighbor_data = new Neighbor_Data[m_nr_of_particles];

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
		m_neighbor_data[i].neighbor = new int[100];
		m_neighbor_data[i].g_value = new float[100];
	}
}

void SPH::update(float dT) 
{
	find_neighborhoods();
	calculate_densities();
	calculate_factors();
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
	for (int i = 0; i < m_nr_of_particles; ++i)
	{
		for (int n = 0; n < m_nr_of_particles; ++n)
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
	for (int i = 0; i < m_nr_of_particles; ++i)
	{	
		n_neighbors = m_neighbor_data[i].n;
		for (int n = 0; n < n_neighbors ; ++n)
		{
			m_particles.dens[i] += m_particles.mass[m_neighbor_data[i].neighbor[n]]; // * kernel function
		}
	}
}

// TODO: Add gradient kernal function
void SPH::calculate_factors()
{
	int nr_neighbors;
	float abs_sum_denom = 0;
	float temp;
	float sum_abs_denom = 0;
	float denom;

	for (int i = 0; i < m_nr_of_particles; i++)
	{
		nr_neighbors = m_neighbor_data[i].n;
		for (int n = 0; n < nr_neighbors; ++n)
		{
			temp = m_particles.mass[m_neighbor_data[i].neighbor[n]];// * gradient kernel function
			sum_abs_denom += temp; 
			abs_sum_denom += temp*temp;
		}
		//this addition is suggested in the report as an alternative to clamping
		denom = D_EPSILON + sum_abs_denom*sum_abs_denom + abs_sum_denom;
		m_particles.alpha[i] = m_particles.dens[i]/denom;
	}
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
	for (int i = 0; i < m_nr_of_particles; i++) 
	{
		m_particles.F_adv.x[i] = 0.f;
		m_particles.F_adv.y[i] = D_GRAVITY;
		m_particles.F_adv.z[i] = 0.f;
	}
}
void SPH::pressure_forces()
{
	// calculate the particle pressure
	for (int i = 0; i < m_nr_of_particles; ++i)
	{
		m_particles.p[i] = m_particles.dens[i] - C_REST_DENS;
	}
}
void SPH::calculate_time_step() 
{
	float v_max_2 = 0;
	float x, y, z;
	for (int i = 0; i < m_nr_of_particles; ++i)
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
	for (unsigned int i = 0; i < m_nr_of_particles; ++i)
	{
		m_particles.vel.x[i] += m_particles.F_adv.x[i] * dT;
		m_particles.vel.y[i] += m_particles.F_adv.y[i] * dT;
		m_particles.pos.z[i] += m_particles.F_adv.z[i] * dT;
	}
}

//TODO: add kernel function
void SPH::correct_density_error()
{

}

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

void SPH::update_function_g()
{
	float kernel_derive;
	//Loop through all particles
	for (int particle = 0; particle < m_nr_of_particles; ++particle)
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
/*
float SPH::kernel(const float q) const
{
	auto kernel = 1 / PI;
	if (q >= 0 || q <= 0.5f)
	{
		kernel *= (1 - 6.f*q*q + 6.f*q*q*q);
	}
	else if (q > 0.5f || q <= 1.0f)
	{
		kernel *= 2.0f*(1 - q)*(1 - q)*(1 - q);
	}
	else
	{
		kernel = 0;
	}

	return kernel;
}
*/

//TODO add kernelgradient
void SPH::correct_divergence_error()
{
	float k_v_i;
	float div_i;
	float dens_i;
	float sum;
	int neighbor_index;
	calculate_kvi();
	for (int i = 0; i < m_nr_of_particles; ++i)
	{
		k_v_i = m_particles.k_v_i[i];
		dens_i = m_particles.dens[i];
		div_i = k_v_i / dens_i;

		for (int j = 0; j < m_neighbor_data[i].n; ++j)
		{
			neighbor_index = m_neighbor_data[i].neighbor[j];
			sum += m_particles.mass[neighbor_index] * (div_i + m_particles.k_v_i[neighbor_index] / m_particles.dens[neighbor_index]); //here
		}
		sum = m_delta_t * sum;
		m_particles.vel.x[i] -= sum;
		m_particles.vel.y[i] -= sum;
		m_particles.vel.z[i] -= sum;
	}
	

}
void SPH::calculate_kvi()
{
	float d_dens;
	float x, y, z;
	for (int i = 0; i < m_nr_of_particles; ++i)
	{
		for (int n = 0; n < m_nr_of_particles; ++n)
		{
			x = (m_particles.vel.x[n] - m_particles.vel.x[i]); // here 
			y = (m_particles.vel.y[n] - m_particles.vel.y[i]); // here
			z = (m_particles.vel.z[n] - m_particles.vel.z[i]); // and here
			d_dens += m_particles.mass[n] * (x + y + z);

		}

		m_particles.k_v_i[i] = (1.f / m_delta_t)* d_dens *  m_particles.alpha[i];
	}
}

//TODO: remake this function using the predicted velocity
void SPH::update_velocities(float dT)
{
	for (unsigned int i = 0; i < m_nr_of_particles; ++i)
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
	delete[] m_particles.alpha;
	delete[] m_particles.dens;
	delete[] m_particles.mass;
	delete[] m_particles.p;

	delete[] m_neighbor_data->g_value;
	delete[] m_neighbor_data->neighbor;
	delete[] m_neighbor_data;
	
}
