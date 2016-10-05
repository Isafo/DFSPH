#include "SPH.h"

#include <math.h>

#define D_GRAVITY -9.82f
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

	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i) {
		m_particles.dens[i] = 0.f;
		m_particles.mass[i] = 1.f;
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
	static float alpha[D_NR_OF_PARTICLES];
	static float scalar_values[D_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS];
	static float kernel_values[D_NR_OF_PARTICLES*D_MAX_NR_OF_NEIGHBORS];

	find_neighborhoods();

	update_scalar_function( &m_particles.pos, m_neighbor_data, scalar_values, C_NEIGHBOR_RAD);

	update_kernel_values(kernel_values, &m_particles.pos, m_neighbor_data, C_NEIGHBOR_RAD);

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, D_NR_OF_PARTICLES, m_neighbor_data, alpha, kernel_values);

	non_pressure_forces();

	calculate_time_step();

	predict_velocities(dT);

	//correct_density_error();

	update_positions(dT);

	find_neighborhoods();  // t + dt

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, D_NR_OF_PARTICLES, m_neighbor_data, alpha, kernel_values);

	// TODO: this function is incorretly implemented correct
	//correct_divergence_error(alpha);

	update_velocities(dT);
}

void SPH::find_neighborhoods() const
{
	float vector_i_n[3];
	float dist_i_n_2;
	const float neigborhod_rad_2{ C_NEIGHBOR_RAD };
	int count{ 0 };
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		for (auto n = 0; n < D_NR_OF_PARTICLES; ++n)
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


inline void update_density_and_factors(float* mass, Float3* pos, float* dens, float* scalar_values, float nr_particles, 
										Neighbor_Data* neighbor_data, float* alpha, float* kernel_values)
{
	int nr_neighbors;
	float abs_sum_denom{ 0 };
	float temp;

	float sum_abs_denom{ 0 };
	float denom;

	float neighbor_mass;
	float dx, dy, dz;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float sqrt_val;
	unsigned int neighbor_index;

	const float min_denom{ 0.000001f };

	for (auto particle = 0; particle < nr_particles; ++particle)
	{
		nr_neighbors = neighbor_data[particle].n;
		for (auto neighbor = 0; neighbor < nr_neighbors; ++neighbor)
		{
			neighbor_index = neighbor_data[particle].neighbor[neighbor];
			neighbor_mass = mass[neighbor_index];

			int linear_ind = neighbor + D_MAX_NR_OF_NEIGHBORS*particle;

			//Update density
			dens[particle] += neighbor_mass*kernel_values[linear_ind];

			dx = pos->x[neighbor_index] - pos->x[particle];
			dy = pos->y[neighbor_index] - pos->y[particle];
			dz = pos->z[neighbor_index] - pos->y[particle];

			kernel_gradient_x = dx * scalar_values[linear_ind];
			kernel_gradient_y = dy * scalar_values[linear_ind];
			kernel_gradient_z = dz * scalar_values[linear_ind];

			sqrt_val = (kernel_gradient_x*kernel_gradient_x + kernel_gradient_y*kernel_gradient_y
				+ kernel_gradient_z*kernel_gradient_z);

			sum_abs_denom = neighbor_mass*neighbor_mass*sqrt(sqrt_val);
			abs_sum_denom += sum_abs_denom*sum_abs_denom;
		}
		denom = sum_abs_denom*sum_abs_denom + abs_sum_denom;

		// set alpha to max(denom,min_denom)
		denom = denom > min_denom ? denom : min_denom;
		alpha[particle] = dens[particle] / denom;
	}
}

void SPH::init_positions(glm::vec3* start_pos, int rows, int cols) const
{
	float dist_between{ 2.f * m_particles.rad };
	float padding_factor{ 1.1f };

	float x, y, z;
	int ind;

	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
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

void SPH::non_pressure_forces() const
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.F_adv.x[i] = 0.f;
		m_particles.F_adv.y[i] = D_GRAVITY;
		m_particles.F_adv.z[i] = 0.f;
	}
}
void SPH::pressure_forces() const
{
	// calculate the particle pressure
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.p[i] = m_particles.dens[i] - C_REST_DENS;
	}
}
void SPH::calculate_time_step()
{
	float v_max_2 = 0;
	float x, y, z;
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		x = m_particles.pos.x[i] * m_particles.pos.x[i];
		y = m_particles.pos.y[i] * m_particles.pos.y[i];
		z = m_particles.pos.z[i] * m_particles.pos.z[i];
		if (v_max_2 < (x + y + z))v_max_2 = (x + y + z);
	}
	m_delta_t = 0.4f*(m_particles.rad * 2) / sqrtf(v_max_2) - D_EPSILON;
}

void SPH::predict_velocities(float dT) const
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.vel.x[i] += m_particles.F_adv.x[i] * m_delta_t;
		m_particles.vel.y[i] += m_particles.F_adv.y[i] * m_delta_t;
		m_particles.pos.z[i] += m_particles.F_adv.z[i] * m_delta_t;
	}
}

//TODO: add kernel function
void SPH::correct_density_error(float* alpha) { }

void SPH::correct_strain_rate_error() {}

void SPH::update_positions(float dT) const
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.pos.x[i] += m_particles.vel.x[i] * m_delta_t;
		m_particles.pos.y[i] += m_particles.vel.y[i] * m_delta_t;
		m_particles.pos.z[i] += m_particles.vel.z[i] * m_delta_t;
	}
}

inline void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* g, const float NEIGHBOR_RADIUS)
{
	float kernel_derive;
	//Loop through all particles
	for (auto particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		//Loop through particle neibors
		for (auto neighbor = 0; neighbor < neighbor_data[particle].n; ++neighbor)
		{
			auto neighbor_ind = neighbor_data[particle].neighbor[neighbor];

			auto dx = pos->x[particle] - pos->x[neighbor_ind];
			auto dy = pos->y[particle] - pos->y[neighbor_ind];
			auto dz = pos->z[particle] - pos->z[neighbor_ind];

			auto dist = std::sqrt(dx*dx + dy*dy + dz*dz);
			auto q = dist / NEIGHBOR_RADIUS;

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
			g[neighbor + particle*D_MAX_NR_OF_NEIGHBORS] = kernel_derive / (dist * NEIGHBOR_RADIUS);
		}
	}
}

inline void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data, const float NEIGHBOR_RAD)
{
	unsigned int ind;
	float x, y, z, q;
	float kernel_val;
	for (auto particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		float particle_pos_x = pos->x[particle];
		float particle_pos_y = pos->y[particle];
		float particle_pos_z = pos->z[particle];
		for (auto neighbor = 0; neighbor < neighbor_data[particle].n; ++neighbor)
		{
			// compute q
			ind = neighbor_data[particle].neighbor[neighbor];
			x = pos->x[ind] - particle_pos_x;
			y = pos->y[ind] - particle_pos_y;
			z = pos->z[ind] - particle_pos_z;

			float len = x*x + y*y + z*z;
			q = sqrt(x*x + y*y + z*z) / NEIGHBOR_RAD;

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
	float sum{ 0.f };
	int neighbor_index;
	float k_v_i[D_NR_OF_PARTICLES];

	calculate_kvi(alpha, &m_particles.vel, m_particles.mass, D_NR_OF_PARTICLES, m_delta_t, k_v_i);
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		dens_i = m_particles.dens[i];
		div_i = k_v_i[i] / dens_i;

		for (auto j = 0; j < m_neighbor_data[i].n; ++j)
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

void SPH::update_velocities(float dT) const
{
}

inline void calculate_kvi(float* alpha, Float3* vel, float* mass, int nr_particles, float delta_t, float* k_v_i)
{
	float d_dens{ 0.0f };
	float x, y, z;
	for (auto i = 0; i < nr_particles; ++i)
	{
		for (auto n = 0; n < nr_particles; ++n)
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
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.vel.x[i] += m_particles.F_adv.x[i] * m_delta_t;
		m_particles.vel.y[i] += m_particles.F_adv.y[i] * m_delta_t;
		m_particles.vel.z[i] += m_particles.F_adv.z[i] * m_delta_t;

		if (abs(m_particles.pos.x[i] + m_particles.vel.x[i] * m_delta_t) >= 0.5)
		{
			m_particles.vel.x[i] = 0.0f;
		}

		if (m_particles.pos.y[i] + m_particles.vel.y[i] * m_delta_t <= -0.5)
		{
			m_particles.vel.y[i] = 0.0f;
		}

		if (abs(m_particles.pos.z[i] + m_particles.vel.z[i] * m_delta_t) >= 0.5)
		{
			m_particles.vel.z[i] = 0.0f;
		}
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
