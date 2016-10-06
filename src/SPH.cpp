#include "SPH.h"

#include <math.h>

#define D_GRAVITY -9.82f
#define D_PI 3.1415926559f;
#define D_EPSILON 0.000000000000001f;

SPH::SPH()
{
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

void SPH::update(float dT)
{
	static float alpha[D_NR_OF_PARTICLES];
	static Float3s k_v_i[D_NR_OF_PARTICLES];
	static Float3s f_tot[D_NR_OF_PARTICLES];
	static float scalar_values[D_NR_OF_PARTICLES * D_NR_OF_PARTICLES];
	static float kernel_values[D_NR_OF_PARTICLES*D_NR_OF_PARTICLES];

	find_neighborhoods();

	update_scalar_function(&m_particles.pos, m_neighbor_data, scalar_values, C_NEIGHBOR_RAD);

	update_kernel_values(kernel_values, &m_particles.pos, m_neighbor_data, C_NEIGHBOR_RAD);

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, D_NR_OF_PARTICLES, m_neighbor_data, alpha, kernel_values);

	calculate_kvi(alpha, &m_particles.vel, &m_particles.pos, m_particles.mass, D_NR_OF_PARTICLES, m_delta_t, k_v_i, m_neighbor_data, scalar_values);

	non_pressure_forces();

	calculate_time_step(dT);

	predict_velocities(dT);

	correct_density_error(alpha, dT, scalar_values, f_tot, k_v_i);

	update_positions(dT);

	find_neighborhoods();  // t + dt

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, D_NR_OF_PARTICLES, m_neighbor_data, alpha, kernel_values);

	// TODO: this function is incorretly, implemented correct
	//correct_divergence_error(alpha);

	update_velocities(dT);
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

void SPH::pressure_forces() const
{
	// calculate the particle pressure
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.p[i] = m_particles.dens[i] - C_REST_DENS;
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

void SPH::calculate_time_step(float dT)
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

	if (dT < m_delta_t)
		m_delta_t = dT;
}

void SPH::predict_velocities(float dT)
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.vel.x[i] += m_particles.F_adv.x[i] * m_delta_t;
		m_particles.vel.y[i] += m_particles.F_adv.y[i] * m_delta_t;
		m_particles.pos.z[i] += m_particles.F_adv.z[i] * m_delta_t;
	}
}

//TODO: add kernel function
void SPH::correct_density_error(float* alpha, float dT, float* g_values, Float3s* f_tot, Float3s* k_v_i)
{
	static Float3s predicted_pressure[D_NR_OF_PARTICLES];
	static Float3s k[D_NR_OF_PARTICLES];

	float dens_i;
	float k_i_x;
	float k_i_y;
	float k_i_z;
	float sum_x{ 0.f };
	float sum_y{ 0.f };
	float sum_z{ 0.f };
	int neighbor_index;

	calculate_pressure_force(f_tot, k_v_i, &m_particles.pos, m_particles.mass, g_values, m_neighbor_data, m_particles.dens);
	//if (f_tot[0].y != 0)std::cout << f_tot[0].y << std::endl;
	calculate_predicted_pressure(predicted_pressure, f_tot, m_particles.mass, m_particles.dens, g_values, m_delta_t, m_neighbor_data, &m_particles.pos, C_REST_DENS);
	//std::cout << f_tot[0].x

	for (int i = 0; i < D_NR_OF_PARTICLES; i++)
	{
		k[i].x = predicted_pressure[i].x * alpha[i] / (m_delta_t*m_delta_t);
		k[i].y = predicted_pressure[i].y * alpha[i] / (m_delta_t*m_delta_t);
		k[i].z = predicted_pressure[i].z * alpha[i] / (m_delta_t*m_delta_t);
	}

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		dens_i = m_particles.dens[i];
		k_i_x = k[i].x / dens_i;
		k_i_y = k[i].y / dens_i;
		k_i_z = k[i].z / dens_i;

		for (int j = 0; j < m_neighbor_data[i].n; ++j)
		{
			neighbor_index = m_neighbor_data[i].neighbor[j];
			sum_x += m_particles.mass[neighbor_index] * (k_i_x + k[neighbor_index].x / m_particles.dens[neighbor_index]); //here
			sum_y += m_particles.mass[neighbor_index] * (k_i_y + k[neighbor_index].y / m_particles.dens[neighbor_index]); //here
			sum_z += m_particles.mass[neighbor_index] * (k_i_z + k[neighbor_index].z / m_particles.dens[neighbor_index]); //here

		}

		m_particles.vel.x[i] -= m_delta_t * sum_x;
		m_particles.vel.y[i] -= m_delta_t * sum_y;
		m_particles.vel.z[i] -= m_delta_t * sum_z;
		sum_x = .0f;
		sum_y = .0f;
		sum_z = .0f;
	}
}

void SPH::correct_strain_rate_error() {}

void SPH::update_positions(float dT) const
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.pos.x[i] += m_particles.vel.x[i] * m_delta_t;
		m_particles.pos.y[i] += m_particles.vel.y[i] * m_delta_t;
		m_particles.pos.z[i] += m_particles.vel.z[i] * m_delta_t;
	}
	//std::cout << "pos particle 0 : " << m_particles.pos.x[0] << " " << m_particles.pos.y[0] << " " << m_particles.pos.z[0] << std::endl;
}

//TODO add kernelgradient
void SPH::correct_divergence_error(float* alpha, Float3s* k_v_i)
{
	float dens_i;
	float div_i_x;
	float div_i_y;
	float div_i_z;
	float sum_x{ 0.f };
	float sum_y{ 0.f };
	float sum_z{ 0.f };

	int neighbor_index;

	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		dens_i = m_particles.dens[i];
		div_i_x = k_v_i[i].x / dens_i;
		div_i_y = k_v_i[i].y / dens_i;
		div_i_z = k_v_i[i].z / dens_i;

		for (auto j = 0; j < m_neighbor_data[i].n; ++j)
		{
			neighbor_index = m_neighbor_data[i].neighbor[j];
			sum_x += m_particles.mass[neighbor_index] * (div_i_x + k_v_i[neighbor_index].x / m_particles.dens[neighbor_index]); //here
			sum_y += m_particles.mass[neighbor_index] * (div_i_y + k_v_i[neighbor_index].y / m_particles.dens[neighbor_index]); //here
			sum_z += m_particles.mass[neighbor_index] * (div_i_z + k_v_i[neighbor_index].z / m_particles.dens[neighbor_index]); //here

		}
		m_particles.vel.x[i] -= m_delta_t * sum_x;
		m_particles.vel.y[i] -= m_delta_t * sum_y;
		m_particles.vel.z[i] -= m_delta_t * sum_z;
		sum_x = .0f;
		sum_y = .0f;
		sum_z = .0f;
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

inline void update_density_and_factors(float* mass, Float3* pos, float* dens, float* scalar_values, float nr_particles,
	Neighbor_Data* neighbor_data, float* alpha, float* kernel_values)
{
	int nr_neighbors;
	float abs_sum_denom{ 0 };
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

			int linear_ind = neighbor_index + D_NR_OF_PARTICLES*particle;

			//Update density
			dens[particle] += neighbor_mass*kernel_values[linear_ind];

			dx = pos->x[neighbor_index] - pos->x[particle];
			dy = pos->y[neighbor_index] - pos->y[particle];
			dz = pos->z[neighbor_index] - pos->y[particle];

			kernel_gradient_x = dx * scalar_values[linear_ind];
			kernel_gradient_y = dy * scalar_values[linear_ind];
			kernel_gradient_z = dz * scalar_values[linear_ind];

			std::cout << scalar_values[linear_ind] << '\n';

			sqrt_val = (kernel_gradient_x*kernel_gradient_x + kernel_gradient_y*kernel_gradient_y
				+ kernel_gradient_z*kernel_gradient_z);

			sum_abs_denom = neighbor_mass*neighbor_mass*sqrt(sqrt_val);
			abs_sum_denom += sum_abs_denom*sum_abs_denom;
		}
		denom = (sum_abs_denom*sum_abs_denom + abs_sum_denom);

		// set alpha to max(denom,min_denom)
		denom = denom > min_denom ? denom : min_denom;
		alpha[particle] = dens[particle] / denom;
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
			q = sqrt(len) / NEIGHBOR_RAD;

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

inline void calculate_pressure_force(Float3s* f_tot, Float3s* k_v_i, Float3* pos, float* mass, float* g_val, Neighbor_Data* neighbor_data, float* dens)
{
	unsigned int neighbor_index;
	unsigned int n_neighbors;
	float sum_x{ 0.f };
	float sum_y{ 0.f };
	float sum_z{ 0.f };
	float x{ 0.f };
	float y{ 0.f };
	float z{ 0.f };
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		n_neighbors = neighbor_data[i].n;

		for (int j = 0; j < n_neighbors; ++j)
		{
			neighbor_index = neighbor_data[i].neighbor[j] + i * D_MAX_NR_OF_NEIGHBORS;

			x = pos->x[neighbor_index];
			y = pos->y[neighbor_index];
			z = pos->z[j];

			kernel_gradient_x = x * g_val[j + D_MAX_NR_OF_NEIGHBORS*i];
			kernel_gradient_y = y * g_val[j + D_MAX_NR_OF_NEIGHBORS*i];
			kernel_gradient_z = z * g_val[j + D_MAX_NR_OF_NEIGHBORS*i];

			sum_x += mass[neighbor_index] * (k_v_i[i].x / (dens[i] + 1.000f) + k_v_i[neighbor_index].x / (dens[neighbor_index] + 1.000f))*kernel_gradient_x;
			sum_y += mass[neighbor_index] * (k_v_i[i].y / (dens[i] + 1.000f) + k_v_i[neighbor_index].y / (dens[neighbor_index] + 1.000f))*kernel_gradient_y;
			sum_z += mass[neighbor_index] * (k_v_i[i].z / (dens[i] + 1.000f) + k_v_i[neighbor_index].z / (dens[neighbor_index] + 1.000f))*kernel_gradient_z;
			//std::cout << dens[i] << " :  dens[i]" << std::endl;
		}
		f_tot[i].x = -mass[i] * sum_x;
		f_tot[i].y = -mass[i] * sum_y;
		f_tot[i].z = -mass[i] * sum_z;
		sum_x = sum_y = sum_z = 0.f;
	}
}

//TODO: add gradient
inline void calculate_predicted_pressure(Float3s* predicted_pressure, Float3s* f_p, float* mass, float_t*dens, float* g_val, float delta_t, Neighbor_Data* n_data, Float3 * pos, const float rest_dens)
{
	int neighbor_length;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float d_t_2 = delta_t*delta_t;
	float res_x{ 0.f };
	float res_y{ 0.f };
	float res_z{ 0.f };;
	unsigned int neighbor_index;

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		neighbor_length = n_data[i].n;
		for (int j = 0; j < neighbor_length; ++j)
		{
			neighbor_index = n_data[i].neighbor[j] + i * D_MAX_NR_OF_NEIGHBORS;

			x = pos->x[j];
			y = pos->y[j];
			z = pos->z[j];

			kernel_gradient_x = x * g_val[j + D_MAX_NR_OF_NEIGHBORS*i];
			kernel_gradient_y = y * g_val[j + D_MAX_NR_OF_NEIGHBORS*i];
			kernel_gradient_z = z * g_val[j + D_MAX_NR_OF_NEIGHBORS*i];

			res_x += mass[j] * (f_p[i].x / dens[i] - f_p[j].x / dens[i])*kernel_gradient_x;
			res_y += mass[j] * (f_p[i].y / dens[i] - f_p[j].y / dens[i])*kernel_gradient_y;
			res_z += mass[j] * (f_p[i].z / dens[i] - f_p[j].z / dens[i])*kernel_gradient_z;
		}
		predicted_pressure[i].x = d_t_2 *res_x + rest_dens;
		predicted_pressure[i].y = d_t_2 *res_y + rest_dens;
		predicted_pressure[i].z = d_t_2 *res_z + rest_dens;
	}
}

//TODO: add gradient kernal
inline void calculate_kvi(float* alpha, Float3* vel, Float3* pos, float* mass, int nr_particles,
	float delta_t, Float3s* k_v_i, Neighbor_Data* neighbor_data, float* scalar_value)
{
	float x, y, z;
	unsigned int neighbor_index;
	unsigned int n_neighbors;
	float particle_mass;
	float kernel_gradient_x;
	float kernel_gradient_y;
	float kernel_gradient_z;
	float d_dens_x{ 0.f };
	float d_dens_y{ 0.f };
	float d_dens_z{ 0.f };

	for (int i = 0; i < nr_particles; ++i)
	{
		n_neighbors = neighbor_data[i].n;
		for (int n = 0; n < n_neighbors; ++n)

		{
			neighbor_index = neighbor_data[i].neighbor[n] + i * D_MAX_NR_OF_NEIGHBORS;

			particle_mass = mass[neighbor_index];
			x = pos->x[neighbor_index];
			y = pos->y[neighbor_index];
			z = pos->z[neighbor_index];

			kernel_gradient_x = x * scalar_value[n + D_MAX_NR_OF_NEIGHBORS*i];
			kernel_gradient_y = y * scalar_value[n + D_MAX_NR_OF_NEIGHBORS*i];
			kernel_gradient_z = z * scalar_value[n + D_MAX_NR_OF_NEIGHBORS*i];

			x = (vel->x[i] - vel->x[n])*kernel_gradient_x;
			y = (vel->y[i] - vel->y[n])*kernel_gradient_y;
			z = (vel->z[i] - vel->z[n])*kernel_gradient_z;

			d_dens_x += mass[n] * (x);
			d_dens_y += mass[n] * (y);
			d_dens_z += mass[n] * (z);
		}
		k_v_i[i].x = (1.f / delta_t)* d_dens_x *  alpha[i];
		k_v_i[i].y = (1.f / delta_t)* d_dens_y *  alpha[i];
		k_v_i[i].z = (1.f / delta_t)* d_dens_z *  alpha[i];

		//if(k_v_i[i].y =! 0.f)std::cout << k_v_i[i].y << "k : i " << i << std::endl;

		d_dens_x = d_dens_y = d_dens_z = 0.f;
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
				kernel_derive = 0.0f;
			}
			g[neighbor_ind + particle*D_NR_OF_PARTICLES] = kernel_derive / (dist * NEIGHBOR_RADIUS);
		}
	}
}