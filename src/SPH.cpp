#include "SPH.h"
#include <iostream>
#include <math.h>
#include <assert.h>

#define D_GRAVITY -9.82f
#define D_PI 3.1415926559f;
#define D_EPSILON 0.000000000000001f;

//since force is low a higher radius is requiered for small number of particles
#define D_RAD 0.08f;

SPH::SPH(int x, int y, int z)
{
	m_particles.dens = new float[D_NR_OF_PARTICLES];
	m_particles.pos.x = new float[D_NR_OF_PARTICLES];
	m_particles.pos.y = new float[D_NR_OF_PARTICLES];
	m_particles.pos.z = new float[D_NR_OF_PARTICLES];
	m_particles.vel.x = new float[D_NR_OF_PARTICLES];
	m_particles.vel.y = new float[D_NR_OF_PARTICLES];
	m_particles.vel.z = new float[D_NR_OF_PARTICLES];

	// predicted velocities
	m_particles.pred_vel.x = new float[D_NR_OF_PARTICLES];
	m_particles.pred_vel.y = new float[D_NR_OF_PARTICLES];
	m_particles.pred_vel.z = new float[D_NR_OF_PARTICLES];

	m_particles.F_adv.x = new float[D_NR_OF_PARTICLES];
	m_particles.F_adv.y = new float[D_NR_OF_PARTICLES];
	m_particles.F_adv.z = new float[D_NR_OF_PARTICLES];
	m_neighbor_data = new Neighbor_Data[D_NR_OF_PARTICLES];

	/* The radius should be read from a
	settings class with static values
	instead of being defined here. */
	m_particles.rad = 0.01f;
	m_particles.mass = 0.0042f;

	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i) {
		m_particles.dens[i] = 100.f;
		m_particles.pos.x[i] = 0.f;
		m_particles.pos.y[i] = 0.f;
		m_particles.pos.z[i] = 0.f;
		m_particles.vel.x[i] = 0.f;
		m_particles.vel.y[i] = 0.f;
		m_particles.vel.z[i] = 0.f;

		// predicted velocities
		m_particles.pred_vel.x[i] = 0.f;
		m_particles.pred_vel.y[i] = 0.f;
		m_particles.pred_vel.z[i] = 0.f;
	}

	init_positions(x, y, z, 5, 5);
}

SPH::~SPH()
{
	delete[] m_particles.F_adv.x;
	delete[] m_particles.F_adv.y;
	delete[] m_particles.F_adv.z;
	delete[] m_particles.vel.x;
	delete[] m_particles.vel.y;
	delete[] m_particles.vel.z;
	delete[] m_particles.pred_vel.x;
	delete[] m_particles.pred_vel.y;
	delete[] m_particles.pred_vel.z;
	delete[] m_particles.pos.x;
	delete[] m_particles.pos.y;
	delete[] m_particles.pos.z;
	delete[] m_particles.dens;
	delete[] m_neighbor_data;
}

void SPH::update(float dT)
{
	static float alpha[D_NR_OF_PARTICLES];
	static float k_v_i[D_NR_OF_PARTICLES];
	static Float3s f_tot[D_NR_OF_PARTICLES];
	static float scalar_values[D_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS];
	static float kernel_values[D_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS];
	
	find_neighborhoods();

	update_scalar_function(&m_particles.pos, m_neighbor_data, scalar_values);

	update_kernel_values(kernel_values, &m_particles.pos, m_neighbor_data);

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, m_neighbor_data, alpha, kernel_values);

	calculate_time_step(dT);

	non_pressure_forces();

	predict_velocities();

	//correct_density_error(alpha, m_delta_t, scalar_values, f_tot, k_v_i);

	update_positions();

	find_neighborhoods();  // t + dt

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, m_neighbor_data, alpha, kernel_values);

	correct_divergence_error(k_v_i, scalar_values, alpha);

	update_velocities();
}

void SPH::init_positions(int x_start, int y_start, int z_start, int rows, int cols) const
{
	float dist_between{ 1.2f * m_particles.rad };
	float padding_factor{ 1.4f };
	float x, y, z;
	int ind;

	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		x = i / (rows*cols);
		ind = i - x*rows*cols;
		y = ind / rows;
		z = ind % rows;

		m_particles.pos.y[i] = x_start + x * dist_between*padding_factor;
		m_particles.pos.x[i] = y_start + y * dist_between*padding_factor;
		m_particles.pos.z[i] = z_start + z * dist_between*padding_factor;
	}
}

void SPH::find_neighborhoods() const
{
	float vector_i_n[3];
	float dist_i_n_2;
	const float neigborhod_rad = D_RAD;
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
				dist_i_n_2 = vector_i_n[0]*vector_i_n[0] + vector_i_n[1]*vector_i_n[1] + vector_i_n[2]*vector_i_n[2];

				if (sqrt(dist_i_n_2) < neigborhod_rad)
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

void SPH::non_pressure_forces() const
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.F_adv.x[i] = 0.f;
		m_particles.F_adv.y[i] = m_particles.mass*D_GRAVITY;
		m_particles.F_adv.z[i] = 0.f;
	}
}

void SPH::calculate_time_step(float dT)
{
	float v_max_2 = 0;
	float x_2, y_2, z_2;

	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		x_2 = m_particles.vel.x[i] * m_particles.vel.x[i];
		y_2 = m_particles.vel.y[i] * m_particles.vel.y[i];
		z_2 = m_particles.vel.z[i] * m_particles.vel.z[i];
		if (v_max_2 < (x_2 + y_2 + z_2))
			v_max_2 = (x_2 + y_2 + z_2);
	}
	m_delta_t = 0.4f * (2.f * m_particles.rad)/sqrtf(v_max_2) - D_EPSILON;

	if (dT < m_delta_t)
		m_delta_t = dT;

	m_delta_t = m_delta_t < 0.00033f ? 0.00033f : m_delta_t;

}

void SPH::predict_velocities()
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.pred_vel.x[i] = m_particles.vel.x[i] + m_particles.F_adv.x[i] * m_delta_t / m_particles.mass;
		m_particles.pred_vel.y[i] = m_particles.vel.y[i] + m_particles.F_adv.y[i] * m_delta_t / m_particles.mass;
		m_particles.pred_vel.z[i] = m_particles.vel.z[i] + m_particles.F_adv.z[i] * m_delta_t / m_particles.mass;

		if (abs(m_particles.pos.x[i] + m_particles.pred_vel.x[i] * m_delta_t) >= 0.1)
		{
			m_particles.pred_vel.x[i] = 0.f;
		}
		if ((m_particles.pos.y[i] + m_particles.pred_vel.y[i] * m_delta_t) <= -0.5)
		{
			m_particles.pred_vel.y[i] = 0.f;
		}
		if (abs(m_particles.pos.z[i] + m_particles.pred_vel.z[i] * m_delta_t) >= 0.1)
		{
			m_particles.pred_vel.z[i] = 0.f;
		}
	}
}

void SPH::correct_density_error(float* alpha, float dT, float* scalar_values, Float3s* f_tot, float* k_v_i)
{
	static Float3s predicted_pressure[D_NR_OF_PARTICLES];
	static Float3s k[D_NR_OF_PARTICLES];

	int neighbor_index;
	int iter = 0;

	float dens_i;
	float dens_avg;
	float pressure_avg;
	float k_i_x, k_i_y, k_i_z;
	float sum_x = 0, sum_y = 0, sum_z = 0;
	float dx, dy, dz;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float ny;
	float dens_min, dens_max;

	do {
		dens_avg = 0.f;
		dens_min = INFINITY;
		dens_max = -1.f;

		calculate_pressure_force(f_tot, k_v_i, &m_particles.pos, m_particles.mass, scalar_values, m_neighbor_data, m_particles.dens);
		calculate_predicted_pressure(predicted_pressure, f_tot, m_particles.mass, m_particles.dens, scalar_values, m_delta_t, m_neighbor_data, &m_particles.pos);
		
		assert(m_delta_t != 0.0f, "deltaT");

		for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
		{
			k[i].x = predicted_pressure[i].x * alpha[i] / (m_delta_t*m_delta_t);
			k[i].y = predicted_pressure[i].y * alpha[i] / (m_delta_t*m_delta_t);
			k[i].z = predicted_pressure[i].z * alpha[i] / (m_delta_t*m_delta_t);
		}

		for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
		{
			dens_i = m_particles.dens[i];
		
			assert(dens_i != 0, "dens");
			
			dens_avg += dens_i;
			if (dens_min > dens_i) dens_min = dens_i;
			if (dens_max < dens_i) dens_max = dens_i;

			k_i_x = k[i].x / dens_i;
			k_i_y = k[i].y / dens_i;
			k_i_z = k[i].z / dens_i;
			for (int j = 0; j < m_neighbor_data[i].n; ++j)
			{
				neighbor_index = m_neighbor_data[i].neighbor[j];

				assert(m_particles.dens[neighbor_index] != 0.0f, "n dens");

				int linear_ind = j + D_MAX_NR_OF_NEIGHBORS*i;

				dx = m_particles.pos.x[neighbor_index] - m_particles.pos.x[i];
				dy = m_particles.pos.y[neighbor_index] - m_particles.pos.y[i];
				dz = m_particles.pos.z[neighbor_index] - m_particles.pos.z[i];

				kernel_gradient_x = dx * scalar_values[linear_ind];
				kernel_gradient_y = dy * scalar_values[linear_ind];
				kernel_gradient_z = dz * scalar_values[linear_ind];

				sum_x += m_particles.mass * (k_i_x + k[neighbor_index].x / m_particles.dens[neighbor_index]) * kernel_gradient_x;
				sum_y += m_particles.mass * (k_i_y + k[neighbor_index].y / m_particles.dens[neighbor_index]) * kernel_gradient_y;
				sum_z += m_particles.mass * (k_i_z + k[neighbor_index].z / m_particles.dens[neighbor_index]) * kernel_gradient_z;
			}
			m_particles.pred_vel.x[i] = m_particles.pred_vel.x[i] - m_delta_t*sum_x;
			m_particles.pred_vel.y[i] = m_particles.pred_vel.y[i] - m_delta_t*sum_y;
			m_particles.pred_vel.z[i] = m_particles.pred_vel.z[i] - m_delta_t*sum_z;
			sum_x = sum_z = sum_y = .0f;
		}

		// condition p_avg - p_0 > ny  ny = 1.01*(p_max_p_min) p_0 = 1000
		pressure_avg = dens_avg / D_NR_OF_PARTICLES - C_REST_DENS;
		ny = 1.01*(dens_max - dens_min);
		++iter;
	} while (pressure_avg - C_REST_DENS > ny || iter < 2);
}

void SPH::correct_strain_rate_error() {}

void SPH::update_positions() const
{
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.pos.x[i] += m_particles.pred_vel.x[i] * m_delta_t;
		m_particles.pos.y[i] += m_particles.pred_vel.y[i] * m_delta_t;
		m_particles.pos.z[i] += m_particles.pred_vel.z[i] * m_delta_t;
	}
}

void SPH::correct_divergence_error(float* stiffenss_velocity, float* scalar_values, float* alpha)
{
	float dens_i;
	float div_i;
	float pressure_acc_x, pressure_acc_y, pressure_acc_z;
	float dx, dy, dz;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	int neighbor_ind;
	static int iter = 0;
	//should not be here, send
	float mass = 0.00418f;
	float delta_dens_avg = calculate_stiffness(alpha, &m_particles.pred_vel, &m_particles.pos, m_delta_t, stiffenss_velocity, m_neighbor_data, scalar_values, m_particles.mass);
	++iter;
	while ((delta_dens_avg) > 1.0f && iter > 1)
	{
		for (auto particle_ind = 0; particle_ind < D_NR_OF_PARTICLES; ++particle_ind)
		{
			dens_i = m_particles.dens[particle_ind];

			assert(dens_i != 0.0f, "dens");

			div_i = stiffenss_velocity[particle_ind] / dens_i;
			pressure_acc_x = pressure_acc_y = pressure_acc_z = 0.0f;
			for (auto j = 0; j < m_neighbor_data[particle_ind].n; ++j)
			{
				int linear_ind = j + D_MAX_NR_OF_NEIGHBORS*particle_ind;
				neighbor_ind = m_neighbor_data[particle_ind].neighbor[j];

				assert(m_particles.dens[neighbor_ind] != 0.0f, "n dens");

				dx = m_particles.pos.x[neighbor_ind] - m_particles.pos.x[particle_ind];
				dy = m_particles.pos.y[neighbor_ind] - m_particles.pos.y[particle_ind];
				dz = m_particles.pos.z[neighbor_ind] - m_particles.pos.z[particle_ind];

				kernel_gradient_x = dx*scalar_values[linear_ind];
				kernel_gradient_y = dy*scalar_values[linear_ind];
				kernel_gradient_z = dz*scalar_values[linear_ind];

				pressure_acc_x += m_particles.mass * (div_i + stiffenss_velocity[neighbor_ind] /
					m_particles.dens[neighbor_ind]) * kernel_gradient_x;

				pressure_acc_y += m_particles.mass * (div_i + stiffenss_velocity[neighbor_ind] /
					m_particles.dens[neighbor_ind]) * kernel_gradient_y;

				pressure_acc_z += m_particles.mass * (div_i + stiffenss_velocity[neighbor_ind] /
					m_particles.dens[neighbor_ind]) * kernel_gradient_z;
			}
			//pressure_force_z/mass is not in report bot it is a force so should be a = F/m *delta_t
			m_particles.pred_vel.x[particle_ind] = m_particles.pred_vel.x[particle_ind] - m_delta_t * pressure_acc_x;
			m_particles.pred_vel.y[particle_ind] = m_particles.pred_vel.y[particle_ind] - m_delta_t * pressure_acc_y;
			m_particles.pred_vel.z[particle_ind] = m_particles.pred_vel.z[particle_ind] - m_delta_t * pressure_acc_z;
		}
		delta_dens_avg = calculate_stiffness(alpha, &m_particles.pred_vel, &m_particles.pos, m_delta_t, stiffenss_velocity, m_neighbor_data, scalar_values, m_particles.mass);

	}
}

void SPH::update_velocities()
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.vel.x[i] = m_particles.pred_vel.x[i];
		m_particles.vel.y[i] = m_particles.pred_vel.y[i];
		m_particles.vel.z[i] = m_particles.pred_vel.z[i];
	}
}

void update_density_and_factors(float mass, Float3* pos, float* dens, float* scalar_values,
	Neighbor_Data* neighbor_data, float* alpha, float* kernel_values)
{
	int nr_neighbors;
	float abs_sum_denom, sum_abs_denom = 0;
	float denom;
	float dx, dy, dz;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float sqrt_val;
	unsigned int neighbor_index;
	float x = 0.f, y = 0.f, z = 0.f;
	const float min_denom{ 0.000001f };
	int ind;

	for (auto particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		nr_neighbors = neighbor_data[particle].n;
		//added 1 * mass as particles own density as kernel is 1 at dist == 0
		//this should atleast be the case, but needs to be checked
		dens[particle] = mass;
		for (auto neighbor = 0; neighbor < nr_neighbors; ++neighbor)
		{
			neighbor_index = neighbor_data[particle].neighbor[neighbor];

			ind = neighbor + D_MAX_NR_OF_NEIGHBORS*particle;

			//Update density
			dens[particle] += kernel_values[ind] * mass;

			dx = pos->x[neighbor_index] - pos->x[particle];
			dy = pos->y[neighbor_index] - pos->y[particle];
			dz = pos->z[neighbor_index] - pos->z[particle];

			kernel_gradient_x = dx * mass * scalar_values[ind];
			kernel_gradient_y = dy * mass * scalar_values[ind];
			kernel_gradient_z = dz * mass * scalar_values[ind];

			x += kernel_gradient_x;
			y += kernel_gradient_y;
			z += kernel_gradient_z;

			sqrt_val = kernel_gradient_x*kernel_gradient_x + kernel_gradient_y*kernel_gradient_y
				+ kernel_gradient_z*kernel_gradient_z;

			float temp = sqrt(sqrt_val);

			sum_abs_denom += temp*temp;;
		}

		abs_sum_denom = sqrt(x*x + y*y + z*z);
		denom = (abs_sum_denom*abs_sum_denom + sum_abs_denom);

		// set alpha to max(denom,min_denom)
		denom = denom > min_denom ? denom : min_denom;
		alpha[particle] = dens[particle] / denom;

		sum_abs_denom = 0.f;
		x = y = z = 0.f;
	}
}

/*
 * Using a Cubic spline kernel
 * Divergence-Free SPH for Incompressible and Viscous Fluids, ref 6
 */
void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data)
{
	int ind;
	float x, y, z, q, q_2, length;
	float kernel_val = 0.f;
	float search_area = D_RAD;
	float pi = D_PI;

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

			length = sqrt(x*x + y*y + z*z);
			q = length / D_RAD;
			q_2 = q*q;

			// length is always equal or smaller to D_RAD => implicit intervall between [0, 1]
			kernel_val = (1.0f/(search_area*pi))*(1.0f - 1.5f*q_2 + 0.75f*q_2*q);

			kernel_values[particle*D_NR_OF_PARTICLES + neighbor] = kernel_val;
		}
	}
}

void calculate_pressure_force(Float3s* f_tot, float* k_v_i, Float3* pos, float mass, float* scalar_values, Neighbor_Data* neighbor_data, float* dens)
{
	int neighbor_ind;
	int linear_ind;
	int nr_of_neighbors;
	float sum_x = 0, sum_y = 0, sum_z = 0;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		nr_of_neighbors = neighbor_data[i].n;

		float kv_dens = k_v_i[i] / dens[i];

		for (int j = 0; j < nr_of_neighbors; ++j)
		{
			neighbor_ind = neighbor_data[i].neighbor[j];
			linear_ind = neighbor_ind + D_MAX_NR_OF_NEIGHBORS * i;

			x = pos->x[neighbor_ind] - pos->x[i];
			y = pos->y[neighbor_ind] - pos->y[i];
			z = pos->z[neighbor_ind] - pos->z[i];

			kernel_gradient_x = x * scalar_values[linear_ind];
			kernel_gradient_y = y * scalar_values[linear_ind];
			kernel_gradient_z = z * scalar_values[linear_ind];

			assert(dens[i] != 0.0f, "dens");
			sum_x += mass * (kv_dens + k_v_i[neighbor_ind] / dens[neighbor_ind]) * kernel_gradient_x;
			sum_y += mass * (kv_dens + k_v_i[neighbor_ind] / dens[neighbor_ind]) * kernel_gradient_y;
			sum_z += mass * (kv_dens + k_v_i[neighbor_ind] / dens[neighbor_ind]) * kernel_gradient_z;
		}
		f_tot[i].x = -mass * sum_x;
		f_tot[i].y = -mass * sum_y;
		f_tot[i].z = -mass * sum_z;
		sum_x = sum_y = sum_z = 0.f;
	}
}
void calculate_predicted_pressure(Float3s* predicted_pressure, Float3s* f_p, float mass, float_t*dens, float* scalar_value,
	float delta_t, Neighbor_Data* neighbor_data, Float3 * pos)
{
	int neighbor_length;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float res_x = 0.f, res_y = 0.f, res_z = 0.f;
	unsigned int neighbor_index;
	float d_t_2 = delta_t * delta_t;

	/*for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
	neighbor_length = neighbor_data[i].n;
	for (int j = 0; j < neighbor_length; ++j)
	{
	neighbor_index = neighbor_data[i].neighbor[j];
	int linear_ind = neighbor_index + D_MAX_NR_OF_NEIGHBORS*i;

	x = pos->x[neighbor_index] - pos->x[i];
	y = pos->y[neighbor_index] - pos->y[i];
	z = pos->z[neighbor_index] - pos->z[i];

	kernel_gradient_x = x * scalar_value[linear_ind];
	kernel_gradient_y = y * scalar_value[linear_ind];
	kernel_gradient_z = z * scalar_value[linear_ind];

	res_x += mass * (pred_vel->x[i] - pred_vel->x[neighbor_index]) * kernel_gradient_x;
	res_y += mass * (pred_vel->y[i] - pred_vel->y[neighbor_index]) * kernel_gradient_y;
	res_z += mass * (pred_vel->z[i] - pred_vel->z[neighbor_index]) * kernel_gradient_z;
	}
	predicted_pressure[i].x += delta_t * res_x;
	predicted_pressure[i].y += delta_t * res_y;
	predicted_pressure[i].z += delta_t * res_z;
	res_x = res_y = res_z = 0.f;
	}
	*/
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		neighbor_length = neighbor_data[i].n;
		for (int j = 0; j < neighbor_length; ++j)
		{
			neighbor_index = neighbor_data[i].neighbor[j];
			int linear_ind = neighbor_index + D_MAX_NR_OF_NEIGHBORS * i;

			x = pos->x[neighbor_index] - pos->x[i];
			y = pos->y[neighbor_index] - pos->y[i];
			z = pos->z[neighbor_index] - pos->z[i];

			kernel_gradient_x = x * scalar_value[linear_ind];
			kernel_gradient_y = y * scalar_value[linear_ind];
			kernel_gradient_z = z * scalar_value[linear_ind];

			res_x += mass * (f_p[i].x / dens[i] - f_p[neighbor_index].x / dens[i]) * kernel_gradient_x;
			res_y += mass * (f_p[i].y / dens[i] - f_p[neighbor_index].y / dens[i]) * kernel_gradient_y;
			res_z += mass * (f_p[i].z / dens[i] - f_p[neighbor_index].z / dens[i]) * kernel_gradient_z;
		}
		predicted_pressure[i].x = d_t_2 *res_x;
		predicted_pressure[i].y = d_t_2 *res_y;
		predicted_pressure[i].z = d_t_2 *res_z;
		res_x = res_y = res_z = 0.f;
	}
}

float calculate_stiffness(float* alpha, Float3* pred_vel, Float3* pos, float delta_t, 
	float* k_v_i, Neighbor_Data* neighbor_data, float* scalar_value, float mass)
{
	unsigned int neighbor_index;
	int neighbor_length;

	float dx = 0.f, dy = 0.f, dz = 0.f;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float d_dens;
	float d_dens_avg = 0.0f;

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		neighbor_length = neighbor_data[i].n;
		for (int j = 0; j < neighbor_length; ++j)
		{
			neighbor_index = neighbor_data[i].neighbor[j];
			int linear_ind = neighbor_index + D_MAX_NR_OF_NEIGHBORS*i;

			x = pos->x[neighbor_index] - pos->x[i];
			y = pos->y[neighbor_index] - pos->y[i];
			z = pos->z[neighbor_index] - pos->z[i];

			kernel_gradient_x = x*scalar_value[linear_ind];
			kernel_gradient_y = y*scalar_value[linear_ind];
			kernel_gradient_z = z*scalar_value[linear_ind];

			dx += kernel_gradient_x*(pred_vel->x[i] - pred_vel->x[neighbor_index]) ;
			dy += kernel_gradient_y*(pred_vel->y[i] - pred_vel->y[neighbor_index]) ;
			dz += kernel_gradient_z*(pred_vel->z[i] - pred_vel->z[neighbor_index]) ;
		}
		d_dens = (dx + dy + dz)*delta_t;
		d_dens_avg += d_dens;
		
		k_v_i[i] = (1.f / delta_t)*d_dens*alpha[i];
		dx = dy = dz = 0.f;
	}
	return d_dens_avg / D_NR_OF_PARTICLES;
}

/*
* Derived the Cubic spline kernel by q
* Divergence-Free SPH for Incompressible and Viscous Fluids, section 4.2
*/
void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values)
{
	int neighbor_ind;

	float kernel_derive, scalar_value;
	float search_area = D_RAD;
	float pi = D_PI;
	float q, q_2, dist;
	float dx, dy, dz;

	//Loop through all particles
	for (auto particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		//Loop through particle neibors
		for (auto neighbor = 0; neighbor < neighbor_data[particle].n; ++neighbor)
		{
			neighbor_ind = neighbor_data[particle].neighbor[neighbor];

			dx = pos->x[particle] - pos->x[neighbor_ind];
			dy = pos->y[particle] - pos->y[neighbor_ind];
			dz = pos->z[particle] - pos->z[neighbor_ind];

			dist = sqrt(dx*dx + dy*dy + dz*dz);
			q = dist / D_RAD;
			q_2 = q*q;

			// length is always equal or smaller to D_RAD => implicit intervall between [0, 1]
			kernel_derive = (1.0f / (search_area*pi))*(-3.0f*q + 2.25f*q_2);

			scalar_value = kernel_derive*(1.0f /  (search_area*dist));

			scalar_values[particle*D_MAX_NR_OF_NEIGHBORS + neighbor] = scalar_value;
		}
	}

}