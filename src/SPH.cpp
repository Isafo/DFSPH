#include "SPH.h"
#include <iostream>
#include <math.h>
#include <assert.h>

//GRAVITY IS SO DAMN WRONG, other forces are to small atm to compensate for the
//acceleration coused by gravity
//#define D_GRAVITY -0.001982f

#define D_GRAVITY -9.82f
#define D_PI 3.1415926559f;
#define D_EPSILON 0.000000000000001f;

//since force is low a higher radius is requiered for small number of particles
#define D_RAD 0.05f;

SPH::SPH(glm::vec3* start_pos)
{
	m_particles.dens = new float[D_NR_OF_PARTICLES];
	m_particles.p = new float[D_NR_OF_PARTICLES];
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
	m_particles.mass = 0.00006f;
	float pi = D_PI;
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i) {
		m_particles.dens[i] = m_particles.mass / (4.f / 3.f * pow(m_particles.rad, 3)*pi);
		m_particles.p[i] = 0.f;
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

	init_positions(start_pos, 5, 5);
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
	delete[] m_particles.p;
	delete[] m_neighbor_data;
}

void SPH::update(float dT)
{
	static float alpha[D_NR_OF_PARTICLES];
	static float k_v_i[D_NR_OF_PARTICLES];
	static Float3s f_tot[D_NR_OF_PARTICLES];
	static float scalar_values[D_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS];
	static float kernel_values[D_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS];
	static float pi = D_PI;

	find_neighborhoods();

	//calculate_time_step(dT);

	update_scalar_function(&m_particles.pos, m_neighbor_data, scalar_values);

	update_kernel_values(kernel_values, &m_particles.pos, m_neighbor_data);

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, m_neighbor_data, alpha, kernel_values);

	pressure_forces();

	calculate_time_step(dT);

	non_pressure_forces();

	predict_velocities();

	correct_density_error(alpha, dT, scalar_values, f_tot, k_v_i);

	//update_positions(dT);

	find_neighborhoods();  // t + dt
	
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i) {
		m_particles.dens[i] = m_particles.mass / (4.f / 3.f * pow(m_particles.rad, 3)*pi);
	}

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, m_neighbor_data, alpha, kernel_values);


	// TODO: this function is incorretly, implemented correct
	correct_divergence_error(k_v_i, scalar_values, alpha);

	update_velocities();

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i) {
		m_particles.dens[i] = m_particles.mass /(4.f / 3.f * pow(m_particles.rad, 3)*pi );
		
	}

	update_positions(dT);

}


void SPH::init_positions(glm::vec3* start_pos, int rows, int cols) const
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

		m_particles.pos.y[i] = start_pos->x + x * dist_between*padding_factor;
		m_particles.pos.x[i] = start_pos->y + y * dist_between*padding_factor;
		m_particles.pos.z[i] = start_pos->z + z * dist_between*padding_factor;
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
				dist_i_n_2 = vector_i_n[0] * vector_i_n[0] + vector_i_n[1] * vector_i_n[1] + vector_i_n[2] * vector_i_n[2];
				
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
		m_particles.F_adv.y[i] = m_particles.mass * D_GRAVITY;
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
	m_delta_t = 0.4f*(m_particles.rad * 2) / sqrtf(v_max_2) - D_EPSILON;

	if (dT < m_delta_t)
		m_delta_t = dT;

	//m_delta_t = m_delta_t < 0.033f ? 0.033f : m_delta_t;

}

void SPH::predict_velocities()
{
	float mass = m_particles.mass;

	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)	
	{
		m_particles.pred_vel.x[i] = m_particles.vel.x[i]  + m_particles.F_adv.x[i] * m_delta_t/mass;
		m_particles.pred_vel.y[i] = m_particles.vel.y[i]  + m_particles.F_adv.y[i] * m_delta_t/mass;
		m_particles.pred_vel.z[i] = m_particles.vel.z[i]  + m_particles.F_adv.z[i] * m_delta_t/mass;

		if (abs(m_particles.pos.x[i] + m_particles.pred_vel.x[i] * m_delta_t) >= 0.05)
		{
			m_particles.pred_vel.x[i] = -0.1f*m_particles.pred_vel.x[i];
			m_particles.vel.x[i] = -0.1f*m_particles.vel.x[i];
		}

		if (abs(m_particles.pos.y[i] + m_particles.pred_vel.y[i] * m_delta_t) >= 0.2)
		{
			m_particles.pred_vel.y[i] = -0.0f*m_particles.pred_vel.y[i];
			m_particles.vel.y[i] = -0.0f*m_particles.vel.y[i];
		}

		if (abs(m_particles.pos.z[i] + m_particles.pred_vel.z[i] * m_delta_t) >= 0.05)
		{
			m_particles.pred_vel.z[i] = -0.1f*m_particles.pred_vel.z[i];
			m_particles.vel.z[i] = -0.1f*m_particles.vel.z[i];
		}
	}
}

void SPH::correct_density_error(float* alpha, float dT, float* scalar_values, Float3s* f_tot, float* k_v_i)
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
	float dx;
	float dy;
	float dz;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	int iter = -1;
	float p0{ 0.f };
	float p_avg{ 0.f };

	

	
		//TODO: why is avrg not used is it old code that should be removed or is this function wrong?
		float avrg = calculate_kv(alpha, &m_particles.vel, &m_particles.pred_vel, &m_particles.pos, m_particles.dens, m_delta_t, k_v_i, m_neighbor_data, scalar_values, m_particles.mass);

	float dens_min = 11111110.f, dens_max = 0.f;
	do{
		//p_avg = 0.f;
		calculate_pressure_force(f_tot, k_v_i, &m_particles.pos, m_particles.mass, scalar_values, m_neighbor_data, m_particles.dens);

		calculate_predicted_pressure(predicted_pressure, f_tot, m_particles.mass, m_particles.dens, scalar_values, m_delta_t, m_neighbor_data, &m_particles.pos);
		++iter;
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
			sum_x = sum_z = sum_y = .0f;

			assert(dens_i != 0, "dens");
			p_avg += dens_i;
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

				sum_x += m_particles.mass * (k_i_x + (k[neighbor_index].x / m_particles.dens[neighbor_index])) * kernel_gradient_x;
				sum_y += m_particles.mass * (k_i_y + (k[neighbor_index].y / m_particles.dens[neighbor_index])) * kernel_gradient_y;
				sum_z += m_particles.mass * (k_i_z + (k[neighbor_index].z / m_particles.dens[neighbor_index])) * kernel_gradient_z;

			}

			m_particles.pred_vel.x[i] -= m_delta_t * sum_x;
			m_particles.pred_vel.y[i] -= m_delta_t * sum_y;
			m_particles.pred_vel.z[i] -= m_delta_t * sum_z;
		}
		if (((p_avg / D_NR_OF_PARTICLES - C_REST_DENS) > (1.01*(dens_max - dens_min))))p_avg = 0.f;
		// condition p_avg - p_0 > ny  ny = 1.01*(p_max_p_min) p_0 = 1000
	} while ( ( (p_avg/D_NR_OF_PARTICLES - C_REST_DENS) > (1.01*(dens_max - dens_min) ) )  );
}

void SPH::correct_strain_rate_error() {}

void SPH::update_positions(float dT) const
{
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.pos.x[i] += m_particles.pred_vel.x[i] * m_delta_t;
		m_particles.pos.y[i] += m_particles.pred_vel.y[i] * m_delta_t;
		m_particles.pos.z[i] += m_particles.pred_vel.z[i] * m_delta_t;
	}
}

void SPH::correct_divergence_error(float* k_v_i, float* scalar_values, float* alpha)
{
	float dens_i;
	float div_i;
	float sum_x{ 0.f };
	float sum_y{ 0.f };
	float sum_z{ 0.f };

	// kernel variables
	float dx;
	float dy;
	float dz;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;

	int neighbor_index;

	float avrg;

	do
	{
		avrg = calculate_kv(alpha, &m_particles.vel, &m_particles.pred_vel, &m_particles.pos, m_particles.dens, m_delta_t, k_v_i, m_neighbor_data, scalar_values, m_particles.mass);
		
		for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
		{
			dens_i = m_particles.dens[i];

			assert(dens_i != 0.0f, "dens");
			div_i = k_v_i[i] / dens_i;

			for (auto j = 0; j < m_neighbor_data[i].n; ++j)
			{
				int linear_ind = j + D_MAX_NR_OF_NEIGHBORS*i;
				neighbor_index = m_neighbor_data[i].neighbor[j];

				assert(m_particles.dens[neighbor_index] != 0.0f, "n dens");

				dx = m_particles.pos.x[neighbor_index] - m_particles.pos.x[i];
				dy = m_particles.pos.y[neighbor_index] - m_particles.pos.y[i];
				dz = m_particles.pos.z[neighbor_index] - m_particles.pos.z[i];

				kernel_gradient_x = dx * scalar_values[linear_ind];
				kernel_gradient_y = dy * scalar_values[linear_ind];
				kernel_gradient_z = dz * scalar_values[linear_ind];

				sum_x += m_particles.mass * (div_i + k_v_i[neighbor_index] /
					m_particles.dens[neighbor_index]) * kernel_gradient_x;

				sum_y += m_particles.mass * (div_i + k_v_i[neighbor_index] /
					m_particles.dens[neighbor_index]) * kernel_gradient_y;

				sum_z += m_particles.mass * (div_i + k_v_i[neighbor_index] /
					m_particles.dens[neighbor_index]) * kernel_gradient_z;
			}

			m_particles.pred_vel.x[i] -= m_delta_t * sum_x;
			m_particles.pred_vel.y[i] -= m_delta_t * sum_y;
			m_particles.pred_vel.z[i] -= m_delta_t * sum_z;
		
			sum_x = sum_y = sum_z = .0f;
		}
	} while ((avrg) > 0.1f );
}

void SPH::update_velocities()
{
	for (auto i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.vel.x[i] = m_particles.pred_vel.x[i];
		m_particles.vel.y[i] = m_particles.pred_vel.y[i];
		m_particles.vel.z[i] = m_particles.pred_vel.z[i];

		if (abs(m_particles.pos.x[i] + m_particles.vel.x[i] * m_delta_t) >= 0.05) {
			m_particles.vel.x[i] = -0.1f*m_particles.vel.x[i];
			m_particles.pred_vel.x[i] = -0.1f*m_particles.pred_vel.x[i];
		}
		if (abs(m_particles.pos.y[i] + m_particles.vel.y[i] * m_delta_t) >= 0.20) {
			m_particles.vel.y[i] = 0.0f;
			m_particles.pred_vel.y[i] = 0.0f;
		}
		if (abs(m_particles.pos.z[i] + m_particles.vel.z[i] * m_delta_t) >= 0.05) {
			m_particles.vel.z[i] = -0.1f*m_particles.vel.z[i];
			m_particles.pred_vel.z[i] = -0.1f*m_particles.pred_vel.z[i];
		}
	}
}

void update_density_and_factors(float mass, Float3* pos, float* dens, float* scalar_values,
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
	float x = 0.f, y = 0.f, z = 0.f;
	const float min_denom{ 0.000001f };
	neighbor_mass = mass;
	for (auto particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		nr_neighbors = neighbor_data[particle].n;
		//if (nr_neighbors == 0)dens[particle] = mass / (4.f / 3.f * pi * pow(0.01, 3));

		for (auto neighbor = 0; neighbor < nr_neighbors; ++neighbor)
		{
			neighbor_index = neighbor_data[particle].neighbor[neighbor];

			int linear_ind = neighbor + D_MAX_NR_OF_NEIGHBORS*particle;

			//Update density
			dens[particle] += scalar_values[linear_ind] * neighbor_mass;// / (4.f / 3.f * pi * pow(0.01, 3));

			dx = pos->x[neighbor_index] - pos->x[particle];
			dy = pos->y[neighbor_index] - pos->y[particle];
			dz = pos->z[neighbor_index] - pos->z[particle];

			kernel_gradient_x = neighbor_mass * dx * scalar_values[linear_ind];
			kernel_gradient_y = neighbor_mass * dy * scalar_values[linear_ind];
			kernel_gradient_z = neighbor_mass * dz * scalar_values[linear_ind];

			x += kernel_gradient_x;
			y += kernel_gradient_y;
			z += kernel_gradient_z;

			sqrt_val = (kernel_gradient_x*kernel_gradient_x + kernel_gradient_y*kernel_gradient_y
				+ kernel_gradient_z*kernel_gradient_z);

			float temp = sqrt(sqrt_val);

			sum_abs_denom += temp*temp;
		}

		abs_sum_denom = sqrt(x*x + y*y + z*z);
		denom = (abs_sum_denom*abs_sum_denom + sum_abs_denom);

		// set alpha to max(denom,min_denom)
		denom = denom > min_denom ? denom : min_denom;
		alpha[particle] = -dens[particle] / denom;

		sum_abs_denom = 0.f;
		x = y = z = 0.f;
	}
}

void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data)
{
	unsigned int ind;
	float x, y, z, q;
	float kernel_val = 0;
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
			len = len == 0.f ? 0.0000000000000000001f : len;
			q = sqrt(len) / D_RAD;

			if (q >= 0 || q <= 0.5f)
				kernel_val *= (1 - 6.f*q*q + 6.f*q*q*q);
			else if (q > 0.5f || q <= 1.0f)
				kernel_val *= 2.0f*(1 - q)*(1 - q)*(1 - q);
			else
				kernel_val = 0;

			kernel_values[particle*D_NR_OF_PARTICLES + neighbor] = kernel_val;
		}
	}
}

void calculate_pressure_force(Float3s* f_tot, float* k_v_i, Float3* pos, float mass, float* scalar_values, Neighbor_Data* neighbor_data, float* dens)
{
	unsigned int neighbor_index;
	unsigned int n_neighbors;
	float sum_x{ 0.f };
	float sum_y{ 0.f };
	float sum_z{ 0.f };
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		n_neighbors = neighbor_data[i].n;

		float kv_dens = k_v_i[i] / dens[i];

		for (int j = 0; j < n_neighbors; ++j)
		{
			neighbor_index = neighbor_data[i].neighbor[j];
			int linear_ind = neighbor_index + D_MAX_NR_OF_NEIGHBORS*i;


			x = pos->x[neighbor_index] - pos->x[i];
			y = pos->y[neighbor_index] - pos->y[i];
			z = pos->z[neighbor_index] - pos->z[i];

			kernel_gradient_x = x * scalar_values[linear_ind];
			kernel_gradient_y = y * scalar_values[linear_ind];
			kernel_gradient_z = z * scalar_values[linear_ind];

			assert(dens[i] != 0.0f, "dens");
			sum_x += mass * (kv_dens + k_v_i[neighbor_index] / dens[neighbor_index])*kernel_gradient_x;
			sum_y += mass * (kv_dens + k_v_i[neighbor_index] / dens[neighbor_index])*kernel_gradient_y;
			sum_z += mass * (kv_dens + k_v_i[neighbor_index] / dens[neighbor_index])*kernel_gradient_z;
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
	//float d_t_2 = delta_t*delta_t;
	float res_x = 0.f, res_y = 0.f, res_z = 0.0f;
	unsigned int neighbor_index;
	float d_t_2 = delta_t*delta_t;
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
			int linear_ind = neighbor_index + D_MAX_NR_OF_NEIGHBORS*i;

			x = pos->x[neighbor_index] - pos->x[i];
			y = pos->y[neighbor_index] - pos->y[i];
			z = pos->z[neighbor_index] - pos->z[i];

			kernel_gradient_x = x * scalar_value[linear_ind];
			kernel_gradient_y = y * scalar_value[linear_ind];
			kernel_gradient_z = z * scalar_value[linear_ind];

			res_x += mass * (f_p[i].x / dens[i] - f_p[neighbor_index].x / dens[i])*kernel_gradient_x;
			res_y += mass * (f_p[i].y / dens[i] - f_p[neighbor_index].y / dens[i])*kernel_gradient_y;
			res_z += mass * (f_p[i].z / dens[i] - f_p[neighbor_index].z / dens[i])*kernel_gradient_z;
		}
		predicted_pressure[i].x = d_t_2 *res_x;
		predicted_pressure[i].y = d_t_2 *res_y;
		predicted_pressure[i].z = d_t_2 *res_z;
		res_x = res_y = res_z = 0.f;
	}
}

float calculate_kv(float* alpha, Float3* vel, Float3* pred_vel, Float3* pos, float* dens,
	float delta_t, float* k_v_i, Neighbor_Data* neighbor_data, float* scalar_value,float mass)
{
	//float x, y, z;
	float dx = 0.f, dy = 0.f, dz= 0.f;
	float x, y, z;
	unsigned int neighbor_index;
	float kernel_gradient_x;
	float kernel_gradient_y;
	float kernel_gradient_z;
	float d_dens{ 0.f };

	float d_dens_avrg = 0.0f;

	//not suppose to be nessesary bot in first iteratin time = 0
	//delta_t = delta_t <= 0 ? 0.00001f : delta_t;
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		int neighbor_length = neighbor_data[i].n;
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

			dx += mass*(pred_vel->x[i] - pred_vel->x[neighbor_index]) *kernel_gradient_x;
			dy += mass*(pred_vel->y[i] - pred_vel->y[neighbor_index]) *kernel_gradient_y;
			dz += mass*(pred_vel->z[i] - pred_vel->z[neighbor_index]) *kernel_gradient_z;

		}
		float particle_dens = dens[i];


		d_dens = particle_dens * -(dx + dy + dz);
		d_dens_avrg += d_dens;
		k_v_i[i] = 1.f*(1.f /(delta_t)	)* d_dens * alpha[i];
		dx = dy = dz = 0.f;
	}

	return d_dens_avrg / D_NR_OF_PARTICLES;
}

void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values)
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
			auto q = dist / D_RAD;
			float d = D_RAD;
			//this is an approximation 
			dist = dist == 0.f ? 0.000000000000000000000000000001f : dist;

			//Compute the derivitive of the kernel function
			if (q >= 0 || q <= 0.5f)
				kernel_derive = (-12.f*q + 18.f*q*q);
			else if (q > 0.5f || q <= 1.0f)
				kernel_derive = -6.0f*(1 - q)*(1 - q);
			else
				kernel_derive = 0.0f;

			scalar_values[particle * D_MAX_NR_OF_NEIGHBORS + neighbor] =( kernel_derive /  (dist / d));
		}
	}
}