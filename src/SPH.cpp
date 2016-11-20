#include "SPH.h"
#include <iostream>
#include <math.h>
#include <assert.h>

#define D_GRAVITY -9.82f
#define D_PI 3.1415926559f;
#define D_EPSILON 0.000000000000001f;

//since force is low a higher radius is requiered for small number of particles
#define D_SEARCH_RANGE 0.05f;

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
	m_particles.mass = 0.00418f;

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

	init_positions(x, y, z, 8, 8);
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
	static float dens_derive[D_NR_OF_PARTICLES];
	static float scalar_values[D_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS];
	static float kernel_values[D_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS];
	
	find_neighborhoods();

	update_scalar_function(&m_particles.pos, m_neighbor_data, scalar_values);

	update_kernel_values(kernel_values, &m_particles.pos, m_neighbor_data);

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, m_neighbor_data, alpha, kernel_values);

	calculate_time_step(dT);

	non_pressure_forces();

	predict_velocities();

	//correct_density_error(alpha, scalar_values, f_tot, k_v_i);

	update_positions();

	find_neighborhoods();

	update_density_and_factors(m_particles.mass, &m_particles.pos, m_particles.dens, scalar_values, m_neighbor_data, alpha, kernel_values);

	correct_divergence_error(dens_derive, scalar_values, alpha);

	update_velocities();
}

void SPH::init_positions(int x_start, int y_start, int z_start, int rows, int cols) const
{
	int ind;

	float dist_between{ 1.2f * m_particles.rad };
	float padding_factor{ 1.8f };
	float x, y, z;

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
	const float neigborhod_rad = D_SEARCH_RANGE;
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
		if (v_max_2 < x_2 + y_2 + z_2)
			v_max_2 = x_2 + y_2 + z_2;
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

		if (abs(m_particles.pos.x[i] + m_particles.pred_vel.x[i] * m_delta_t) >= 0.2f)
		{
			m_particles.pred_vel.x[i] = 0.f;
		}

		if (m_particles.pos.y[i] + m_particles.pred_vel.y[i] * m_delta_t <= -0.3f)
		{
			m_particles.pred_vel.y[i] = 0.f;
		} 
		else if(abs(m_particles.pos.y[i] + m_particles.pred_vel.y[i] * m_delta_t) >= 0.27f) // Ceil
		{
			m_particles.pred_vel.y[i] = m_particles.F_adv.y[i] * m_delta_t / m_particles.mass;
		}

		if (abs(m_particles.pos.z[i] + m_particles.pred_vel.z[i] * m_delta_t) >= 0.2f)
		{
			m_particles.pred_vel.z[i] = 0.f;
		}
	}
}


// TODO rewrite 
/*
void SPH::correct_density_error(float* alpha, float* scalar_values, Float3s* f_tot, float* k_v_i)
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
		//calculate_predicted_pressure(predicted_pressure, f_tot, m_particles.mass, m_particles.dens, scalar_values, m_delta_t, m_neighbor_data, &m_particles.pos);
		
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
*/

//void SPH::correct_strain_rate_error() {}

void SPH::update_positions() const
{
	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		m_particles.pos.x[i] += m_particles.pred_vel.x[i] * m_delta_t;
		m_particles.pos.y[i] += m_particles.pred_vel.y[i] * m_delta_t;
		m_particles.pos.z[i] += m_particles.pred_vel.z[i] * m_delta_t;
	}
}

/*
 * ViscousDFSPH, Algorithm 2
 */
void SPH::correct_divergence_error(float* dens_derive, float* scalar_values, float* alpha)
{
	int neighbor_ind;
	int iter = 0;

	float dens_derive_avg;
	float k_v_i, k_v_j, div_i, div_j;
	float pressure_acc_x, pressure_acc_y, pressure_acc_z;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	
	do
	{
		dens_derive_avg = calculate_derived_density(dens_derive, &m_particles.pred_vel, m_particles.mass, scalar_values, m_neighbor_data, &m_particles.pos);
		
		for (auto particle_ind = 0; particle_ind < D_NR_OF_PARTICLES; ++particle_ind)
		{

			k_v_i = (1.f / m_delta_t)*dens_derive[particle_ind]*alpha[particle_ind];

			pressure_acc_x = pressure_acc_y = pressure_acc_z = 0.0f;
			
			for (auto j = 0; j < m_neighbor_data[particle_ind].n; ++j)
			{
				int linear_ind = j + D_MAX_NR_OF_NEIGHBORS*particle_ind;
				neighbor_ind = m_neighbor_data[particle_ind].neighbor[j];

				assert(m_particles.dens[neighbor_ind] != 0.0f, "n dens");

				k_v_j = (1.f / m_delta_t)*dens_derive[neighbor_ind] * alpha[neighbor_ind];

				x = m_particles.pos.x[particle_ind] - m_particles.pos.x[neighbor_ind];
				y = m_particles.pos.y[particle_ind] - m_particles.pos.y[neighbor_ind];
				z = m_particles.pos.z[particle_ind] - m_particles.pos.z[neighbor_ind];

				kernel_gradient_x = x*scalar_values[linear_ind];
				kernel_gradient_y = y*scalar_values[linear_ind];
				kernel_gradient_z = z*scalar_values[linear_ind];

				div_i = k_v_i / m_particles.dens[particle_ind];
				div_j = k_v_j / m_particles.dens[neighbor_ind];

				pressure_acc_x += m_particles.mass * (div_i + div_j) * kernel_gradient_x;
				pressure_acc_y += m_particles.mass * (div_i + div_j) * kernel_gradient_y;
				pressure_acc_z += m_particles.mass * (div_i + div_j) * kernel_gradient_z;
			}
			//pressure_force_z is not in report but it is a force and it is = F/m *delta_t
			m_particles.pred_vel.x[particle_ind] = m_particles.pred_vel.x[particle_ind] - m_delta_t * pressure_acc_x;
			m_particles.pred_vel.y[particle_ind] = m_particles.pred_vel.y[particle_ind] - m_delta_t * pressure_acc_y;
			m_particles.pred_vel.z[particle_ind] = m_particles.pred_vel.z[particle_ind] - m_delta_t * pressure_acc_z;
		}
	++iter;
	} while (dens_derive_avg > 1.f || iter < 1);

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
	int nr_neighbors, neighbor_index, linear_ind;

	float abs_sum_denom, sum_abs_denom = 0;
	float denom;
	float dx, dy, dz;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float sqrt_val;
	float x = 0.f, y = 0.f, z = 0.f;

	const float min_denom{ 0.000001f };

	for (auto particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		nr_neighbors = neighbor_data[particle].n;
		//added 1 * mass as particles own density as kernel is 1 at dist == 0
		//this should atleast be the case, but needs to be checked
		// => Does not seem to cause a problem when it is 0. So i followed the report
		dens[particle] = 0.f;
		for (auto neighbor = 0; neighbor < nr_neighbors; ++neighbor)
		{
			neighbor_index = neighbor_data[particle].neighbor[neighbor];

			linear_ind = neighbor + D_MAX_NR_OF_NEIGHBORS*particle;

			//Update density
			dens[particle] += mass*kernel_values[linear_ind];

			dx = pos->x[particle] - pos->x[neighbor_index];
			dy = pos->y[particle] - pos->y[neighbor_index];
			dz = pos->z[particle] - pos->z[neighbor_index];

			kernel_gradient_x = mass*dx*scalar_values[linear_ind];
			kernel_gradient_y = mass*dy*scalar_values[linear_ind];
			kernel_gradient_z = mass*dz*scalar_values[linear_ind];

			x += kernel_gradient_x;
			y += kernel_gradient_y;
			z += kernel_gradient_z;

			sqrt_val = sqrt(kernel_gradient_x*kernel_gradient_x + kernel_gradient_y*kernel_gradient_y
				+ kernel_gradient_z*kernel_gradient_z);

			sum_abs_denom += sqrt_val*sqrt_val;;
		}

		abs_sum_denom = sqrt(x*x + y*y + z*z);
		denom = abs_sum_denom*abs_sum_denom + sum_abs_denom;

		// set alpha to max(denom,min_denom)
		denom = denom < min_denom ? min_denom : denom;
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

	float x, y, z;
	float q, q_2; 
	float length;
	float kernel_val;
	float search_area = D_SEARCH_RANGE;
	float pi = D_PI;
	float particle_pos_x, particle_pos_y, particle_pos_z;

	for (auto particle = 0; particle < D_NR_OF_PARTICLES; ++particle)
	{
		particle_pos_x = pos->x[particle];
		particle_pos_y = pos->y[particle];
		particle_pos_z = pos->z[particle];

		for (auto neighbor = 0; neighbor < neighbor_data[particle].n; ++neighbor)
		{
			ind = neighbor_data[particle].neighbor[neighbor];

			x = pos->x[ind] - particle_pos_x;
			y = pos->y[ind] - particle_pos_y;
			z = pos->z[ind] - particle_pos_z;

			length = sqrt(x*x + y*y + z*z);
			q = length / D_SEARCH_RANGE;
			q_2 = q*q;

			// length is always equal or smaller to D_SEARCH_RANGE => implicit intervall between [0, 1]
			kernel_val = (1.0f/(search_area*pi))*(1.0f - 1.5f*q_2 + 0.75f*q_2*q);

			kernel_values[particle*D_NR_OF_PARTICLES + neighbor] = kernel_val;
		}
	}
}

/*
 * ViscousDFSPH, eq 9
 */
float calculate_derived_density(float* derived_density, Float3* pred_vel, float mass, float* scalar_value, Neighbor_Data* neighbor_data, Float3* pos)
{
	int neighbor_index, linear_ind;
	int neighbor_length;

	float pressure_derived, dens_derive_sum = 0.f;
	float pressure_derived_x = 0.f, pressure_derived_y = 0.f, pressure_derived_z = 0.f;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;

	for (int i = 0; i < D_NR_OF_PARTICLES; ++i)
	{
		neighbor_length = neighbor_data[i].n;
		for (int j = 0; j < neighbor_length; ++j)
		{
			neighbor_index = neighbor_data[i].neighbor[j];
			linear_ind = neighbor_index + D_MAX_NR_OF_NEIGHBORS*i;

			x = pos->x[i] - pos->x[neighbor_index];
			y = pos->y[i] - pos->y[neighbor_index];
			z = pos->z[i] - pos->z[neighbor_index];

			kernel_gradient_x = x*scalar_value[linear_ind];
			kernel_gradient_y = y*scalar_value[linear_ind];
			kernel_gradient_z = z*scalar_value[linear_ind];

			//equation 9 in DFSPH, changed 16-11-18
			pressure_derived_x += mass*kernel_gradient_x*(pred_vel->x[i] - pred_vel->x[neighbor_index]);
			pressure_derived_y += mass*kernel_gradient_y*(pred_vel->y[i] - pred_vel->y[neighbor_index]);
			pressure_derived_z += mass*kernel_gradient_z*(pred_vel->z[i] - pred_vel->z[neighbor_index]);
		}

		pressure_derived = pressure_derived_x + pressure_derived_y + pressure_derived_z;

		derived_density[i] = pressure_derived;

		dens_derive_sum += derived_density[i];
		pressure_derived_x = pressure_derived_y = pressure_derived_z = 0.f;
	}
	return dens_derive_sum / D_NR_OF_PARTICLES;
}


/*
* Derived the Cubic spline kernel by q
* Divergence-Free SPH for Incompressible and Viscous Fluids, section 4.2
*/
void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values)
{
	int neighbor_ind;

	float kernel_derive, scalar_value;
	float search_area = D_SEARCH_RANGE;
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
			q = dist / D_SEARCH_RANGE;
			q_2 = q*q;

			// length is always equal or smaller to D_SEARCH_RANGE => implicit intervall between [0, 1]
			kernel_derive = (1.0f / (search_area*pi))*(-3.0f*q + 2.25f*q_2);

			scalar_value = kernel_derive*(1.0f /  (search_area*dist));

			scalar_values[particle*D_MAX_NR_OF_NEIGHBORS + neighbor] = scalar_value;

		}
	}

}