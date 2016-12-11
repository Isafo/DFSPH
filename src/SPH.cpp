#include "SPH.h"
#include <math.h>
#include <assert.h>
#include <vector>
#include <array>

//Dan Koschier, https://github.com/InteractiveComputerGraphics/CompactNSearch
#include "CompactNSearch/include/CompactNSearch.h"
#include "CompactNSearch/include/DataStructures.h"
#include "imconfig.h"
#include "imgui.h"

#define D_PI 3.1415926559f;
#define D_EPSILON 10e-6f;
//since force is low a higher radius is requiered for small number of particles
#define D_SEARCH_RANGE 0.035f;


SPH::SPH()
{
	// parentheses initialize array to zero 
	m_particles.pos.x = new float[D_MAX_NR_OF_PARTICLES]();
	m_particles.pos.y = new float[D_MAX_NR_OF_PARTICLES]();
	m_particles.pos.z = new float[D_MAX_NR_OF_PARTICLES]();

	m_particles.vel.x = new float[D_MAX_NR_OF_PARTICLES]();
	m_particles.vel.y = new float[D_MAX_NR_OF_PARTICLES]();
	m_particles.vel.z = new float[D_MAX_NR_OF_PARTICLES]();

	m_particles.pred_vel.x = new float[D_MAX_NR_OF_PARTICLES]();
	m_particles.pred_vel.y = new float[D_MAX_NR_OF_PARTICLES]();
	m_particles.pred_vel.z = new float[D_MAX_NR_OF_PARTICLES]();

	m_particles.F_adv.x = new float[D_MAX_NR_OF_PARTICLES]();
	m_particles.F_adv.y = new float[D_MAX_NR_OF_PARTICLES]();
	m_particles.F_adv.z = new float[D_MAX_NR_OF_PARTICLES]();

	m_particles.dens = new float[D_MAX_NR_OF_PARTICLES]();

	m_neighbor_data = new Neighbor_Data[D_MAX_NR_OF_PARTICLES];

	m_alpha = new float[D_MAX_NR_OF_PARTICLES]();
	m_dens_derive = new float[D_MAX_NR_OF_PARTICLES]();
	m_pred_dens = new float[D_MAX_NR_OF_PARTICLES]();

	m_scalar_values = new float[D_MAX_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS]();
	m_kernel_values = new float[D_MAX_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS]();

	m_rad = 0.01f;
	m_mass = 0.004218f;

	m_wind.x = m_wind.y = m_wind.z = 0.0f;
	m_gravity = 9.82f;
}

void SPH::reset()
{
	memset(m_particles.pos.x, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.pos.x));
	memset(m_particles.pos.y, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.pos.x));
	memset(m_particles.pos.z, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.pos.x));

	memset(m_particles.vel.x, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.vel.x));
	memset(m_particles.vel.y, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.vel.x));
	memset(m_particles.vel.z, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.vel.x));

	memset(m_particles.pred_vel.x, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.pred_vel.x));
	memset(m_particles.pred_vel.y, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.pred_vel.x));
	memset(m_particles.pred_vel.z, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.pred_vel.x));

	memset(m_particles.dens, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_particles.dens));

	memset(m_alpha, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_alpha));
	memset(m_dens_derive, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_dens_derive));
	memset(m_pred_dens, 0, D_MAX_NR_OF_PARTICLES * sizeof(*m_pred_dens));

	memset(m_scalar_values, 0, D_MAX_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS * sizeof(*m_scalar_values));
	memset(m_kernel_values, 0, D_MAX_NR_OF_PARTICLES * D_MAX_NR_OF_NEIGHBORS * sizeof(*m_kernel_values));
}

SPH::~SPH()
{
	delete[] m_particles.pos.x;
	delete[] m_particles.pos.y;
	delete[] m_particles.pos.z;

	delete[] m_particles.vel.x;
	delete[] m_particles.vel.y;
	delete[] m_particles.vel.z;

	delete[] m_particles.pred_vel.x;
	delete[] m_particles.pred_vel.y;
	delete[] m_particles.pred_vel.z;

	delete[] m_particles.F_adv.x;
	delete[] m_particles.F_adv.y;
	delete[] m_particles.F_adv.z;

	delete[] m_particles.dens;
	delete[] m_neighbor_data;

	delete[] m_alpha;
	delete[] m_dens_derive;
	delete[] m_pred_dens;

	delete[] m_scalar_values;
	delete[] m_kernel_values;
}


void SPH::update()
{
	find_neighborhoods();

	update_scalar_function(&m_particles.pos, m_neighbor_data, m_scalar_values, m_simulated_particles);

	update_kernel_values(m_kernel_values, &m_particles.pos, m_neighbor_data, m_simulated_particles);

	update_density_and_factors(m_neighbor_data);

	calculate_time_step();

	non_pressure_forces();

	predict_velocities();

	correct_density_error();

	update_positions();

	find_neighborhoods();

	update_density_and_factors(m_neighbor_data);

	correct_divergence_error();

	update_velocities();
}

void SPH::init_positions(Float3s start_pos, int rows, int cols) const
{
	int ind;

	float x, y, z;
	float dist_between = 2.25f * m_rad;

	#pragma omp parallel
	#pragma omp for
	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		x = particle_ind / (rows * cols);
		ind = particle_ind - x * rows * cols;
		y = ind / rows;
		z = ind % rows;

		m_particles.pos.x[particle_ind] = start_pos.x + x * dist_between;
		m_particles.pos.y[particle_ind] = start_pos.y + y * dist_between;
		m_particles.pos.z[particle_ind] = start_pos.z + z * dist_between;
	}
}

void SPH::find_neighborhoods() const
{
	int count = 0;

	const float neigborhod_rad = D_SEARCH_RANGE;

	CompactNSearch::NeighborhoodSearch nsearch(neigborhod_rad);
	std::vector<std::array<double, 3>> positions(m_simulated_particles);

	#pragma omp parallel
	#pragma omp for
	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		positions.at(particle_ind) = { m_particles.pos.x[particle_ind], m_particles.pos.y[particle_ind], m_particles.pos.z[particle_ind] };
	}

	unsigned int point_set_id = nsearch.add_point_set(positions.front().data(), positions.size());

	nsearch.find_neighbors();

	CompactNSearch::PointSet const& ps = nsearch.point_set(point_set_id);
	
	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		for (int neighbor_ind = 0; neighbor_ind < ps.n_neighbors(particle_ind); ++neighbor_ind)
		{
			CompactNSearch::PointID const& pid = ps.neighbor(particle_ind, neighbor_ind);

			m_neighbor_data[particle_ind].neighbor[count] = pid.point_id;
			++count;
		}
		//save nr of neighbor to first position 
		m_neighbor_data[particle_ind].n = count;
		count = 0;
	}
}

void SPH::non_pressure_forces()
{
	Float3s drag;

	float Cd = 1.0f, rho = 1.0f;
	float pi = D_PI;
	float rad = m_rad;
	float sphere_area = 4.0f * pi * rad*rad;

	if(m_wind.x > 0.0f)
	{
		drag.x = rho * Cd * m_wind.x*m_wind.x * sphere_area;
	}
	else
	{
		drag.x = -rho * Cd * m_wind.x*m_wind.x * sphere_area;
	}

	if (m_wind.y > 0.0f)
	{
		drag.y = rho * Cd * m_wind.y*m_wind.y * sphere_area;
	}
	else
	{
		drag.y = -rho * Cd * m_wind.y*m_wind.y * sphere_area;
	}

	if (m_wind.z > 0.0f)
	{
		drag.z = rho * Cd * m_wind.z*m_wind.z * sphere_area;
	}
	else
	{
		drag.z = -rho * Cd * m_wind.z*m_wind.z * sphere_area;
	}

	#pragma omp parallel
	#pragma omp for
	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		m_particles.F_adv.x[particle_ind] = m_mass * drag.x;
		m_particles.F_adv.y[particle_ind] = m_mass * (drag.y - m_gravity);
		m_particles.F_adv.z[particle_ind] = m_mass * drag.z;
	}
}

void SPH::calculate_time_step()
{
	float v_max_2 = 0;
	float x_2, y_2, z_2;

	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		x_2 = m_particles.vel.x[particle_ind] * m_particles.vel.x[particle_ind];
		y_2 = m_particles.vel.y[particle_ind] * m_particles.vel.y[particle_ind];
		z_2 = m_particles.vel.z[particle_ind] * m_particles.vel.z[particle_ind];

		if (v_max_2 < x_2 + y_2 + z_2)
		{
			v_max_2 = x_2 + y_2 + z_2;
		}
	}

	m_delta_t = m_time_factor * (2.0f * m_rad) / (sqrtf(v_max_2) + 1e-6f);

	if (m_delta_t > 0.005f)
	{
		m_delta_t = 0.005f;
	}
	else if (m_delta_t < 0.0005f)
	{
		m_delta_t = 0.0005f;
	}

}


void SPH::predict_velocities()
{
	float dist_2, pos_x, pos_y, pos_z;

	#pragma omp parallel
	#pragma omp for
	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		dist_2 = (m_particles.pos.x[particle_ind] - sc.center.x) * (m_particles.pos.x[particle_ind] - sc.center.x) +
			(m_particles.pos.y[particle_ind] - sc.center.y) * (m_particles.pos.y[particle_ind] - sc.center.y) +
			(m_particles.pos.z[particle_ind] - sc.center.z) * (m_particles.pos.z[particle_ind] - sc.center.z);

		// Sphere
		if (dist_2 == sc.radius_2)
		{
			m_particles.pred_vel.x[particle_ind] = 4.f * m_particles.pos.x[particle_ind];
			m_particles.pred_vel.y[particle_ind] = 2.f * m_particles.pos.y[particle_ind];
			m_particles.pred_vel.z[particle_ind] = 4.f * m_particles.pos.z[particle_ind];
		}
		else if (dist_2 < sc.radius_2)
		{
			m_particles.pred_vel.x[particle_ind] = 4.f * m_particles.pos.x[particle_ind];
			m_particles.pred_vel.y[particle_ind] = 0.f * m_particles.pos.y[particle_ind];
			m_particles.pred_vel.z[particle_ind] = 4.f * m_particles.pos.z[particle_ind];
		}
		else
		{
			pos_x = m_particles.vel.x[particle_ind] + m_particles.F_adv.x[particle_ind] * m_delta_t / m_mass;
			pos_y = m_particles.vel.y[particle_ind] + m_particles.F_adv.y[particle_ind] * m_delta_t / m_mass;
			pos_z = m_particles.vel.z[particle_ind] + m_particles.F_adv.z[particle_ind] * m_delta_t / m_mass;

			if (pos_x < -2.0f)
			{
				m_particles.pred_vel.x[particle_ind] = -2.0f;
			}
			else if (pos_x > 2.0f)
			{
				m_particles.pred_vel.x[particle_ind] = 2.0f;
			}
			else
			{
				m_particles.pred_vel.x[particle_ind] = m_particles.vel.x[particle_ind] + m_particles.F_adv.x[particle_ind] * m_delta_t / m_mass;
			}

			if (pos_y < -2.0f)
			{
				m_particles.pred_vel.y[particle_ind] = -2.0f;
			}
			else if (pos_y > 2.0f)
			{
				m_particles.pred_vel.y[particle_ind] = 2.0f;
			}
			else
			{
				m_particles.pred_vel.y[particle_ind] = m_particles.vel.y[particle_ind] + m_particles.F_adv.y[particle_ind] * m_delta_t / m_mass;
			}

			if (pos_z < -2.0f)
			{
				m_particles.pred_vel.z[particle_ind] = -2.0f;
			}
			else if (pos_z > 2.0f)
			{
				m_particles.pred_vel.z[particle_ind] = 2.0f;
			}
			else
			{
				m_particles.pred_vel.z[particle_ind] = m_particles.vel.z[particle_ind] + m_particles.F_adv.z[particle_ind] * m_delta_t / m_mass;
			}
		}

		pos_x = m_particles.pos.x[particle_ind] + m_particles.pred_vel.x[particle_ind] * m_delta_t;
		pos_y = m_particles.pos.y[particle_ind] + m_particles.pred_vel.y[particle_ind] * m_delta_t;
		pos_z = m_particles.pos.z[particle_ind] + m_particles.pred_vel.z[particle_ind] * m_delta_t;

		// Box
		if (pos_x >= 0.5f && (m_wind.x > 0.0f || m_particles.pred_vel.x[particle_ind] > 0.0f))
		{
			m_particles.pred_vel.x[particle_ind] = 0.0f;
		}
		else if (pos_x <= -0.5f && (m_wind.x < 0.0f || m_particles.pred_vel.x[particle_ind] < 0.0f))
		{
			m_particles.pred_vel.x[particle_ind] = 0.0f;
		}

		if (pos_y <= -0.5f && (m_wind.y < 0.0f || m_particles.pred_vel.y[particle_ind] < 0.0f))
		{
			m_particles.pred_vel.y[particle_ind] = 0.0f;
		}

		if (pos_z >= 0.5f && (m_wind.z > 0.0f || m_particles.pred_vel.z[particle_ind] > 0.0f))
		{
			m_particles.pred_vel.z[particle_ind] = 0.0f;

		}
		else if (pos_z <= -0.5f && (m_wind.z < 0.0f || m_particles.pred_vel.z[particle_ind] < 0.0f))
		{
			m_particles.pred_vel.z[particle_ind] = 0.0f;
		}
	}
}

void SPH::correct_density_error()
{
	int neighbor_ind, linear_ind;
	int iter = 0;

	float eta;
	float scalar_value;
	float x, y, z;
	float kappa_i, kappa_j, div_i, div_j, div_sum;
	float pressure_acc_x, pressure_acc_y, pressure_acc_z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float inv_delta_t_2 = 1.0f / (m_delta_t * m_delta_t);

	calculate_derived_density_pred_dens(m_neighbor_data);

	do
	{
		for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
		{
			kappa_i = fmax(inv_delta_t_2 * (m_pred_dens[particle_ind] - C_REST_DENS) * m_alpha[particle_ind], 0.25f);

			div_i = kappa_i / m_particles.dens[particle_ind];

			pressure_acc_x = pressure_acc_y = pressure_acc_z = 0.0f;

			for (int neighbor = 0; neighbor < m_neighbor_data[particle_ind].n; ++neighbor)
			{
				linear_ind = neighbor + D_MAX_NR_OF_NEIGHBORS * particle_ind;
				neighbor_ind = m_neighbor_data[particle_ind].neighbor[neighbor];

				kappa_j = fmax(inv_delta_t_2 * (m_pred_dens[neighbor_ind] - C_REST_DENS) * m_alpha[neighbor_ind], 0.25f);

				x = m_particles.pos.x[particle_ind] - m_particles.pos.x[neighbor_ind];
				y = m_particles.pos.y[particle_ind] - m_particles.pos.y[neighbor_ind];
				z = m_particles.pos.z[particle_ind] - m_particles.pos.z[neighbor_ind];

				scalar_value = m_scalar_values[linear_ind];

				kernel_gradient_x = x * scalar_value;
				kernel_gradient_y = y * scalar_value;
				kernel_gradient_z = z * scalar_value;

				div_j = kappa_j / m_particles.dens[neighbor_ind];

				div_sum = div_i + div_j;

				pressure_acc_x += m_mass * div_sum * kernel_gradient_x;
				pressure_acc_y += m_mass * div_sum * kernel_gradient_y;
				pressure_acc_z += m_mass * div_sum * kernel_gradient_z;
			}
			m_particles.pred_vel.x[particle_ind] = m_particles.pred_vel.x[particle_ind] - m_delta_t * pressure_acc_x;
			m_particles.pred_vel.y[particle_ind] = m_particles.pred_vel.y[particle_ind] - m_delta_t * pressure_acc_y;
			m_particles.pred_vel.z[particle_ind] = m_particles.pred_vel.z[particle_ind] - m_delta_t * pressure_acc_z;
		}

		calculate_derived_density_pred_dens(m_neighbor_data);

		eta = 0.01f * m_density_error * C_REST_DENS;

		++iter;
	} while ((m_pred_dens_avg - C_REST_DENS > eta || iter < 2) && iter < m_Viter_max);
}

void SPH::update_positions() const
{
	#pragma omp parallel
	#pragma omp for
	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		m_particles.pos.x[particle_ind] += m_particles.pred_vel.x[particle_ind] * m_delta_t;
		m_particles.pos.y[particle_ind] += m_particles.pred_vel.y[particle_ind] * m_delta_t;
		m_particles.pos.z[particle_ind] += m_particles.pred_vel.z[particle_ind] * m_delta_t;
	}
}

/*
* ViscousDFSPH, Algorithm 2
*/
void SPH::correct_divergence_error()
{
	int neighbor_ind, linear_ind;
	int iter = 0;

	float eta;
	float scalar_value;
	float x, y, z;
	float kappa_v_i, kappa_v_j, div_i, div_j, div_sum;
	float pressure_acc_x, pressure_acc_y, pressure_acc_z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float inv_delta_t = 1.f / m_delta_t;

	calculate_derived_density_pred_dens(m_neighbor_data);

	do
	{
		for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
		{
			kappa_v_i = 0.5f * fmax(inv_delta_t * m_dens_derive[particle_ind] * m_alpha[particle_ind], 0.5f);

			div_i = kappa_v_i / m_particles.dens[particle_ind];

			pressure_acc_x = pressure_acc_y = pressure_acc_z = 0.0f;

			for (int neighbor = 0; neighbor < m_neighbor_data[particle_ind].n; ++neighbor)
			{
				linear_ind = neighbor + D_MAX_NR_OF_NEIGHBORS * particle_ind;
				neighbor_ind = m_neighbor_data[particle_ind].neighbor[neighbor];

				kappa_v_j = 0.5f * fmax(inv_delta_t * m_dens_derive[neighbor_ind] * m_alpha[neighbor_ind], 0.5f);

				x = m_particles.pos.x[particle_ind] - m_particles.pos.x[neighbor_ind];
				y = m_particles.pos.y[particle_ind] - m_particles.pos.y[neighbor_ind];
				z = m_particles.pos.z[particle_ind] - m_particles.pos.z[neighbor_ind];

				scalar_value = m_scalar_values[linear_ind];

				kernel_gradient_x = x * scalar_value;
				kernel_gradient_y = y * scalar_value;
				kernel_gradient_z = z * scalar_value;

				div_j = kappa_v_j / m_particles.dens[neighbor_ind];

				div_sum = div_i + div_j;

				pressure_acc_x += m_mass * div_sum * kernel_gradient_x;
				pressure_acc_y += m_mass * div_sum * kernel_gradient_y;
				pressure_acc_z += m_mass * div_sum * kernel_gradient_z;
			}
			m_particles.pred_vel.x[particle_ind] = m_particles.pred_vel.x[particle_ind] - m_delta_t * pressure_acc_x;
			m_particles.pred_vel.y[particle_ind] = m_particles.pred_vel.y[particle_ind] - m_delta_t * pressure_acc_y;
			m_particles.pred_vel.z[particle_ind] = m_particles.pred_vel.z[particle_ind] - m_delta_t * pressure_acc_z;
		}

		calculate_derived_density_pred_dens(m_neighbor_data);
		eta = 0.01f * m_divergence_error * C_REST_DENS * 1.0f / m_delta_t;
		++iter;

	} while (m_dens_derive_avg > eta && iter < m_iter_max); // implicit condition: iter < 1 
}

void SPH::update_velocities()
{
	#pragma omp parallel
	#pragma omp for
	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		m_particles.vel.x[particle_ind] = m_particles.pred_vel.x[particle_ind];
		m_particles.vel.y[particle_ind] = m_particles.pred_vel.y[particle_ind];
		m_particles.vel.z[particle_ind] = m_particles.pred_vel.z[particle_ind];
	}
}


void SPH::update_density_and_factors(Neighbor_Data* neighbor_data)
{
	int neighbor_ind, linear_ind;

	float denom;
	float abs_sum_denom, sum_abs_denom = 0;
	float x, y, z;
	float scalar_value_mul_mass;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z, kernel_gradient_abs;
	float kernel_gradient_x_sum = 0.0f, kernel_gradient_y_sum = 0.0f, kernel_gradient_z_sum = 0.0f;
	
	for (auto particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		//added 1 * mass as particles own density as kernel is 1 at dist == 0
		//this should atleast be the case, but needs to be checked
		m_particles.dens[particle_ind] = 12.000656593f;
		for (auto neighbor = 0; neighbor < neighbor_data[particle_ind].n; ++neighbor)
		{
			neighbor_ind = neighbor_data[particle_ind].neighbor[neighbor];

			linear_ind = neighbor + D_MAX_NR_OF_NEIGHBORS * particle_ind;

			//Update density
			m_particles.dens[particle_ind] += m_mass * m_kernel_values[linear_ind];

			x = m_particles.pos.x[particle_ind] - m_particles.pos.x[neighbor_ind];
			y = m_particles.pos.y[particle_ind] - m_particles.pos.y[neighbor_ind];
			z = m_particles.pos.z[particle_ind] - m_particles.pos.z[neighbor_ind];

			scalar_value_mul_mass = m_mass * m_scalar_values[linear_ind];

			kernel_gradient_x = x * scalar_value_mul_mass;
			kernel_gradient_y = y * scalar_value_mul_mass;
			kernel_gradient_z = z * scalar_value_mul_mass;

			kernel_gradient_x_sum += kernel_gradient_x;
			kernel_gradient_y_sum += kernel_gradient_y;
			kernel_gradient_z_sum += kernel_gradient_z;

			kernel_gradient_abs = sqrt(kernel_gradient_x * kernel_gradient_x + kernel_gradient_y * kernel_gradient_y + kernel_gradient_z * kernel_gradient_z);

			sum_abs_denom += kernel_gradient_abs * kernel_gradient_abs;
		}

		abs_sum_denom = sqrt(kernel_gradient_x_sum*kernel_gradient_x_sum + kernel_gradient_y_sum*kernel_gradient_y_sum + kernel_gradient_z_sum*kernel_gradient_z_sum);

		denom = abs_sum_denom * abs_sum_denom + sum_abs_denom;
		denom = fmax(denom, 1e-6f);

		m_alpha[particle_ind] = m_particles.dens[particle_ind] / denom;

		sum_abs_denom = 0.f;
		kernel_gradient_x_sum = kernel_gradient_y_sum = kernel_gradient_z_sum = 0.f;
	}	
}

/*
* Using a Cubic spline kernel
* Divergence-Free SPH for Incompressible and Viscous Fluids, ref 6
*/
void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data, const int N_PARTICLES)
{
	int neighbor_ind;

	float x, y, z;
	float h, h_2;
	float dist;
	float search_area = D_SEARCH_RANGE;
	float pi = D_PI;
	float div = 1.0f / (search_area * search_area * search_area * pi);

	for (int particle_ind = 0; particle_ind < N_PARTICLES; ++particle_ind)
	{
		for (int neighbor = 0; neighbor < neighbor_data[particle_ind].n; ++neighbor)
		{
			neighbor_ind = neighbor_data[particle_ind].neighbor[neighbor];

			x = pos->x[particle_ind] - pos->x[neighbor_ind];
			y = pos->y[particle_ind] - pos->y[neighbor_ind];
			z = pos->z[particle_ind] - pos->z[neighbor_ind];

			dist = sqrt(x*x + y*y + z*z);
			h = dist / D_SEARCH_RANGE;
			h_2 = h * h;

			// length is always equal or smaller to D_SEARCH_RANGE => implicit intervall between [0, 1]
			kernel_values[particle_ind * D_MAX_NR_OF_NEIGHBORS + neighbor] = div * (1.0f - 1.5f*h_2 + 0.75f*h_2*h);
		}
	}
}


/*
* ViscousDFSPH, eq 9
*/
void SPH::calculate_derived_density_pred_dens(Neighbor_Data* neighbor_data)
{
	int neighbor_ind, linear_ind;

	float pressure_derived;
	float dens_derive_sum = 0.0f, pred_dens_sum = 0.0f;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float pressure_derived_x = 0.0f, pressure_derived_y = 0.0f, pressure_derived_z = 0.0f;

	for (int particle_ind = 0; particle_ind < m_simulated_particles; ++particle_ind)
	{
		for (int neighbor = 0; neighbor < neighbor_data[particle_ind].n; ++neighbor)
		{
			neighbor_ind = neighbor_data[particle_ind].neighbor[neighbor];
			linear_ind = neighbor + D_MAX_NR_OF_NEIGHBORS * particle_ind;

			x = m_particles.pos.x[particle_ind] - m_particles.pos.x[neighbor_ind];
			y = m_particles.pos.y[particle_ind] - m_particles.pos.y[neighbor_ind];
			z = m_particles.pos.z[particle_ind] - m_particles.pos.z[neighbor_ind];

			kernel_gradient_x = x * m_scalar_values[linear_ind];
			kernel_gradient_y = y * m_scalar_values[linear_ind];
			kernel_gradient_z = z * m_scalar_values[linear_ind];

			//equation 9 in DFSPH, changed 16-11-18
			pressure_derived_x += m_mass * kernel_gradient_x * (m_particles.pred_vel.x[particle_ind] - m_particles.pred_vel.x[neighbor_ind]);
			pressure_derived_y += m_mass * kernel_gradient_y * (m_particles.pred_vel.y[particle_ind] - m_particles.pred_vel.y[neighbor_ind]);
			pressure_derived_z += m_mass * kernel_gradient_z * (m_particles.pred_vel.z[particle_ind] - m_particles.pred_vel.z[neighbor_ind]);
		}

		pressure_derived = pressure_derived_x + pressure_derived_y + pressure_derived_z;

		m_dens_derive[particle_ind] = pressure_derived;

		m_pred_dens[particle_ind] = m_particles.dens[particle_ind] + m_delta_t * m_dens_derive[particle_ind];

		dens_derive_sum += m_dens_derive[particle_ind];
		pred_dens_sum += m_pred_dens[particle_ind];

		pressure_derived_x = pressure_derived_y = pressure_derived_z = 0.0f;
	}
	m_dens_derive_avg = dens_derive_sum / m_simulated_particles;
	m_pred_dens_avg = pred_dens_sum / m_simulated_particles;
}



/*
* Derived the Cubic spline kernel by q
* Divergence-Free SPH for Incompressible and Viscous Fluids, section 4.2
*/
void update_scalar_function(Float3* pos, Neighbor_Data* neighbor_data, float* scalar_values, const int N_PARTICLES)
{
	int neighbor_ind;

	float kernel_derive, scalar_value;
	float search_area = D_SEARCH_RANGE;
	float pi = D_PI;
	float h, dist;
	float x, y, z;
	float inv_range = 1.0f / D_SEARCH_RANGE;
	float div = 1.0f / (search_area * search_area * search_area * pi);

	for (int particle_ind = 0; particle_ind < N_PARTICLES; ++particle_ind)
	{
		for (int neighbor = 0; neighbor < neighbor_data[particle_ind].n; ++neighbor)
		{
			neighbor_ind = neighbor_data[particle_ind].neighbor[neighbor];

			x = pos->x[particle_ind] - pos->x[neighbor_ind];
			y = pos->y[particle_ind] - pos->y[neighbor_ind];
			z = pos->z[particle_ind] - pos->z[neighbor_ind];

			dist = sqrt(x*x + y*y + z*z);

			h = dist * inv_range;

			// length is always equal or smaller to D_SEARCH_RANGE => implicit intervall between [0, 1]
			kernel_derive = div * (-3.0f*h + 2.25f*h*h);

			scalar_value = kernel_derive / (search_area * dist);

			scalar_values[particle_ind * D_MAX_NR_OF_NEIGHBORS + neighbor] = scalar_value;
		}
	}
}
