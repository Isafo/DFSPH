#include "SPH.h"
#include <math.h>
#include <assert.h>
#include <vector>
#include <array>

//Compact Search was written by
//Dan Koschier, https://github.com/InteractiveComputerGraphics/CompactNSearch
#include "CompactNSearch/include/CompactNSearch.h"
#include "CompactNSearch/include/DataStructures.h"
#include "imconfig.h"

#define D_GRAVITY -9.82f
#define D_PI 3.1415926559f;
#define D_EPSILON 10e-6f;
//since force is low a higher radius is requiered for small number of particles
#define D_SEARCH_RANGE 0.035f;

SPH::SPH(int n_particles) : C_N_PARTICLES{ n_particles } {}

void SPH::init()
{
	m_particles.dens = new float[C_N_PARTICLES];
	m_particles.pos.x = new float[C_N_PARTICLES];
	m_particles.pos.y = new float[C_N_PARTICLES];
	m_particles.pos.z = new float[C_N_PARTICLES];
	m_particles.vel.x = new float[C_N_PARTICLES];
	m_particles.vel.y = new float[C_N_PARTICLES];
	m_particles.vel.z = new float[C_N_PARTICLES];

	// predicted velocities
	m_particles.pred_vel.x = new float[C_N_PARTICLES];
	m_particles.pred_vel.y = new float[C_N_PARTICLES];
	m_particles.pred_vel.z = new float[C_N_PARTICLES];

	m_particles.F_adv.x = new float[C_N_PARTICLES];
	m_particles.F_adv.y = new float[C_N_PARTICLES];
	m_particles.F_adv.z = new float[C_N_PARTICLES];
	m_neighbor_data = new Neighbor_Data[C_N_PARTICLES];

	m_alpha = new float[C_N_PARTICLES];
	m_dens_derive = new float[C_N_PARTICLES];
	m_pred_dens = new float[C_N_PARTICLES];
	m_scalar_values = new float[C_N_PARTICLES * D_MAX_NR_OF_NEIGHBORS];
	m_kernel_values = new float[C_N_PARTICLES * D_MAX_NR_OF_NEIGHBORS];

	m_rad = 0.01f;
	m_mass = 0.004218f;

	for (auto i = 0; i < C_N_PARTICLES; ++i)
	{
		m_particles.pos.x[i] = 0.f;
		m_particles.pos.y[i] = 0.f;
		m_particles.pos.z[i] = 0.f;

		m_particles.vel.x[i] = 0.f;
		m_particles.vel.y[i] = 0.f;
		m_particles.vel.z[i] = 0.f;

		m_particles.pred_vel.x[i] = 0.f;
		m_particles.pred_vel.y[i] = 0.f;
		m_particles.pred_vel.z[i] = 0.f;

		m_particles.F_adv.x[i] = 0.f;
		m_particles.F_adv.y[i] = m_mass*D_GRAVITY;
		m_particles.F_adv.z[i] = 0.f;

		m_particles.dens[i] = 0.f;

		m_alpha[i] = 0.0f;
		m_dens_derive[i] = 0.0f;
		m_pred_dens[i] = 0.0f;
		m_scalar_values[i] = 0.0f;
		m_kernel_values[i] = 0.0f;
	}
}

void SPH::reset() {
	for (auto i = 0; i < C_N_PARTICLES; ++i)
	{
		m_particles.dens[i] = 0.f;
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

		m_alpha[i] = 0.0f;
		m_dens_derive[i] = 0.0f;
		m_pred_dens[i] = 0.0f;
		m_scalar_values[i] = 0.0f;
		m_kernel_values[i] = 0.0f;
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
	delete[] m_particles.pred_vel.x;
	delete[] m_particles.pred_vel.y;
	delete[] m_particles.pred_vel.z;
	delete[] m_particles.pos.x;
	delete[] m_particles.pos.y;
	delete[] m_particles.pos.z;
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

	update_scalar_function(&m_particles.pos, m_neighbor_data, m_scalar_values, current_n_particles);

	update_kernel_values(m_kernel_values, &m_particles.pos, m_neighbor_data, current_n_particles);

	update_density_and_factors(m_neighbor_data);

	calculate_time_step();

	//non_pressure_forces();

	predict_velocities();

	correct_density_error();

	update_positions();

	find_neighborhoods();

	update_density_and_factors(m_neighbor_data);

	correct_divergence_error();

	update_velocities();
}

void SPH::init_positions(int x_start, int y_start, int z_start, int rows, int cols) const
{
	int ind;

	float dist_between{ 2.25f * m_rad };
	//float padding_factor{ 1.8f };
	float x, y, z;

	#pragma omp parallel
	#pragma omp for
	for (auto i = 0; i < current_n_particles; ++i)
	{
		x = i / (rows*cols);
		ind = i - x*rows*cols;
		y = ind / rows;
		z = ind % rows;

		m_particles.pos.y[i] = x_start + x * dist_between;
		m_particles.pos.x[i] = y_start + y * dist_between;
		m_particles.pos.z[i] = z_start + z * dist_between;
	}
}

void SPH::find_neighborhoods() const
{
	const float neigborhod_rad = D_SEARCH_RANGE;
	int count{ 0 };

	CompactNSearch::NeighborhoodSearch nsearch(neigborhod_rad);
	std::vector<std::array<double, 3>> positions(current_n_particles);

	#pragma omp parallel
	#pragma omp for
	for (int i = 0; i < current_n_particles; ++i)
	{
		positions.at(i) = { m_particles.pos.x[i], m_particles.pos.y[i], m_particles.pos.z[i] };
	}

	unsigned int point_set_id = nsearch.add_point_set(positions.front().data(), positions.size());
	nsearch.find_neighbors();
	CompactNSearch::PointSet const& ps = nsearch.point_set(point_set_id);

	//#pragma omp parallel
	//#pragma omp for
	for (auto i = 0; i < current_n_particles; ++i)
	{
		for (auto n = 0; n < ps.n_neighbors(i); ++n)
		{
			CompactNSearch::PointID const& pid = ps.neighbor(i, n);

			m_neighbor_data[i].neighbor[count] = pid.point_id;
			++count;
		}
		//save nr of neighbor to first position 
		m_neighbor_data[i].n = count;
		count = 0;
	}
}

void SPH::non_pressure_forces() const
{
	#pragma omp parallel
	#pragma omp for
	for (auto i = 0; i < current_n_particles; ++i)
	{
		m_particles.F_adv.x[i] = 0.f;
		m_particles.F_adv.y[i] = m_mass*D_GRAVITY;
		m_particles.F_adv.z[i] = 0.f;
	}
}

void SPH::calculate_time_step()
{
	float v_max_2 = 0;
	float x_2, y_2, z_2;

	for (auto i = 0; i < current_n_particles; ++i)
	{
		x_2 = m_particles.vel.x[i] * m_particles.vel.x[i];
		y_2 = m_particles.vel.y[i] * m_particles.vel.y[i];
		z_2 = m_particles.vel.z[i] * m_particles.vel.z[i];
		if (v_max_2 < x_2 + y_2 + z_2)
			v_max_2 = x_2 + y_2 + z_2;
	}

	m_delta_t = time_factor * (2.f * m_rad) / (sqrtf(v_max_2) + 0.000001f);
	//m_delta_t = 0.003f;

	if (m_delta_t > 0.005f)
		m_delta_t = 0.005f;
	else if (m_delta_t < 0.0005f)
		m_delta_t = 0.0005f;

}

void SPH::predict_velocities()
{

	static float dist;
	#pragma omp parallel
	#pragma omp for
	for (auto i = 0; i < current_n_particles; ++i)
	{
		dist = ((m_particles.pos.x[i] - sc.center.x) * (m_particles.pos.x[i] - sc.center.x) +
			(m_particles.pos.y[i] - sc.center.y) * (m_particles.pos.y[i] - sc.center.y) +
			(m_particles.pos.z[i] - sc.center.z) * (m_particles.pos.z[i] - sc.center.z));

		if (dist == sc.radius_2 )
		{
			m_particles.pred_vel.x[i] = 4.f*(m_particles.pos.x[i]);// +m_particles.F_adv.x[i] * m_delta_t / m_mass;
			m_particles.pred_vel.y[i] = 2.f*(m_particles.pos.y[i]);// +m_particles.F_adv.y[i] * m_delta_t / m_mass;
			m_particles.pred_vel.z[i] = 4.f*(m_particles.pos.z[i]);// + m_particles.F_adv.z[i] * m_delta_t / m_mass;
		}
		else if(dist < sc.radius_2)
		{
			m_particles.pred_vel.x[i] = 4.f*(m_particles.pos.x[i]);// +m_particles.F_adv.x[i] * m_delta_t / m_mass;
			m_particles.pred_vel.y[i] = 0.f*(m_particles.pos.y[i]);// + m_particles.F_adv.y[i] * m_delta_t / m_mass;
			m_particles.pred_vel.z[i] = 4.f*(m_particles.pos.z[i]);// + m_particles.F_adv.z[i] * m_delta_t / m_mass;
		}
		else 
		{
			if(m_particles.vel.x[i] + m_particles.F_adv.x[i] * m_delta_t / m_mass < -2.0f)
				m_particles.pred_vel.x[i] = -2.0f;
			else if(m_particles.vel.x[i] + m_particles.F_adv.x[i] * m_delta_t / m_mass > 2.0f)
				m_particles.pred_vel.x[i] = 2.0f;
			else m_particles.pred_vel.x[i] = m_particles.vel.x[i] + m_particles.F_adv.x[i] * m_delta_t / m_mass;

			if (m_particles.vel.y[i] + m_particles.F_adv.y[i] * m_delta_t / m_mass < -2.0f)
				m_particles.pred_vel.y[i] = -2.0f;
			else if (m_particles.vel.y[i] + m_particles.F_adv.y[i] * m_delta_t / m_mass > 2.0f)
				m_particles.pred_vel.y[i] = 2.0f;
			else m_particles.pred_vel.y[i] = m_particles.vel.y[i] + m_particles.F_adv.y[i] * m_delta_t / m_mass;

			if (m_particles.vel.z[i] + m_particles.F_adv.z[i] * m_delta_t / m_mass < -2.0f)
				m_particles.pred_vel.z[i] = -2.0f;
			else if (m_particles.vel.z[i] + m_particles.F_adv.z[i] * m_delta_t / m_mass > 2.0f)
				m_particles.pred_vel.z[i] = 2.0f;
			else m_particles.pred_vel.z[i] = m_particles.vel.z[i] + m_particles.F_adv.z[i] * m_delta_t / m_mass;

		}

		
		if (abs(m_particles.pos.x[i] + m_particles.pred_vel.x[i] * m_delta_t) >= 0.5f)
		{
			m_particles.pred_vel.x[i] = 0.0f;
		
		}

		if (abs(m_particles.pos.y[i] + m_particles.pred_vel.y[i] * m_delta_t) >= 0.5f)
		{
			m_particles.pred_vel.y[i] = 0.0f;
		}

		if (abs(m_particles.pos.z[i] + m_particles.pred_vel.z[i] * m_delta_t) >= 0.5f)
		{
			m_particles.pred_vel.z[i] = 0.0f;
		
		}
	}
}

void SPH::correct_density_error()
{
	int neighbor_ind;
	int iter{ 0 };
	float eta;
	float k_i, k_j, div_i, div_j, div_sum;
	float pressure_acc_x, pressure_acc_y, pressure_acc_z;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float scalar_value;

	float inv_delta_t_2 = 1.f / (m_delta_t*m_delta_t);

	calculate_derived_density_pred_dens(m_neighbor_data);
	do
	{
		for (auto particle_ind = 0; particle_ind < current_n_particles; ++particle_ind)
		{

			k_i = fmax(inv_delta_t_2*(m_pred_dens[particle_ind] - C_REST_DENS) * m_alpha[particle_ind], 0.25f);
			div_i = k_i / m_particles.dens[particle_ind];

			pressure_acc_x = pressure_acc_y = pressure_acc_z = 0.0f;

			for (auto j = 0; j < m_neighbor_data[particle_ind].n; ++j)
			{
				int linear_ind = j + D_MAX_NR_OF_NEIGHBORS*particle_ind;
				neighbor_ind = m_neighbor_data[particle_ind].neighbor[j];

				//assert(m_particles.dens[neighbor_ind] != 0.0f, "n dens");

				k_j = fmax(inv_delta_t_2*(m_pred_dens[neighbor_ind] - C_REST_DENS) * m_alpha[neighbor_ind], 0.25f);

				x = m_particles.pos.x[particle_ind] - m_particles.pos.x[neighbor_ind];
				y = m_particles.pos.y[particle_ind] - m_particles.pos.y[neighbor_ind];
				z = m_particles.pos.z[particle_ind] - m_particles.pos.z[neighbor_ind];

				scalar_value = m_scalar_values[linear_ind];

				kernel_gradient_x = x*scalar_value;
				kernel_gradient_y = y*scalar_value;
				kernel_gradient_z = z*scalar_value;

				div_j = k_j / m_particles.dens[neighbor_ind];

				div_sum = (div_i + div_j);

				pressure_acc_x += m_mass * div_sum * kernel_gradient_x;
				pressure_acc_y += m_mass * div_sum * kernel_gradient_y;
				pressure_acc_z += m_mass * div_sum * kernel_gradient_z;
			}
			//pressure_force_z is not in report but it is a force and it is = F/m *delta_t
			m_particles.pred_vel.x[particle_ind] = m_particles.pred_vel.x[particle_ind] - m_delta_t * pressure_acc_x;
			m_particles.pred_vel.y[particle_ind] = m_particles.pred_vel.y[particle_ind] - m_delta_t * pressure_acc_y;
			m_particles.pred_vel.z[particle_ind] = m_particles.pred_vel.z[particle_ind] - m_delta_t * pressure_acc_z;
		}
		calculate_derived_density_pred_dens(m_neighbor_data);
		++iter;
		eta = 0.01f * density_error * C_REST_DENS;
	} while ((m_pred_dens_avg - C_REST_DENS > eta || iter < 2) && iter < Viter_max);
}

void SPH::update_positions() const
{
#pragma omp parallel
#pragma omp for
	for (int i = 0; i < current_n_particles; ++i)
	{
		m_particles.pos.x[i] += m_particles.pred_vel.x[i] * m_delta_t;
		m_particles.pos.y[i] += m_particles.pred_vel.y[i] * m_delta_t;
		m_particles.pos.z[i] += m_particles.pred_vel.z[i] * m_delta_t;
	}
}

/*
 * ViscousDFSPH, Algorithm 2
 */
void SPH::correct_divergence_error()
{
	int neighbor_ind;
	int iter{ 0 };

	float k_v_i, k_v_j, div_i, div_j, div_sum;
	float pressure_acc_x, pressure_acc_y, pressure_acc_z;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float scalar_value;
	float eta;
	float inv_delta_t = 1.f / m_delta_t;

	calculate_derived_density_pred_dens(m_neighbor_data);

	do
	{
		for (auto particle_ind = 0; particle_ind < current_n_particles; ++particle_ind)
		{
			k_v_i = 0.5f * fmax(inv_delta_t * m_dens_derive[particle_ind] * m_alpha[particle_ind], 0.5f);
			div_i = k_v_i / m_particles.dens[particle_ind];

			pressure_acc_x = pressure_acc_y = pressure_acc_z = 0.0f;

			for (auto j = 0; j < m_neighbor_data[particle_ind].n; ++j)
			{
				int linear_ind = j + D_MAX_NR_OF_NEIGHBORS*particle_ind;
				neighbor_ind = m_neighbor_data[particle_ind].neighbor[j];

				//assert(m_particles.dens[neighbor_ind] != 0.0f, "n dens");

				k_v_j = 0.5f * fmax(inv_delta_t * m_dens_derive[neighbor_ind] * m_alpha[neighbor_ind], 0.5f);

				x = m_particles.pos.x[particle_ind] - m_particles.pos.x[neighbor_ind];
				y = m_particles.pos.y[particle_ind] - m_particles.pos.y[neighbor_ind];
				z = m_particles.pos.z[particle_ind] - m_particles.pos.z[neighbor_ind];

				scalar_value = m_scalar_values[linear_ind];

				kernel_gradient_x = x*scalar_value;
				kernel_gradient_y = y*scalar_value;
				kernel_gradient_z = z*scalar_value;

				div_j = k_v_j / m_particles.dens[neighbor_ind];

				div_sum = div_i + div_j;

				pressure_acc_x += m_mass * div_sum * kernel_gradient_x;
				pressure_acc_y += m_mass * div_sum * kernel_gradient_y;
				pressure_acc_z += m_mass * div_sum * kernel_gradient_z;
			}
			//pressure_force_z is not in report but it is a force and it is = F/m *delta_t
			m_particles.pred_vel.x[particle_ind] = m_particles.pred_vel.x[particle_ind] - m_delta_t * pressure_acc_x;
			m_particles.pred_vel.y[particle_ind] = m_particles.pred_vel.y[particle_ind] - m_delta_t * pressure_acc_y;
			m_particles.pred_vel.z[particle_ind] = m_particles.pred_vel.z[particle_ind] - m_delta_t * pressure_acc_z;
		}

		calculate_derived_density_pred_dens(m_neighbor_data);

		eta = 0.01f * divergence_error * C_REST_DENS * 1.f / m_delta_t;
		++iter;
	} while (m_dens_derive_avg > eta && iter < iter_max); // implicit condition: iter < 1 
}

void SPH::update_velocities()
{
#pragma omp parallel
#pragma omp for
	for (auto i = 0; i < current_n_particles; ++i)
	{
		m_particles.vel.x[i] = m_particles.pred_vel.x[i];
		m_particles.vel.y[i] = m_particles.pred_vel.y[i];
		m_particles.vel.z[i] = m_particles.pred_vel.z[i];
	}
}

void SPH::update_density_and_factors(Neighbor_Data* neighbor_data)
{
	int nr_neighbors, neighbor_index, linear_ind;

	float abs_sum_denom, sum_abs_denom = 0;
	float denom;
	float dx, dy, dz;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;
	float scalar_value_mul_mass;
	float x = 0.f, y = 0.f, z = 0.f;
	float temporary_sum_abs;
	const float min_denom{ 10e-6f };


	//#pragma omp parallel default(shared)
	{
		
		//#pragma omp for schedule(static)  
		for (auto particle = 0; particle < current_n_particles; ++particle)
		{
			nr_neighbors = neighbor_data[particle].n;
			//added 1 * mass as particles own density as kernel is 1 at dist == 0
			//this should atleast be the case, but needs to be checked
			// => Does not seem to cause a problem when it is 0. So i followed the report

			m_particles.dens[particle] = 12.000656593f;
			for (auto neighbor = 0; neighbor < nr_neighbors; ++neighbor)
			{
				neighbor_index = neighbor_data[particle].neighbor[neighbor];

				linear_ind = neighbor + D_MAX_NR_OF_NEIGHBORS*particle;

				//Update density
				m_particles.dens[particle] += m_mass*m_kernel_values[linear_ind];

				dx = m_particles.pos.x[particle] - m_particles.pos.x[neighbor_index];
				dy = m_particles.pos.y[particle] - m_particles.pos.y[neighbor_index];
				dz = m_particles.pos.z[particle] - m_particles.pos.z[neighbor_index];

				scalar_value_mul_mass = m_mass*m_scalar_values[linear_ind];

				kernel_gradient_x = dx * scalar_value_mul_mass;
				kernel_gradient_y = dy * scalar_value_mul_mass;
				kernel_gradient_z = dz * scalar_value_mul_mass;

				x += kernel_gradient_x;
				y += kernel_gradient_y;
				z += kernel_gradient_z;

				temporary_sum_abs = sqrt(kernel_gradient_x*kernel_gradient_x + kernel_gradient_y*kernel_gradient_y + kernel_gradient_z*kernel_gradient_z);

				sum_abs_denom += temporary_sum_abs*temporary_sum_abs;
			}

			abs_sum_denom = sqrt(x*x + y*y + z*z);
			denom = abs_sum_denom*abs_sum_denom + sum_abs_denom;

			// set alpha to max(denom,min_denom)
			//denom = denom < min_denom ? min_denom : denom;
			denom = fmax(denom, 1.0e-6);
			m_alpha[particle] = m_particles.dens[particle] / denom;

			sum_abs_denom = 0.f;
			x = y = z = 0.f;
		}
	}
}

/*
 * Using a Cubic spline kernel
 * Divergence-Free SPH for Incompressible and Viscous Fluids, ref 6
*/
void update_kernel_values(float* kernel_values, Float3* pos, Neighbor_Data* neighbor_data, const int N_PARTICLES)
{
	int ind;

	float x, y, z;
	float q, q_2;
	float length;
	float search_area = D_SEARCH_RANGE;
	float pi = D_PI;
	float particle_pos_x, particle_pos_y, particle_pos_z;

	float div = 1.0f / (search_area*search_area*search_area*pi);

	//#pragma omp parallel default(shared)
	//{

	//#pragma omp for schedule(static)  
	for (auto particle = 0; particle < N_PARTICLES; ++particle)
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
			kernel_values[particle*D_MAX_NR_OF_NEIGHBORS + neighbor] = div*(1.0f - 1.5f*q_2 + 0.75f*q_2*q);;
		}
	//}
	}
}


/*
 * ViscousDFSPH, eq 9
 */
void SPH::calculate_derived_density_pred_dens(Neighbor_Data* neighbor_data)
{
	int neighbor_index, linear_ind;
	int neighbor_length;

	float pressure_derived;
	float dens_derive_sum = 0.f, pred_dens_sum = 0.f;
	float pressure_derived_x = 0.f, pressure_derived_y = 0.f, pressure_derived_z = 0.f;
	float x, y, z;
	float kernel_gradient_x, kernel_gradient_y, kernel_gradient_z;

	//#pragma omp parallel default(shared)
	//{
	//	#pragma omp for schedule(static)  
		for (int i = 0; i < current_n_particles; ++i)
		{
			neighbor_length = neighbor_data[i].n;
			for (int j = 0; j < neighbor_length; ++j)
			{
				neighbor_index = neighbor_data[i].neighbor[j];
				linear_ind = neighbor_index + D_MAX_NR_OF_NEIGHBORS*i;

				x = m_particles.pos.x[i] - m_particles.pos.x[neighbor_index];
				y = m_particles.pos.y[i] - m_particles.pos.y[neighbor_index];
				z = m_particles.pos.z[i] - m_particles.pos.z[neighbor_index];

				kernel_gradient_x = x * m_scalar_values[linear_ind];
				kernel_gradient_y = y * m_scalar_values[linear_ind];
				kernel_gradient_z = z * m_scalar_values[linear_ind];

				//equation 9 in DFSPH, changed 16-11-18
				pressure_derived_x += m_mass * kernel_gradient_x * (m_particles.pred_vel.x[i] - m_particles.pred_vel.x[neighbor_index]);
				pressure_derived_y += m_mass*kernel_gradient_y*(m_particles.pred_vel.y[i] - m_particles.pred_vel.y[neighbor_index]);
				pressure_derived_z += m_mass*kernel_gradient_z*(m_particles.pred_vel.z[i] - m_particles.pred_vel.z[neighbor_index]);
			}

			pressure_derived = pressure_derived_x + pressure_derived_y + pressure_derived_z;

			m_dens_derive[i] = pressure_derived;
			m_pred_dens[i] = m_particles.dens[i] + m_delta_t * m_dens_derive[i];

			dens_derive_sum += m_dens_derive[i];
			pred_dens_sum += m_pred_dens[i];
			pressure_derived_x = pressure_derived_y = pressure_derived_z = 0.f;
		}
	//}

	m_dens_derive_avg = dens_derive_sum / current_n_particles;
	m_pred_dens_avg = pred_dens_sum / current_n_particles;
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
	float q, dist;
	float dx, dy, dz;

	float inv_range = 1.0f / D_SEARCH_RANGE;
	float div = 1.0f / (search_area*search_area*search_area*pi);

	//#pragma omp parallel default(shared)
	//{
	
	//Loop through all particles
	//#pragma omp for schedule(static)  
	for (auto particle = 0; particle < N_PARTICLES; ++particle)
	{
		//Loop through particle neibors
		for (auto neighbor = 0; neighbor < neighbor_data[particle].n; ++neighbor)
		{
			neighbor_ind = neighbor_data[particle].neighbor[neighbor];

			dx = pos->x[particle] - pos->x[neighbor_ind];
			dy = pos->y[particle] - pos->y[neighbor_ind];
			dz = pos->z[particle] - pos->z[neighbor_ind];

			dist = sqrt(dx*dx + dy*dy + dz*dz);

			q = dist * inv_range;

			// length is always equal or smaller to D_SEARCH_RANGE => implicit intervall between [0, 1]
			kernel_derive = div*(-3.0f*q + 2.25f*q*q);
			
			scalar_values[particle*D_MAX_NR_OF_NEIGHBORS + neighbor] = kernel_derive / (search_area * dist);

		}
	}
	//}
}