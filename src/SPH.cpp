#include "SPH.h"

#define gravity 9.82

SPH::SPH()
{
}

void SPH::non_pressure_forces(){}
void SPH::calculate_time_step() {}
void SPH::predict_velocities() {}
void SPH::correct_density_error() {}
void SPH::correct_strain_rate_error() {}
void SPH::update_positions() {}
void SPH::update_neighbors() {}
//
//void find_neighborhoods();
//
void SPH::correct_divergence_error() {}
void SPH::update_velocities() {}

SPH::~SPH()
{
}
