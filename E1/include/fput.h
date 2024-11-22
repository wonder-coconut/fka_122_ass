#pragma once

void calculate_acceleration_fput(double *accelerations, double *positions, double alpha, const unsigned int N);

void velocity_verlet_one_step_fput(double *accelerations, double *positions, double *velocities,
                              double alpha, double timestep, const unsigned int N);