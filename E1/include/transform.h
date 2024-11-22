#pragma once

void gen_transform_mat(double **transform_mat, const unsigned int N);

void transform_to_normal_modes(double *positions, double **transform_mat, double *Q, const unsigned int N);

void calculate_normal_mode_energies(double *energies, double *positions, double **transform_mat, double *velocities, const unsigned int N);
