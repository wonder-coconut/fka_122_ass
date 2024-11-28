#pragma once
double calc_alpha_T(double timestep, double time_constant_T, double T_inst, double T_eq);

double calc_alpha_P(double timestep, double k_T, double time_constant_P, double P_inst, double P_eq);

void velocity_update(double **velocities, double timestep, double time_constant_T, double T_inst, double T_eq, double n_atoms);

void pressure_update(double **positions, double k_T, double timestep, double time_constant_P, double P_eq, double P_inst, double n_atoms, double *lp);
