#include <math.h>
#include "kinetic_energy.h"

double calc_alpha_T(double timestep, double time_constant_T, double T_inst, double T_eq)
{
    return  1 + 2 * timestep * (T_eq - T_inst) / (time_constant_T * T_inst);
}

double calc_alpha_P(double timestep, double k_T, double time_constant_P, double P_inst, double P_eq)
{
    return 1 - k_T * timestep * (P_eq - P_inst) / time_constant_P;
}

void velocity_update(double **velocities, double timestep, double time_constant_T, double T_inst, double T_eq, double n_atoms)
{
    double alpha = calc_alpha_T(timestep, time_constant_T, T_inst, T_eq);
    double alpha_scale = pow(alpha,0.5);
        
    for(int i = 0; i < n_atoms; i++)
    {
        velocities[i][0] = alpha_scale * velocities[i][0];
        velocities[i][1] = alpha_scale * velocities[i][1];
        velocities[i][2] = alpha_scale * velocities[i][2];
    }
}

void pressure_update(double **positions, double k_T, double timestep, double time_constant_P, double P_inst, double P_eq, double n_atoms)
{
    double alpha = calc_alpha_P(timestep, k_T, time_constant_P, P_inst, P_eq);
    double alpha_scale = pow(alpha, 1.0/3);
    for(int i = 0; i < n_atoms; i++)
    {
        positions[i][0] = alpha_scale * positions[i][0];
        positions[i][1] = alpha_scale * positions[i][1];
        positions[i][2] = alpha_scale * positions[i][2];
    }
}