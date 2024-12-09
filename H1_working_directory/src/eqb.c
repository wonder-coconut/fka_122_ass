#include <math.h>

#include "potential.h"
#include "vv.h"
#include "tools.h"

#define KB 8.617341273302701E-05
#define KT 2.102594007874016

//calculate instantaneous temperature as a function of kinetic energy
double inst_temp_calc(double kin_e, double n_atoms)
{
    return 2/(3 * n_atoms * KB) * kin_e;
}

//calculate instantaneous pressure as a function of the virial, volume, and kinetic energy
double inst_press_calc(double kin_e, double virial, double volume)
{
    return (2*kin_e/3 + virial)/volume;
}

//alpha factor calculation for T equilibration
double alpha_T_calc(double T_eq , double T_inst, double T_const, double timestep)
{
    return 1 + 2 * timestep * (T_eq - T_inst) / (T_inst * T_const);
}

//alpha factor calculation for P equilibration
double alpha_P_calc(double P_eq, double P_inst, double P_const, double timestep)
{
    return 1 - KT * timestep * (P_eq - P_inst) / P_const;
}

//temperature equilibration 1 step
void eqb_temp_1_step(double **velocities, double T_eq, double T_const, double timestep, double kin_e, double mass, int n_atoms)
{
    double T_inst = inst_temp_calc(kin_e, n_atoms);//instantaneous temp

    double alpha_T = alpha_T_calc(T_eq, T_inst, T_const, timestep);
    double alpha_T_scaled = pow(alpha_T,0.5);//alpha ^ 1/2

    for(int i = 0; i < n_atoms; i++)
    {
        //velocity scaling
        velocities[i][0] *= alpha_T_scaled;
        velocities[i][1] *= alpha_T_scaled;
        velocities[i][2] *= alpha_T_scaled;
    }
}

//pressure equilibration 1 step
void eqb_press_1_step(double **positions, double virial, double P_eq, double P_const, double timestep, double kin_e, double mass, int n_atoms, double *lp)
{
    double L = 4 * *lp; //box size
    double volume = L * L * L;
    double P_inst = inst_press_calc(kin_e, virial, volume); //instantaneous pressure

    double alpha_P = alpha_P_calc(P_eq, P_inst, P_const, timestep);
    double alpha_P_scaled = pow(alpha_P, 1/3.0);//alpha ^ 1/3

    *lp = *lp * alpha_P_scaled; //lattice parameter scaling
    for(int i = 0; i < n_atoms; i++)
    {
        //position scaling
        positions[i][0] *= alpha_P_scaled;
        positions[i][1] *= alpha_P_scaled;
        positions[i][2] *= alpha_P_scaled;
    }
}