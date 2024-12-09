#include <stdio.h>
#include <stdlib.h>

#include "potential.h"

//kinetic energy calculation using velocities
double calculate_kinetic_energy(double **velocities, double mass, int n_atoms)
{
    double k_e = 0;

    for(int i = 0; i < n_atoms; i++)
    {
        k_e +=  velocities[i][0] * velocities[i][0];
        k_e +=  velocities[i][1] * velocities[i][1];
        k_e +=  velocities[i][2] * velocities[i][2];
    }

    return k_e * 0.5 * mass;
}

//single step of the velocity verlet iteration
void vv_one_step(double **forces, double **velocities, double **positions, double *pot_e, double *virial, double timestep, double cell_length, int n_atoms, double mass)
{
    int i;
    for(i = 0; i < n_atoms; i++)
    {
        velocities[i][0] = velocities[i][0] + 0.5 * (forces[i][0] / mass) * timestep;
        velocities[i][1] = velocities[i][1] + 0.5 * (forces[i][1] / mass) * timestep;
        velocities[i][2] = velocities[i][2] + 0.5 * (forces[i][2] / mass) * timestep;
    }
    
    for(i = 0; i < n_atoms; i++)
    {
        positions[i][0] += velocities[i][0] * timestep;
        positions[i][1] += velocities[i][1] * timestep;
        positions[i][2] += velocities[i][2] * timestep;
    }
    
    //update forces and potential e
    calculate(pot_e, virial, forces, positions, cell_length, n_atoms);

    for(i = 0; i < n_atoms; i++)
    {
        velocities[i][0] += 0.5 * (forces[i][0] / mass) * timestep;
        velocities[i][1] += 0.5 * (forces[i][1] / mass) * timestep;
        velocities[i][2] += 0.5 * (forces[i][2] / mass) * timestep;
    }

}
