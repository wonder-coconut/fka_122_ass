#include <stddef.h>

#include "potential.h"

void velocity_verlet_one_step_lattice(double **forces ,double **velocities ,double **positions ,double timestep, int N, double lp ,int n_atoms, double al_mass, double *pot_e)
{
    int i;
    for(i = 0; i < n_atoms; i++)
    {
        velocities[i][0] = velocities[i][0] + 0.5 * forces[i][0] * timestep / al_mass;
        velocities[i][1] = velocities[i][1] + 0.5 * forces[i][1] * timestep / al_mass;
        velocities[i][2] = velocities[i][2] + 0.5 * forces[i][2] * timestep / al_mass;
    }

    for(i = 0; i < n_atoms; i++)
    {
        positions[i][0] = positions[i][0] + velocities[i][0] * timestep;
        positions[i][1] = positions[i][1] + velocities[i][1] * timestep;
        positions[i][2] = positions[i][2] + velocities[i][2] * timestep;
    }

    //get_forces_AL(forces, positions, N*lp, n_atoms);
    calculate(pot_e,NULL,forces,positions,N*lp,n_atoms);

    for(i = 0; i < n_atoms; i++)
    {
        velocities[i][0] = velocities[i][0] + 0.5 * forces[i][0] * timestep / al_mass;
        velocities[i][1] = velocities[i][1] + 0.5 * forces[i][1] * timestep / al_mass;
        velocities[i][2] = velocities[i][2] + 0.5 * forces[i][2] * timestep / al_mass;
    }
}

void velocity_verlet_one_step_lattice_no_pot_e(double **forces ,double **velocities ,double **positions ,double timestep, int N, double lp ,int n_atoms, double al_mass)
{
    int i;
    for(i = 0; i < n_atoms; i++)
    {
        velocities[i][0] = velocities[i][0] + 0.5 * forces[i][0] * timestep / al_mass;
        velocities[i][1] = velocities[i][1] + 0.5 * forces[i][1] * timestep / al_mass;
        velocities[i][2] = velocities[i][2] + 0.5 * forces[i][2] * timestep / al_mass;
    }

    for(i = 0; i < n_atoms; i++)
    {
        positions[i][0] = positions[i][0] + velocities[i][0] * timestep;
        positions[i][1] = positions[i][1] + velocities[i][1] * timestep;
        positions[i][2] = positions[i][2] + velocities[i][2] * timestep;
    }

    get_forces_AL(forces, positions, N*lp, n_atoms);
    //calculate(pot_e,NULL,forces,positions,N*lp,n_atoms);

    for(i = 0; i < n_atoms; i++)
    {
        velocities[i][0] = velocities[i][0] + 0.5 * forces[i][0] * timestep / al_mass;
        velocities[i][1] = velocities[i][1] + 0.5 * forces[i][1] * timestep / al_mass;
        velocities[i][2] = velocities[i][2] + 0.5 * forces[i][2] * timestep / al_mass;
    }
}