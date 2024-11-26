#pragma once
/*
Gives the instantaneous temperature of the lattice
Input:
kin_e : kinetic energy of lattice (double)   
*/

double get_temperature(double);

/*
Gives the kinetic energy of the lattice
Input:
velocities :    N * 3 double matrix
al_asu :        Mass of Al in ASU (double)
n_atoms :       number of atoms in lattice (int)
*/
double get_kinetic_energy_al(double**, double , int);