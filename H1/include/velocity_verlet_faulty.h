#pragma once

/* function that exectutes one step of the velocity verlet algorithm
input params:
forces :        n_atoms * 3 double matrix 
velocities :    n_atoms * 3 double matrix
positions :     n_atoms * 3 double matrix
timestep :      double
N :             int
lp :            double
n_atoms :       int
al_mass :       double
pot_e :         double pointer*/

void velocity_verlet_one_step_lattice(double **,double **,double **,double,int,double,int, double, double *);