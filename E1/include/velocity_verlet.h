void calculate_acceleration(double *accelerations, double *positions, double *masses, double kappa);

double calculate_potential_energy(double *positions, double kappa);

double calculate_kinetic_energy(double *velocities, double *masses);

void velocity_verlet_one_step(double *accelerations, double *positions, double *velocities,
                              double *masses, double kappa, double timestep);