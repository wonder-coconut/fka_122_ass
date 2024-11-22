/**
 * Calculate the acceleration.
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  masses - Vector of masses
 *  kappa - Spring constant
 *
*/
void calculate_acceleration(double *accelerations, double *positions, double *masses, double kappa)
{
    //const unsigned int N = 3;  // There are three particles

    accelerations[0] = kappa * (positions[1] - positions[0]) / masses[0];
    accelerations[1] = kappa * (positions[2] - 2 * positions[1] + positions[0]) / masses[1];
    accelerations[2] = kappa * (-positions[2] + positions[1]) / masses[2];
    
}

/**
 * Calculate the potential energy
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  positions - Vector of positions
 *  kappa - Spring constant
 * Returns
 * -------
 *  Potential energy
 *
*/
double calculate_potential_energy(double *positions, double kappa)
{
    const unsigned int N = 3;  // There are three particles

    double p_e = 0;
    for(int i = 0; i < N - 1; i++)
        p_e += (positions[i+1] - positions[i])*(positions[i+1] - positions[i]) * kappa * 0.5;

    return p_e;
}

/**
 * Calculate the kinetic energy
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  velocities - Vector of velocities
 *  masses - Vector of masses
 *
 * Returns
 * -------
 *  Kinetic energy
 *
*/
double calculate_kinetic_energy(double *velocities, double *masses)
{
    const unsigned int N = 3;  // There are three particles

    double k_e = 0;
    for(int i = 0; i < N; i++)
        k_e += masses[i] * velocities[i] * velocities[i] * 0.5;

    return k_e;
}

/**
 * Perform one velocity Verlet step
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  masses - Vector of masses
 *  kappa - Spring constant
 *  timestep - Time step
 *
*/
void velocity_verlet_one_step(double *accelerations, double *positions, double *velocities,
                              double *masses, double kappa, double timestep)
{
    const unsigned int N = 3;  // There are three particles

    // Write to accelerations, positions and velocities vectors
    int i;
    for(i = 0; i < N; i++)
        velocities[i] = velocities[i] + 0.5 * accelerations[i] * timestep;

    for(i = 0; i < N; i++)
        positions[i] = positions[i] + velocities[i] * timestep;
    
    calculate_acceleration(accelerations,positions,masses,kappa);

    for(i = 0; i < N; i++)
        velocities[i] = velocities[i] + 0.5 * accelerations[i] * timestep;
}
