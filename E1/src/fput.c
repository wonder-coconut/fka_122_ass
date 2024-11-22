/**
 * Calculate the acceleration.
 *
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  alpha - Anharmonicity constant
 *  N - Number of atoms
 *
*/
void calculate_acceleration_fput(double *accelerations, double *positions, double alpha, const unsigned int N)
{
    int i;
    for(i = 0; i < N; i++)
    {
        if(i == 0)
            accelerations[i] = (positions[i + 1] - 2 * positions[i]) * (1 + alpha * positions[i + 1]);
        else if(i == N - 1)
            accelerations[i] = (positions[i - 1] - 2 * positions[i]) * (1 - alpha * positions[i - 1]);
        else
            accelerations[i] = (positions[i + 1] - 2 * positions[i] + positions[i - 1]) * (1 + alpha * (positions[i+1] - positions[i - 1]));
    }
}

/**
 * Perform one velocity Verlet step
 *
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  alpha - Anharmonicity constant
 *  timestep - Time step
 *  N - Number of atoms
 *
*/
void velocity_verlet_one_step_fput(double *accelerations, double *positions, double *velocities,
                              double alpha, double timestep, const unsigned int N)
{
    // Write to accelerations, positions and velocities vectors
    int i;
    for(i = 0; i < N; i++)
        velocities[i] = velocities[i] + 0.5 * accelerations[i] * timestep;

    for(i = 0; i < N; i++)
        positions[i] = positions[i] + velocities[i] * timestep;
    
    calculate_acceleration_fput(accelerations,positions,alpha,N);

    for(i = 0; i < N; i++)
        velocities[i] = velocities[i] + 0.5 * accelerations[i] * timestep;
}
