#include <math.h>

double get_kinetic_energy_al(double **velocities, double al_asu, int n_atoms)
{
    double e_kin = 0;
    for(int i = 0; i < n_atoms; i++)
    {
        double v_norm = sqrt(velocities[i][0]*velocities[i][0] + velocities[i][1]*velocities[i][1] + velocities[i][2]*velocities[i][2]);
        e_kin += 0.5 * al_asu * v_norm * v_norm;
    }
    return e_kin;
}