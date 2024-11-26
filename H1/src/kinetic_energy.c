#include <math.h>

double get_temperature(double kin_e)
{
    double kb_in_ev = 8.618011058E-05;
    double T = kin_e * 2 / (3 * 256 * kb_in_ev);
    return T;
}

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