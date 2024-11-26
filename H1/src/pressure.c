#include "tools.h"

double get_pressure_al(double volume, double kin_e, double **positions, double virial, double n_atoms)
{
    double pressure = 0;
    double p1 = 2*kin_e/(3*volume);
    double p2 = virial/volume;
    pressure = p1+p2;
    return pressure;
}