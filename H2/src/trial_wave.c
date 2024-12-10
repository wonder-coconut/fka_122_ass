#include <math.h>

#include "tools.h"

double trial_wave_function(double *r1, double *r2, double alpha)
{
    //interelectron distance
    double m_r12 = distance_between_vectors(r2,r1,3);
    
    double m_r1 = vector_norm(r1,3);
    double m_r2 = vector_norm(r2,3);

    double exp_param = (-2*m_r1) + (-2*m_r2) + (m_r12/(2 * (1 + alpha * m_r12)));
    //printf("exp param:\t%f\n",exp_param);
    return exp(exp_param);
}