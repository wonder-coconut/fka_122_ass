#include <math.h>

#include "tools.h"

double ln_psi_gradient(double *r1, double *r2, double alpha)
{
    double m_r12 = distance_between_vectors(r1,r2,3);

    double res = -0.5 * (m_r12*m_r12)/((1 + alpha * m_r12) * (1 + alpha * m_r12));

    return res;
}

double alpha_descent_one_step(double alpha, double E_gradient, int p, double A, double beta)
{
    //learning rate
    double gamma = A * pow(p, -beta);
    //double gamma = 0.1;
    //alpha descent
    double alpha_res = alpha - gamma * E_gradient;
    return alpha_res;
}