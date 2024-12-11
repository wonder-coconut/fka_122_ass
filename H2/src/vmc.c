#include <math.h>

#include "tools.h"

double energy_gradient(double *r1, double *r2, double alpha)
{ 
    double unit_r1[3];
    double unit_r2[3];

    double unit_r21[3];
    double r21[3];

    double m_r12 = distance_between_vectors(r1,r2,3);

    //copying ri to unit_ri
    addition_with_constant(unit_r1, r1, 0, 3);
    addition_with_constant(unit_r2, r2, 0, 3);

    //unit position vecotr
    normalize_vector(unit_r1,3);
    normalize_vector(unit_r1,3);

    //init vectors
    elementwise_subtraction(unit_r21,r1,r2,3);
    elementwise_subtraction(r21,r1,r2,3);

    //energy
    double unit_r21_dot_r21 = dot_product(unit_r21,r21,3);

    //gradient of E wrt alpha
    double e2 = -2 * m_r12 * pow(1 + alpha * m_r12,-3);
    e2 *= unit_r21_dot_r21/m_r12;
    double e3 = -3 * m_r12 * pow(1 + alpha * m_r12,-4);
    e3 *= 1/m_r12;
    double e4 = -1 * m_r12 * pow(1 + alpha * m_r12,-5);

    double e_grad = e2 - e3 - e4;
    return e_grad;
}

double alpha_descent_one_step(double alpha, double E_gradient, int p, double A, double beta)
{
    //learning rate
    double gamma = A * pow(p, -beta);
    //alpha descent
    double alpha_res = alpha - gamma * E_gradient;
    return alpha_res;
}