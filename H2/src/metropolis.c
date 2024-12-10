#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "tools.h"
#include "mcmc_step.h"
#include "trial_wave.h"

double q_calc(double *r1_new, double *r2_new, double *r1, double *r2, double alpha)
{

    //new parameters
    double m_r12_new = distance_between_vectors(r2_new,r1_new,3);

    double m_r1_new = vector_norm(r1_new,3);
    double m_r2_new = vector_norm(r2_new,3);

    //original parameters
    double m_r12 = distance_between_vectors(r2,r1,3);

    double m_r1 = vector_norm(r1,3);
    double m_r2 = vector_norm(r2,3);

    double exp_param_new = (-2*m_r1_new) + (-2*m_r2_new) + (m_r12_new/(2 * (1 + alpha * m_r12_new)));
    double exp_param = (-2*m_r1) + (-2*m_r2) + (m_r12/(2 * (1 + alpha * m_r12)));

    double q_exp_param = 2*(exp_param_new - exp_param);
    return exp(q_exp_param);
}

double pdf(double *r1, double *r2, double alpha)
{
    double psi = trial_wave_function(r1, r2, alpha);
    return psi*psi;
}

double energy(double *r1, double *r2, double alpha)
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
    double unit_r2_dot_r21 = dot_product(unit_r21,r21,3);

    double e_t1 = unit_r2_dot_r21/(m_r12 * (1+alpha*m_r12)*(1+alpha*m_r12));
    double e_t2 = 1/(m_r12 * (1+alpha*m_r12)*(1+alpha*m_r12)*(1+alpha*m_r12));
    double e_t3 = 0.25/((1+alpha*m_r12)*(1+alpha*m_r12)*(1+alpha*m_r12)*(1+alpha*m_r12));
    double e_t4 = 1/m_r12;

    return -4 + e_t1 - e_t2 - e_t3 +e_t4;
}

mcmc_step mcmc_displace_all(double *r1, double *r2, double d, double alpha, gsl_rng *r)
{   //struct def
    mcmc_step iteration;

    //new walker step init
    double *r1_new = (double *)malloc(3 * sizeof(double));
    double *r2_new = (double *)malloc(3 * sizeof(double));

    //mcmc walk
    r1_new[0] = r1[0] + d * (gsl_rng_uniform(r) - 0.5);
    r1_new[1] = r1[1] + d * (gsl_rng_uniform(r) - 0.5);
    r1_new[2] = r1[2] + d * (gsl_rng_uniform(r) - 0.5);

    r2_new[0] = r2[0] + d * (gsl_rng_uniform(r) - 0.5);
    r2_new[1] = r2[1] + d * (gsl_rng_uniform(r) - 0.5);
    r2_new[2] = r2[2] + d * (gsl_rng_uniform(r) - 0.5);

    //acceptance rubric
    double q = q_calc(r1_new,r2_new,r1,r2,alpha);
    //double q = pdf(r1_new,r2_new,alpha)/pdf(r1,r2,alpha);
    //printf("q:\t%f\n",q);
    double epsilon = gsl_rng_uniform(r);

    if(q >= epsilon)
    {//confirm walk
        r1[0] = r1_new[0];
        r1[1] = r1_new[1];
        r1[2] = r1_new[2];

        r2[0] = r2_new[0];
        r2[1] = r2_new[1];
        r2[2] = r2_new[2];

        iteration.accepted = 1;
    }
    else
        iteration.accepted = 0;

    iteration.probability = pdf(r1,r2,alpha);
    iteration.energy = energy(r1,r2,alpha);

    free(r1_new);
    free(r2_new);

    return iteration;
}