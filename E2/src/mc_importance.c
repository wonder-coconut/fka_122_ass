#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "tools.h"
#include "result_t_temp.h"

#define PI 3.141592653589

double f(double x)
{
    return x - x*x;
}

double p(double x)
{
    return PI * 0.5 *sin(PI * x);
}

double g(double x)
{
    return f(x)/p(x);
}

double F_1(double x)
{
    return acos(1 - 2 * x)/PI;
}

/*
 * The following struct is defined in the header file
 * Use it to store the result

typedef struct {
    double integral;
    double error;
} result_t;
 */


/* **********************************************
 *
 * Perform Monte Carlo integration of the given
 * integral without using importance sampling.
 *
 * Parameters
 * ----------
 *  N - Number of points to sample
 *  k - GSL random number generator object
 *
 * Returns
 * -------
 *  Struct with the result
 *
 * **********************************************/
result_t_temp MC_with_importance_sampling(int N, gsl_rng *k)
{
    //func = x * (1 - x)
    result_t_temp result;

    double *x = (double *)malloc(N * sizeof(double));

    int i;
    double rtemp;
    for(i = 0 ; i < N; i++)
    {
        rtemp = gsl_rng_uniform(k);
        x[i] = F_1(rtemp);
    }
    double temp,integral_temp,integral_error,integral_temp_squared, variance;
    temp = integral_temp = integral_error = integral_temp_squared = variance = 0;

    for(i = 0 ; i < N ; i++)
    {
        temp = g(x[i]);
        integral_temp += temp;
        integral_temp_squared += temp*temp; 
    }
    
    integral_temp_squared /= N;
    integral_temp /= N;

    variance = integral_temp_squared - (integral_temp * integral_temp);

    result.integral = integral_temp; // Write the integral here
    result.error = sqrt(variance/N);// Write the error here

    free(x);

    return result;
}
