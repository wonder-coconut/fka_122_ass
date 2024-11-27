#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "tools.h"
#include "result_t_temp.h"

#define PI 3.141592653589

double f_x(double x)
{
    return x - x*x;
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
result_t_temp MC_without_importance_sampling(int N, gsl_rng *k)
{
    //func = x * (1 - x)
    result_t_temp result;

    double *x = (double *)malloc(N * sizeof(double));

    int i;
    for(i = 0 ; i < N; i++)
        x[i] = gsl_rng_uniform(k);
    
    double temp,integral_temp,integral_error,integral_temp_squared, variance;
    temp = integral_temp = integral_error = integral_temp_squared = variance = 0;

    for(i = 0 ; i < N ; i++)
    {
        temp = f_x(x[i]);
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
