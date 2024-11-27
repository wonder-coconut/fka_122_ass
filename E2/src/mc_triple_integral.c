#include <math.h>
#include <gsl/gsl_rng.h>

#include "result_t.h"
//#include "tools.h"

#define PI 3.141592653589

/*
 * The following struct is defined in the header file
 * Use it to store the result

typedef struct{
    double probability;
    double function_value;
    int accepted;
} result_t;
*/

/* ***************************************
*
* Perfmore the calculation of the
* unormalized probability function, i.e.,
* the weight used in the MCMC routine
*
* Parameters
* ----------
*  x - current walker position, size = 3
*
* Returns
* -------
* The weight of the current position
*
* ***************************************/
double weight(double *x)
{
    return exp(-(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
}

/* ***************************************
*
* Perfmore the calculation of the
* function to be sampled
*
* Parameters
* ----------
*  x - current walker position, size = 3
*
* Returns
* -------
* The function value
*
* ***************************************/
double function(double *x)
{
    return x[0]*x[0] + x[0]*x[0]*x[1]*x[1] + x[0]*x[0]*x[1]*x[1]*x[2]*x[2];
}

/* ***************************************
*
* Perfmore the calculation of the
* function to be sampled
*
* - Use gsl_rng_uniform to displace the
*   walker.
* - The walker should be displaced
*   in the order x,y,z.
* - The draw for the acceptance condition
*   should be done after the draws
*   for the displacing of the walker
*
*
* Parameters
* ----------
*  x - current walker position, size = 3
*  delta - step size
*  k - GSL random number generator object
*
* Returns
* -------
* - The function should update the x parameter
*   to reflect if the move was accepted or
*   rejected.
* - result should contain the probability
*   of the "exiting" x parameter.
* - result should should contain the function
*   value of the "exiting" x parameter.
* - result should be 1 if the move was accepted
*   and 0 if it was rejected.
*
* ***************************************/
 result_t MCMC_step_displace_all(
                                 double *x,
                                 double delta,
                                 gsl_rng *k)
{
    result_t result;

    double *x_new = (double *)malloc(3 * sizeof(double));
    double scale = pow(PI, -1.5);
    
    x_new[0] = x[0] + delta * (gsl_rng_uniform(k) - 0.5);
    x_new[1] = x[1] + delta * (gsl_rng_uniform(k) - 0.5);
    x_new[2] = x[2] + delta * (gsl_rng_uniform(k) - 0.5);

    double q = weight(x_new) / weight(x);
    double epsilon = gsl_rng_uniform(k);

    if(q >= epsilon)
    {
        x[0] = x_new[0];
        x[1] = x_new[1];
        x[2] = x_new[2];

        result.accepted = 1;
    }
    else
        result.accepted = 0;

    result.probability = scale * weight(x);
    result.function_value = function(x);

    free(x_new);
    return result;
}