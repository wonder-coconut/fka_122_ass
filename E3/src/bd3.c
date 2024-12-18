#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

#include "result_t.h"

#define PI 3.141592653589


double acceleration(double w0, double x)
{
    return -w0*w0*x;
}

/*
 * The following struct is defined in the header file
 * Use it to store the result

typedef struct {
    double position;
    double velocity; 
} result_t;
 */


/* **********************************************
 *
 * Perform one step of the BD3 algorithm and
 * store the positions and velocities in result
 *
 * Parameters
 * ----------
 *  initial_position  - the initial position of the 
 *                      particle
 *  initial_veclotity - the initial velocity of the 
 *                      particle
 *  w0                - the angular frequency of the
 *                      optical trap
 *  dt                - the timestep
 *  eta               - 1 / relaxation time
 *  kB                - Boltzmann constant
 *  mass              - mass of particle
 *  T                 - The temperature
 *  k                 - GSL random number generator object
 *
 * Returns
 * -------
 *  Struct with the result
 *
 * **********************************************/
result_t BD3(double initial_position, double initial_velocity, double w0, double dt, double eta, double kB, double mass, double T, gsl_rng *k)
{
    result_t result;

    double c_0 = exp(-eta*(dt));
    double vth = sqrt(kB*T/mass);

    //initial calc
    double a = acceleration(w0, initial_position);

    double v_temp = 0.5*a*dt + sqrt(c_0)*initial_velocity + vth*sqrt(1-c_0)*gsl_ran_gaussian(k,1);
    double x_t = initial_position + v_temp*dt;

    //accl calc
    double a_t = acceleration(w0, x_t);

    double v_t = 0.5*sqrt(c_0)*a_t*dt + sqrt(c_0)*v_temp + vth*sqrt(1-c_0)*gsl_ran_gaussian(k,1);

    result.position = x_t;
    result.velocity = v_t;
    return result;
}


