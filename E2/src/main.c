#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>

#include "mc_uniform.h"
#include "mc_importance.h"
#include "result_t_temp.h"
#include "result_t.h"
#include "mc_triple_integral.h"

void part2()
{
    //file output
    //FILE *f_out = fopen("op_text/triple_integral_stats.txt","w");

    //rng setup
    const gsl_rng_type * T;
    gsl_rng * r;

    // Initialize RNG
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Set seed
    gsl_rng_set(r, 299792458);

    result_t res;

    double *x = (double *)malloc(3 * sizeof(double));
    
    x[0] = gsl_rng_uniform(r);
    x[1] = gsl_rng_uniform(r);
    x[2] = gsl_rng_uniform(r);

    int N = 10;
    double delta = 0.5;
    int i,j;
    double integral_sum = 0;

    for(i = 0; i < 10; i++)
    {
        printf("%d\t",i);
        N = pow(10,i);
        for(j = 0; j < N ; j++)
        {
            res = MCMC_step_displace_all(x, delta, r);
            integral_sum += res.function_value;
        }
        printf("%f\n",integral_sum/N);
    }
}

void part1()
{

    //file output
    FILE *f_out = fopen("op_text/integral_stats.txt","w");

    //rng setup
    const gsl_rng_type * T;
    gsl_rng * r;

    // Initialize RNG
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Set seed
    gsl_rng_set(r, 299792458);

    //integral calc
    result_t_temp res;
    int N,i;

    for(i = 0 ; i < 5; i++)
    {
        N = pow(10,i+1);
        res = MC_without_importance_sampling(N, r);
        fprintf(f_out,"%f\t%f\t",res.integral,res.error);
        
        res = MC_with_importance_sampling(N, r);
        fprintf(f_out,"%f\t%f\n",res.integral,res.error);
    }
    gsl_rng_free(r);
}

int main(int argc, char *argv[])
{
    int choice = atoi(argv[1]);
    if(choice == 1)
        part1();
    else if(choice == 2)
        part2();

    return 0;
}