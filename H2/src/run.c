#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "tools.h"
#include "trial_wave.h"
#include "metropolis.h"
#include "mcmc_step.h"
#include "inefficiency.h"


void task2_block_avg_inefficiency()
{
    //energy samples input
    FILE *ip_energy = fopen("op_text/mcmc_energy_eqb_t2.txt","r");
    //block inefficiency output
    FILE *op_blocks = fopen("op_text/block_avg_ineff_t2.txt","w");

    //n samples
    int n;
    fscanf(ip_energy,"%d\n",&n);

    double *data = (double *)malloc(n * sizeof(double));
    //read data
    int i;
    for(i = 0; i < n ; i++)
        fscanf(ip_energy,"%lf\n",&data[i]);
    
    double ineff = 0;
    
    //block averaging method of inefficiency calculation
    for(int b = 2; b < 1000; b++)
    {
        ineff = block_average(data, n, b);
        fprintf(op_blocks,"%f\n",ineff);
    }

    free(data);
}

void task2_autocorrelation_inefficiency()
{
    //energy samples input
    FILE *ip_energy = fopen("op_text/mcmc_energy_eqb_t2.txt","r");
    //block inefficiency output
    FILE *op_autocorr = fopen("op_text/autocorrelation_ineff_t2.txt","w");

    //n samples
    int n;
    fscanf(ip_energy,"%d\n",&n);

    double *data = (double *)malloc(n * sizeof(double));
    //read data
    int i;
    for(i = 0; i < n ; i++)
        fscanf(ip_energy,"%lf\n",&data[i]);
    
    //mean bias
    double f_mean = average(data,n);
    for(i = 0; i < n ; i++)
        data[i] -= f_mean;

    //autocorrelation
    for(i = 1; i < 200; i++)
        fprintf(op_autocorr,"%f\n",autocorrelation(data,n,i));

    free(data);

}

void task2_sampling(int argc, char *argv[], gsl_rng *r)
{
    //file output
    FILE *op_energy = fopen("op_text/mcmc_energy_t2.txt","w");
    FILE *op_samples = fopen("op_text/samples.txt","w");

    int N = atoi(argv[2]);
    double d = atof(argv[3]);

    //initial states
    double r1[] = {-10000,10000,-10000};
    double r2[] = {10000,-10000,10000};
    double alpha = 0.1;

    mcmc_step iteration;

    int i;
    int accepted_iterations = 0;
    double energy = 0;

    //mcmc sampling
    for(i = 0; i < N; i++)
    {
        iteration = mcmc_displace_all(r1,r2,d,alpha,r);
        fprintf(op_energy,"%f\n",iteration.energy);
        //if(i%100 == 0)
            fprintf(op_samples,"%f\t%f\t%f\t%f\t%f\t%f\n",r1[0],r1[1],r1[2],r2[0],r2[1],r2[2]);
        
        accepted_iterations += iteration.accepted;
        energy += iteration.energy;
    }

    printf("Acceptance ratio:\t%f\n",accepted_iterations*1.0/N);
    printf("Energy:\t\t%f\n",(energy/N));
}

void task1_correlation()
{
    //file input
    FILE *samples = fopen("op_text/samples.txt","r");
    //file output
    FILE *op_theta = fopen("op_text/theta_t1.txt","w");
    FILE *op_uniform = fopen("op_text/x_uniform_t1.txt","w");

    int n;
    fscanf(samples,"%d\n",&n);

    //initial states
    double r1[] = {0,0,0};
    double r2[] = {0,0,0};
    
    double dot_r1r2 = 0;
    double m_r1 = 0;
    double m_r2 = 0;
    double theta = 0;
    double x = 0;

    //reading the file iteratively
    int i;
    for(i = 0; i < n; i++)
    {
        fscanf(samples,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&r1[0],&r1[1],&r1[2],&r2[0],&r2[1],&r2[2]);
        dot_r1r2 = dot_product(r1,r2,3);
        m_r1 = vector_norm(r1,3);
        m_r2 = vector_norm(r2,3);
        x = dot_r1r2/(m_r1 * m_r2);
        theta = acos(x);
        fprintf(op_theta,"%f\n",theta);
        fprintf(op_uniform,"%f\n",x);
    }
}

void task1_sampling(int argc,char *argv[], gsl_rng *r)
{
    //file output
    FILE *samples = fopen("op_text/samples.txt","w");

    int N = atoi(argv[2]);
    double d = atof(argv[3]);

    //initial states
    double r1[] = {0,0,0};
    double r2[] = {0,0,0};
    double alpha = 0.1;

    mcmc_step iteration;

    int i;
    int init_limit = pow(10,5);
    int accepted_iterations = 0;
    double energy = 0;

    //init run
    for(i = 0; i < init_limit; i++)
        iteration = mcmc_displace_all(r1,r2,d,alpha,r);
    
    //print length
    fprintf(samples,"%d\n",N);

    //sample run
    for(i = 0; i < N; i++)
    {
        iteration = mcmc_displace_all(r1,r2,d,alpha,r);

        fprintf(samples,"%f\t%f\t%f\t",r1[0],r1[1],r1[2]);
        fprintf(samples,"%f\t%f\t%f\n",r2[0],r2[1],r2[2]);

        energy += iteration.energy;
        accepted_iterations += iteration.accepted;
    }
    printf("Acceptance ratio:\t%f\n",accepted_iterations*1.0/N);
    printf("Average energy:\t\t%f\n",energy/N);
}

int
run(int argc, char *argv[])
{
    unsigned long epochs = (unsigned long) time (NULL); 

    //gsl setup
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_rng_set(r, epochs);

    //task split
    int choice = atoi(argv[1]);
    if(choice == 1)
        task1_sampling(argc,argv, r);
    else if(choice == 2)
        task1_correlation();
    else if(choice == 3)
        task2_sampling(argc, argv, r);
    else if(choice == 4)
        task2_autocorrelation_inefficiency();
    else if(choice == 5)
        task2_block_avg_inefficiency();

    gsl_rng_free(r);

    return 0;
}
