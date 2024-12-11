#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "tools.h"
#include "trial_wave.h"
#include "metropolis.h"
#include "mcmc_step.h"
#include "inefficiency.h"
#include "vmc.h"

void task5(int argc, char *argv[], gsl_rng *r)
{
    //FILE output
    //FILE *energy_alpha_t4 = fopen("op_text/energy_evol_alpha_t4.txt","w");
    //hard coded statistical inefficiency
    //double n_s = 36;
    //limits
    int N = atoi(argv[2]);
    double d = atof(argv[3]);

    //initial states
    double r1[] = {0,0,0};
    double r2[] = {0,0,0};

    mcmc_step iteration;

    //simulation parameters
    int i;
    int eqb_limit = 1E5;
    double e = 0;
    double hardcoded_alpha = 0.249966;

    //init random positions
    r2[0] = (100*d) * (gsl_rng_uniform(r) - 0.5);
    r2[1] = (100*d) * (gsl_rng_uniform(r) - 0.5);
    r2[2] = (100*d) * (gsl_rng_uniform(r) - 0.5);

    r2[0] = (100*d) * (gsl_rng_uniform(r) - 0.5);
    r2[1] = (100*d) * (gsl_rng_uniform(r) - 0.5);
    r2[2] = (100*d) * (gsl_rng_uniform(r) - 0.5);

    //eqb init run
    for(i = 0 ; i < eqb_limit; i++)
        iteration = mcmc_displace_all(r1,r2,d,hardcoded_alpha,r);
        
    //sampling run
    for(i = 0; i < N; i++)
    {
        if(i%(N/100) == 0)
            printf("%d\n",(int)(i*100.0/N));
        iteration = mcmc_displace_all(r1,r2,d,hardcoded_alpha,r);
        e += iteration.energy;
    }

    //energy calculation
    e = e/N;
    printf("%f\t%f\n",hardcoded_alpha,e);

} 


//input format : ./H2 7 N   d   i_limit alpha_init  beta
void task4(int argc, char *argv[], gsl_rng *r)
{
    //FILE output
    //FILE *energy_alpha_t4 = fopen("op_text/energy_evol_alpha_t4.txt","w");
    //hard coded statistical inefficiency
    //double n_s = 36;
    //limits
    int N = atoi(argv[2]);
    double d = atof(argv[3]);

    //energy array
    double *energy = (double *)malloc(N * sizeof(double));
    double *energy_grad = (double *)malloc(N * sizeof(double));

    //initial states
    double r1[] = {0,0,0};
    double r2[] = {0,0,0};
    double alpha;

    mcmc_step iteration;

    //simulation parameters
    int i,p;
    int eqb_limit = 1E5;
    int p_limit = atoi(argv[4]);
    double e = 0;
    double e_grad = 0;
    //double e_err = 0;
    //double sd = 0;
    double alpha_init = atof(argv[5]);
    double A = 1;
    double beta = atof(argv[6]);

    alpha = alpha_init;

    //init energy arr
    for(i = 0; i < N; i++)
    {
        energy[i] = 0;
        energy_grad[i] = 0;
    }

    //iterating through different alpha values 
    for(p = 0; p < p_limit; p++)
    {
        //printf("%d\t%d\n",p,p_limit);
        //init random positions
        r2[0] = (100*d) * (gsl_rng_uniform(r) - 0.5);
        r2[1] = (100*d) * (gsl_rng_uniform(r) - 0.5);
        r2[2] = (100*d) * (gsl_rng_uniform(r) - 0.5);

        r2[0] = (100*d) * (gsl_rng_uniform(r) - 0.5);
        r2[1] = (100*d) * (gsl_rng_uniform(r) - 0.5);
        r2[2] = (100*d) * (gsl_rng_uniform(r) - 0.5);

        //energy calc using MCMC for alpha
        //eqb init run
        for(i = 0 ; i < eqb_limit; i++)
            iteration = mcmc_displace_all(r1,r2,d,alpha,r);
        
        //sampling run
        for(i = 0; i < N; i++)
        {
            iteration = mcmc_displace_all(r1,r2,d,alpha,r);
            energy[i] = iteration.energy;
            energy_grad[i] = energy_gradient(r1,r2,alpha);
        }

        //energy calculation
        e = average(energy,N);
        printf("%f\t%f\n",alpha,e);

        //gradient descent
        e_grad = average(energy_grad,N);
        printf("gradient expectation:\t%f\n",e_grad);
        alpha = alpha_descent_one_step(alpha,e_grad,i,A,beta);
    }

    free(energy);
    free(energy_grad);
}

void task3(int argc, char *argv[], gsl_rng *r)
{
    //FILE output
    FILE *energy_alpha_t3 = fopen("op_text/energy_evol_alpha_t3.txt","w");
    //hard coded statistical inefficiency
    double n_s = 36;
    //limits
    int N = atoi(argv[2]);
    double d = atof(argv[3]);

    //energy array
    double *energy = (double *)malloc(N * sizeof(double));

    //initial states
    double r1[] = {0,0,0};
    double r2[] = {0,0,0};
    double alpha;

    mcmc_step iteration;

    //simulation parameters
    int i;
    int eqb_limit = 1E5;
    double e = 0;
    double e_err = 0;
    double sd = 0;
    double alpha_ll = 0.05;
    double alpha_ul = 0.25;

    //init energy arr
    for(i = 0; i < N; i++)
        energy[i] = 0;

    //iterating through different alpha values 
    for(alpha = alpha_ll; alpha <=  alpha_ul; alpha += 0.01)
    {
        //init random positions
        r2[0] = (100*d) * (gsl_rng_uniform(r) - 0.5);
        r2[1] = (100*d) * (gsl_rng_uniform(r) - 0.5);
        r2[2] = (100*d) * (gsl_rng_uniform(r) - 0.5);

        r2[0] = (100*d) * (gsl_rng_uniform(r) - 0.5);
        r2[1] = (100*d) * (gsl_rng_uniform(r) - 0.5);
        r2[2] = (100*d) * (gsl_rng_uniform(r) - 0.5);

        //energy calc using MCMC for alpha
        //eqb init run
        for(i = 0 ; i < eqb_limit; i++)
            iteration = mcmc_displace_all(r1,r2,d,alpha,r);
        
        //sampling run
        for(i = 0; i < N; i++)
        {
            iteration = mcmc_displace_all(r1,r2,d,alpha,r);
            energy[i] = iteration.energy;
        }

        //energy calculation
        e = average(energy,N);
        //error calculation
        sd = standard_deviation(energy,N);
        printf("%f\n",sd);
        e_err = n_s/N * (sd*sd);
        fprintf(energy_alpha_t3,"%f\t%f\t%.30f\n",alpha,e,e_err);
    }

    free(energy);
}

void task2_block_avg_inefficiency()
{
    //energy samples input
    FILE *ip_energy = fopen("op_text/mcmc_energy_eqb_t2.txt","r");
    //block inefficiency output
    FILE *op_blocks = fopen("op_text/block_avg_ineff_t2.txt","w");

    //n samples
    int N;
    fscanf(ip_energy,"%d\n",&N);

    double *data = (double *)malloc(N * sizeof(double));
    //read data
    int i;
    for(i = 0; i < N ; i++)
        fscanf(ip_energy,"%lf\n",&data[i]);
    
    double ineff = 0;
    
    //block averaging method of inefficiency calculation
    for(int b = 2; b < 1000; b++)
    {
        ineff = block_average(data, N, b);
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
    int N;
    fscanf(ip_energy,"%d\n",&N);

    double *data = (double *)malloc(N * sizeof(double));
    //read data
    int i;
    for(i = 0; i < N ; i++)
        fscanf(ip_energy,"%lf\n",&data[i]);
    
    //mean bias
    double f_mean = average(data,N);
    for(i = 0; i < N ; i++)
        data[i] -= f_mean;

    //autocorrelation
    for(i = 1; i < 200; i++)
        fprintf(op_autocorr,"%f\n",autocorrelation(data,N,i));

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
    else if(choice == 6)
        task3(argc, argv, r);
    else if(choice == 7)
        task4(argc, argv, r);
    else if(choice == 8)
        task5(argc, argv, r);

    gsl_rng_free(r);

    return 0;
}
