#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#include "tools.h"
#include "result_t.h"
#include "bd3.h"

#define KB 1.3806503E-08 //micrometer^2 microgram/milisecond^2 K
#define MASS 3.0134022351607534e-05 //microgram
#define V_THERM 0.36888544303442244 //micrometer/milisecond
#define W0 19.477874452256717 //rad/milisecond
#define VOLUME 11.371329189285863 //micrometer^3
#define T 297 //kelvin

void t_corr(FILE *fin,FILE *fout, double dt)
{

    //signal length
    int n = 0;

    double temp = 0;

    //read file
    fscanf(fin,"%lf\n",&temp);
    fscanf(fin,"%d\n",&n);

    //signal arrays
    double *v = (double*)malloc(sizeof(double) * n);


    //read file
    int i;
    for(i = 0; i < n; i++)
        fscanf(fin,"%lf\t%lf\t%lf\n",&temp,&temp,&v[i]);

    //sim params
    double corr_dt = 2*dt;
    double sum = 0;
    double M = 1000;
    double l = 0;
    double m = 0;
    double a_m = 0;
    double a_m_l = 0;
    double time = 0;
    int idx = 0;
    
    for(l = 0; l < M; l++)
    {
        sum = 0;
        for(m = 0; m < M - l; m++)
        {
            time = m*corr_dt;
            idx = (int)time/dt;

            a_m = v[idx];
            
            time = (m+l)*corr_dt;
            idx = (int)time/dt;

            a_m_l = v[idx];

            sum+= a_m * a_m_l;
        }
        fprintf(fout,"%f\n",sum/(M - l));
    }

    free(v);
}

void task5()
{
    //file in
    FILE *traj_in_low = fopen("op_text/traj_out_low_0.001.txt","r");
    //FILE *traj_in_high = fopen("op_text/traj_out_high_0.001.txt","r");
    
    //file out
    FILE *t_corr_low = fopen("op_text/t_corr_low.txt","w");
    //FILE *t_corr_high = fopen("op_text/t_corr_high.txt","w");
    //signal length

    t_corr(traj_in_low,t_corr_low,0.001);
}

void fft(FILE *fin, FILE *fout, double dt)
{
    //signal length
    int n = 0;
    double temp = 0;

    //read file
    int i;
    fscanf(fin,"%lf\n",&temp);
    fscanf(fin,"%d\n",&n);

    //signal arrays
    double *v = (double*)malloc(sizeof(double) * n);
    double *p = (double*)malloc(sizeof(double) * n);
    double *f = (double*)malloc(sizeof(double) * n);

    for(i = 0; i < n; i++)
        fscanf(fin,"%lf\t%lf\t%lf\n",&temp,&temp,&v[i]);

    //fft
    fft_freq(f,n,dt);
    powerspectrum(p,v,n,dt);

    //write to file
    for(i = 0; i < n; i++)
        fprintf(fout,"%f\t%f\n",f[i],p[i]);

    //free dynamic array
    free(v);
    free(p);
    free(f);
}

void task4()
{
    //file in
    FILE *traj_in_low = fopen("op_text/traj_out_low_0.001.txt","r");
    FILE *traj_in_high = fopen("op_text/traj_out_high_0.001.txt","r");
    ////file out
    FILE *pspec_low_5 = fopen("op_text/pspec_low_0.05.txt","w");
    FILE *pspec_high_5 = fopen("op_text/pspec_high_0.05.txt","w");
    FILE *pspec_low_25 = fopen("op_text/pspec_low_0.025.txt","w");
    FILE *pspec_high_25 = fopen("op_text/pspec_high_0.025.txt","w");
    
    double dt_fft = 0;

    //case 1 for dt = 0.05
    dt_fft = 0.05;
    fft(traj_in_low,pspec_low_5,dt_fft);

    traj_in_low = fopen("op_text/traj_out_low_0.001.txt","r");
    //case 2 for dt = 0.025
    dt_fft = 0.025;
    fft(traj_in_low,pspec_low_25,dt_fft);

    //case 3 for dt = 0.05
    dt_fft = 0.05;
    fft(traj_in_high,pspec_high_5,dt_fft);

    traj_in_high = fopen("op_text/traj_out_high_0.001.txt","r");
    //case 4 for dt = 0.025
    dt_fft = 0.025;
    fft(traj_in_high,pspec_high_25,dt_fft);
    
}

void run_sim(FILE* f_out, result_t res , int eqb_iterations, int sim_iterations, double dt, double eta, gsl_rng *k)
{
    int i;
    double t;
    //eqb run
    for(i = 0; i < eqb_iterations; i++)
        res = BD3(res.position, res.velocity, W0, dt, eta, KB, MASS, T, k);

    //sim run
    t = 0;
    for(i = 0; i < sim_iterations; i++)
    {
        t += dt; 
        res = BD3(res.position, res.velocity, W0, dt, eta, KB, MASS, T, k);
        fprintf(f_out,"%f\t%f\t%f\n",t,res.position,res.velocity);
    }
}

void task3(gsl_rng *k)
{
    result_t res;

    FILE *traj_out_low_1 = fopen("op_text/traj_out_low_0.001.txt","w");
    FILE *traj_out_low_5 = fopen("op_text/traj_out_low_0.005.txt","w");
    FILE *traj_out_high_1 = fopen("op_text/traj_out_high_0.001.txt","w");
    FILE *traj_out_high_5 = fopen("op_text/traj_out_high_0.005.txt","w");

    //all time units are chosen to be in milliseconds
    double tau_low = 147.3 * 1E-3;
    double tau_high = 48.5 * 1E-3;

    double eta_low = 1/tau_low;
    double eta_high = 1/tau_high;

    double dt_2 = 0.005;
    double dt_1 = 0.001;

    //arbitrary initial conditions (distance is in micrometers)
    res.position = 0.001;
    res.velocity = 0.001;
    //simulation settings
    double dt;
    double eqb_time_low = tau_low * 1000;
    double eqb_time_high = tau_high * 1000;
    double sim_time = 10;
    int eqb_iterations;
    int sim_iterations;
    
    
    //sim 1 : P = 99.8 KPa, dt = 0.001
    dt = dt_1;
    eqb_iterations = (int)eqb_time_low/dt;
    sim_iterations = (int)sim_time/dt;
    fprintf(traj_out_low_1,"99.8\n%d\n",sim_iterations);
    run_sim(traj_out_low_1,res,eqb_iterations,sim_iterations,dt,eta_low, k);

    //sim 2 : P = 99.8 KPa, dt = 0.001
    dt = dt_2;
    eqb_iterations = (int)eqb_time_low/dt;
    sim_iterations = (int)sim_time/dt;
    fprintf(traj_out_low_5,"99.8\n%d\n",sim_iterations);
    run_sim(traj_out_low_5,res,eqb_iterations,sim_iterations,dt,eta_low, k);

    //sim 3 : P = 2.75 KPa, dt = 0.001
    dt = dt_1;
    eqb_iterations = (int)eqb_time_high/dt;
    sim_iterations = (int)sim_time/dt;
    fprintf(traj_out_high_1,"2.75\n%d\n",sim_iterations);
    run_sim(traj_out_high_1,res,eqb_iterations,sim_iterations,dt,eta_high, k);

    //sim 4 : P = 2.75 KPa, dt = 0.005
    dt = dt_2;
    eqb_iterations = (int)eqb_time_high/dt;
    sim_iterations = (int)sim_time/dt;
    fprintf(traj_out_high_5,"2.75\n%d\n",sim_iterations);
    run_sim(traj_out_high_5,res,eqb_iterations,sim_iterations,dt,eta_high, k);
    
}

int main(int argc, char *argv[])
{
    int epochs = (unsigned long)time(NULL);
    //rng setup
    const gsl_rng_type * t;
    gsl_rng * r;

    // init gsl
    gsl_rng_env_setup();
    t = gsl_rng_default;
    r = gsl_rng_alloc(t);

    // Set seed
    gsl_rng_set(r, epochs);

    //task choice
    int choice = atoi(argv[1]);

    if(choice == 1)
        task3(r);
    else if(choice == 2)
        task4();
    else if(choice == 3)
        task5();

    gsl_rng_free(r);
}