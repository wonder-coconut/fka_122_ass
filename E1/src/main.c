#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tools.h"
#include "velocity_verlet.h"
#include "transform.h"
#include "fput.h"


#define C_MASS 0.0012437130574206093
#define O_MASS 0.0016582840765608123
#define KAPPA 99.86414519137223

void part1()
{
    int n_samples = 10000;

    //constants init

    //initial conditions
    
    //positions (in angstroms)
    double positions[3] = {0.01,0.0005,-0.005};
    double **position_samples = create_2D_array(3,n_samples);

    //velocities
    double velocities[3] = {0,0,0};

    //masses (O-C-O chain)
    double masses[3] = {O_MASS, C_MASS, O_MASS};

    //accelerations (init)
    double accelerations[3] = {0,0,0};

    calculate_acceleration(accelerations, positions, masses, KAPPA);

    double timestep = 0.00025;

    double pot_e, kin_e;
    pot_e = kin_e = 0;
    //pot_e = calculate_potential_energy(positions,kappa);

    FILE *fptr,*fptr2;
    fptr = fopen("op_text/output.txt","w"); //position and energies
    fptr2 = fopen("op_text/power_spectrum.txt","w"); //power spectrum

    fprintf(fptr,"%d\t%f\n",n_samples,timestep);

    fprintf(fptr,"x0\t\t\tx1\t\t\tx2\t\t\tpot_e\t\t\tkin_e");
    fprintf(fptr,"\n");

    fprintf(fptr,"%f\t%f\t%f\t",positions[0],positions[1],positions[2]);
    fprintf(fptr,"%f\t%f",pot_e,kin_e);
    fprintf(fptr,"\n");

    int i;

    for(i = 0; i < n_samples; i++)
    {
        velocity_verlet_one_step(accelerations,positions,velocities,masses,KAPPA,timestep); //one iteration

        position_samples[0][i] = positions[0];
        position_samples[1][i] = positions[1];
        position_samples[2][i] = positions[2];
        
        pot_e = calculate_potential_energy(positions,KAPPA);
        kin_e = calculate_kinetic_energy(velocities,masses);

        fprintf(fptr,"%f\t%f\t%f\t",positions[0],positions[1],positions[2]); // write to file
        fprintf(fptr,"%f\t%f",pot_e,kin_e);
        fprintf(fptr,"\n");
    }

    //fourier transform
    double *f = (double *)malloc(n_samples * sizeof(double));
    double **P_n = create_2D_array(3,n_samples);

    fft_freq(f,n_samples,timestep);

    powerspectrum(P_n[0],position_samples[0],n_samples,timestep);
    powerspectrum(P_n[1],position_samples[1],n_samples,timestep);
    powerspectrum(P_n[2],position_samples[2],n_samples,timestep);

    fprintf(fptr2,"%d\t%f\n",n_samples,timestep);

    for(i = 0; i < n_samples; i++)
        //scaling up for better visibility
        fprintf(fptr2,"%f\t%f\t%f\t%f\n",f[i],P_n[0][i] * 1000000,P_n[1][i] * 1000000,P_n[2][i] * 1000000);

    destroy_2D_array(position_samples);
    destroy_2D_array(P_n);
    free(f);

    fclose(fptr);
    fclose(fptr2);
}

void part2()
{
    const unsigned int N = 32; //number of particles
    double E_0 = N; //initial energy

    double *positions = (double *)malloc(sizeof(double) * N);
    double *velocities = (double *)malloc(sizeof(double) * N);
    double *accelerations = (double *)malloc(sizeof(double) * N);

    double *P_k = (double *)malloc(sizeof(double) * N); //normal mode momentums/velocities
    double *Q_k = (double *)malloc(sizeof(double) * N); //normal mode positions
    double *E_k = (double *)malloc(sizeof(double) * N); //normal mode energies

    double **transform_mat = create_2D_array(N,N);

    FILE* fptr = fopen("op_text/fput_energies.txt","w");

    int i;
    //initialization

    for(i = 0; i < N; i++)
        P_k[i] = Q_k[i] = E_k[i] = positions[i] = velocities[i] = accelerations[i] = 0;

    gen_transform_mat(transform_mat, N);

    E_k[0] = E_0; //initial energy conditions
    P_k[0] = sqrt(2 * E_0);

    transform_to_normal_modes(Q_k, transform_mat, positions, N);
    transform_to_normal_modes(P_k, transform_mat, velocities, N);

    //sim setup
    float timestep = 0.1;
    int t_max = 25000;
    float i_lim = t_max/timestep;

    fprintf(fptr,"%f\n",timestep);

    for(i = 0; i < i_lim; i++)
    {
        velocity_verlet_one_step_fput(accelerations, positions, velocities, 0, timestep, N);
        calculate_normal_mode_energies(E_k,positions,transform_mat,velocities,N);
        fprintf(fptr,"%f\t%f\t%f\t%f\t%f\n",E_k[0],E_k[1],E_k[2],E_k[3],E_k[4]);
    }

    destroy_2D_array(transform_mat);

    free(positions);
    free(velocities);
    free(accelerations);
    
    free(P_k);
    free(Q_k);
    free(E_k);

    fclose(fptr);
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

