#include <stdio.h>
#include <stdlib.h>

#include "tools.h"
#include "velocity_verlet.h"

int main()
{   
    int n_samples = 10000;

    //constants init
    double C_mass = 0.0012437130574206093;
    double O_mass = 0.0016582840765608123;
    double kappa = 99.86414519137223;

    //initial conditions
    
    //positions (in angstroms)
    double positions[3] = {0.01,0.0005,-0.005};
    double **position_samples = create_2D_array(3,n_samples);

    //velocities
    double velocities[3] = {0,0,0};

    //masses (O-C-O chain)
    double masses[3] = {O_mass, C_mass, O_mass};

    //accelerations (init)
    double accelerations[3] = {0,0,0};

    calculate_acceleration(accelerations, positions, masses, kappa);

    double timestep = 0.00025;

    double pot_e, kin_e;
    pot_e = kin_e = 0;
    //pot_e = calculate_potential_energy(positions,kappa);

    FILE *fptr,*fptr2;
    fptr = fopen("op_text/output.txt","w");
    fptr2 = fopen("op_text/power_spectrum.txt","w");

    fprintf(fptr,"%d\t%f\n",n_samples,timestep);

    fprintf(fptr,"x0\t\t\tx1\t\t\tx2\t\t\tpot_e\t\t\tkin_e");
    fprintf(fptr,"\n");

    fprintf(fptr,"%f\t%f\t%f\t",positions[0],positions[1],positions[2]);
    fprintf(fptr,"%f\t%f",pot_e,kin_e);
    fprintf(fptr,"\n");

    int i;

    for(i = 0; i < n_samples; i++)
    {
        velocity_verlet_one_step(accelerations,positions,velocities,masses,kappa,timestep);

        position_samples[0][i] = positions[0];
        position_samples[1][i] = positions[1];
        position_samples[2][i] = positions[2];
        
        pot_e = calculate_potential_energy(positions,kappa);
        kin_e = calculate_kinetic_energy(velocities,masses);

        fprintf(fptr,"%f\t%f\t%f\t",positions[0],positions[1],positions[2]);
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

    int idx0,idx1,idx2;
    idx0 = idx1 = idx2 = 0;

    fprintf(fptr2,"%d\t%f\n",n_samples,timestep);

    for(i = 0; i < n_samples; i++)
        fprintf(fptr2,"%f\t%f\t%f\t%f\n",f[i],P_n[0][i] * 1000000,P_n[1][i] * 1000000,P_n[2][i] * 1000000);

    for(i = 0; i < n_samples; i++)
    {
        idx0 = P_n[0][i] > P_n[0][idx0] ? i : idx0;
        idx1 = P_n[1][i] > P_n[1][idx1] ? i : idx1;
        idx2 = P_n[2][i] > P_n[2][idx2] ? i : idx2;
    }

    //printf("%f\n%f\n%f\n",f[idx0],f[idx1],f[idx2]);

    destroy_2D_array(position_samples);
    destroy_2D_array(P_n);
    free(f);

    return 0;
}

