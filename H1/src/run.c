#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include<string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "velocity_verlet.h"
#include "kinetic_energy.c"

void task3()
{
    
}

void task2()
{
    FILE* f_in = fopen("op_text/lp_0k.txt","r");

    //gsl setup
    const gsl_rng_type * T;
    gsl_rng * r;

    // Initialize RNG
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Set seed
    gsl_rng_set(r, 299792458);

    //lattice parameters init
    int N = 4;
    int fcc_n_atoms = 4;
    int n_atoms = N*N*N*fcc_n_atoms;
    double al_mass_num = 26.981539;
    int asu_scale = 9649;
    double al_asu = al_mass_num/asu_scale;

    //input lp from file
    double lp = 0;
    fscanf(f_in,"%lf",&lp);
    double bias = -0.5; //rng scale for +- 6.5%
    double scale = 0.13 * lp;
    
    double pot_e = 0;
    double kin_e = 0;

    double **positions = create_2D_array(n_atoms, 3); //atomic properties storage
    double **deviations = create_2D_array(n_atoms, 3);
    double **velocities = create_2D_array(n_atoms, 3);
    double **accelerations = create_2D_array(n_atoms, 3);
    double **forces = create_2D_array(n_atoms, 3);

    int i,j;

    //init to 0
    for(i = 0; i < n_atoms; i++)
    {
        for(j = 0; j < 3; j++)
        {
            positions[i][j] = 0;
            deviations[i][j] = 0;
            velocities[i][j] = 0;
            forces[i][j] = 0;
        }
    }

    //lattice init
    init_fcc(positions,N,lp);
    
    //random deviations
    for(i = 0; i < n_atoms; i++)
    {
        positions[i][0] += (gsl_rng_uniform(r) - bias) * scale;
        positions[i][1] += (gsl_rng_uniform(r) - bias) * scale;
        positions[i][2] += (gsl_rng_uniform(r) - bias) * scale;
    }
    //simulation setup

    double timestep = 0.0001;
    double t_max = 10;
    double max_iterations = t_max/timestep;


    char filename[256]; 
	sprintf(filename, "op_text/energy_evol_%f.txt", timestep);
    printf("%s\n",filename);

    FILE* f_out = fopen(filename,"w");

    fprintf(f_out,"%f\n",timestep);

    pot_e = get_energy_AL(positions, N * lp, n_atoms);
    kin_e = get_kinetic_energy_al(velocities, al_asu, n_atoms);
    fprintf(f_out,"%f\t%f\n",pot_e,kin_e);

    for(i = 0; i < max_iterations; i++)
    {
        velocity_verlet_one_step_lattice(forces, velocities, positions,timestep, N, lp, n_atoms, al_asu);
        pot_e = get_energy_AL(positions, N * lp, n_atoms);
        kin_e = get_kinetic_energy_al(velocities, al_asu, n_atoms);
        fprintf(f_out,"%f\t%f\n",pot_e,kin_e);
    }

    //free dynamic memory
    destroy_2D_array(positions, n_atoms);
    destroy_2D_array(deviations, n_atoms);
    destroy_2D_array(velocities, n_atoms);
    destroy_2D_array(accelerations, n_atoms);
    destroy_2D_array(forces, n_atoms);

    gsl_rng_free(r);

}

void task1()
{
    FILE* fptr = fopen("op_text/pot_e_dist.txt","w");
    
    //lattice parameters init
    int N = 4;
    int fcc_n_atoms = 4;
    int n_atoms = N*N*N*fcc_n_atoms;

    double pot_e = 0;

    double lp_ll = 4;
    double lp_ul = 4.1;
    double lp = 0;
    double lp_step = 0.01;
    double iterations = (lp_ul - lp_ll) / lp_step;
    
    double **positions = create_2D_array(n_atoms, 3);

    int i;
    //fcc lattices for varying lp
    for(i = 0; i < iterations; i++)
    {
        lp = lp_ll + i*lp_step;
        init_fcc(positions, N, lp);

        pot_e = get_energy_AL(positions, N * lp, n_atoms);
        fprintf(fptr,"%f\t%f\n",lp,pot_e);
    }

    destroy_2D_array(positions,n_atoms);
}

int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code
    int choice = atoi(argv[1]);
    if(choice == 1)
        task1();
    else if(choice == 2)
        task2(argc, argv);
    return 0;
}
