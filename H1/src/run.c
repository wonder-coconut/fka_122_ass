#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"

void task2(int argc, char *argv[])
{

    const gsl_rng_type * T;
    gsl_rng * r;

    // Initialize RNG
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Set seed
    gsl_rng_set(r, 299792458);
    
    FILE* fptr = fopen("op_text/temp.xyz","w");

    //initialize lattice constants
    int N = atoi(argv[2]);
    int fcc_n_atoms = 4;
    int n_atoms = N*N*N*fcc_n_atoms;

    //double lp = 4.03031139135711; //lattice parameter at 0K
    double lp = 4.053;
    double bias = -0.5; //rng scale for +- 6.5%
    double scale = 0.13 * lp;

    //double asu_scale = 9649; //converting aluminium to asu
    //double al_amu =  26.981539;
    //double al_asu = al_amu / asu_scale;

    double **positions = create_2D_array(n_atoms, 3); //atomic properties storage
    double **deviations = create_2D_array(n_atoms, 3);
    double **velocities = create_2D_array(n_atoms, 3);
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

    printf("fcc initialization\n");
    //lattice position init for lp
    init_fcc(positions,N,lp);
    printf("fcc init successful\n");

    //random deviations
    for(i = 0; i < n_atoms; i++)
    {
        positions[i][0] += (gsl_rng_uniform(r) - bias) * scale;
        positions[i][1] += (gsl_rng_uniform(r) - bias) * scale;
        positions[i][2] += (gsl_rng_uniform(r) - bias) * scale;
    }
    
    printf("force generation\n");
    get_forces_AL(forces, positions, lp, n_atoms);
    printf("force gen successful\n");
    //multiplication_with_constant(velocities[0], forces[0], al_asu, n_atoms);
    //multiplication_with_constant(velocities[1], forces[1], al_asu, n_atoms);
    //multiplication_with_constant(velocities[2], forces[2], al_asu, n_atoms);

    write_xyz(fptr, "Al", positions, velocities, lp, n_atoms);

    destroy_2D_array(positions, n_atoms);
    destroy_2D_array(deviations, n_atoms);
    destroy_2D_array(velocities, n_atoms);
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
