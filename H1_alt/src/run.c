#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "vv.h"
#include "eqb.h"

#define N 4
#define N_ATOMS 256
#define AL_MASS 0.0027965852524259004
#define GPA_SCALE 6.2415076486555486E-03

void task4(int argc, char *argv[])
{
    FILE *out = fopen("op_text/eqb_evol_t4.txt","w");
    FILE *ip_lattice = fopen("op_text/energy_evol_lattice.xyz","r");
    
    double **positions = create_2D_array(N_ATOMS, 3);
    double **velocities = create_2D_array(N_ATOMS, 3);
    double **forces = create_2D_array(N_ATOMS, 3);

    double kin_e = 0;
    double virial = 0;
    
    double lp_0k = 0;

    char symbol[256];
    read_xyz(ip_lattice, symbol, positions, velocities, &lp_0k);

    int i;
    for(i = 0 ; i < N_ATOMS; i++)
    {
        forces[i][0] = 0;
        forces[i][1] = 0;
        forces[i][2] = 0;
    }

    //simulation settings
    double timestep = atof(argv[2]);
    double max_time = atof(argv[3]);

    double T_eq_temp = 1173.15;
    double T_eq = 973.15;
    double T_const = 500 * timestep;
    double P_eq = 1E-4 * GPA_SCALE;
    double P_const = 500 * timestep;

    double lp = lp_0k; //initial

    //eqb1
    for(double time = 0; time < max_time; time += timestep)
    {
        vv_one_step(forces, velocities, positions, NULL, &virial, timestep, N*lp, N_ATOMS, AL_MASS);
        kin_e = calculate_kinetic_energy(velocities, AL_MASS, N_ATOMS);

        eqb_temp_1_step(velocities, T_eq_temp, T_const, timestep, kin_e, AL_MASS, N_ATOMS);
        eqb_press_1_step(velocities, virial, P_eq, P_const, timestep, kin_e, AL_MASS, N_ATOMS, &lp);

        fprintf(out,"%f\t%f\n",inst_temp_calc(kin_e, N_ATOMS), inst_press_calc(kin_e, virial, N*N*N*lp*lp*lp));
    }

    //eqb2
    for(double time = 0; time < max_time; time += timestep)
    {
        vv_one_step(forces, velocities, positions, NULL, &virial, timestep, N*lp, N_ATOMS, AL_MASS);
        kin_e = calculate_kinetic_energy(velocities, AL_MASS, N_ATOMS);

        eqb_temp_1_step(velocities, T_eq, T_const, timestep, kin_e, AL_MASS, N_ATOMS);
        eqb_press_1_step(velocities, virial, P_eq, P_const, timestep, kin_e, AL_MASS, N_ATOMS, &lp);

        fprintf(out,"%f\t%f\n",inst_temp_calc(kin_e, N_ATOMS), inst_press_calc(kin_e, virial, N*N*N*lp*lp*lp));
    }

    //production run
    for(double time = 0; time < max_time * 2; time += timestep)
    {

        vv_one_step(forces, velocities, positions, NULL, &virial, timestep, N*lp, N_ATOMS, AL_MASS);
        kin_e = calculate_kinetic_energy(velocities, AL_MASS, N_ATOMS);

        fprintf(out,"%f\t%f\n",inst_temp_calc(kin_e, N_ATOMS), inst_press_calc(kin_e, virial, N*N*N*lp*lp*lp));
    }

    destroy_2D_array(positions, N_ATOMS);
    destroy_2D_array(velocities, N_ATOMS);
    destroy_2D_array(forces, N_ATOMS);

}

void task3(int argc, char *argv[])
{
    FILE *out = fopen("op_text/eqb_evol.txt","w");
    FILE *ip_lattice = fopen("op_text/energy_evol_lattice.xyz","r");
    
    double **positions = create_2D_array(N_ATOMS, 3);
    double **velocities = create_2D_array(N_ATOMS, 3);
    double **forces = create_2D_array(N_ATOMS, 3);

    double kin_e = 0;
    double virial = 0;
    
    double lp_0k = 0;

    char symbol[256];
    read_xyz(ip_lattice, symbol, positions, velocities, &lp_0k);

    int i;
    for(i = 0 ; i < N_ATOMS; i++)
    {
        forces[i][0] = 0;
        forces[i][1] = 0;
        forces[i][2] = 0;
    }

    //simulation settings
    double timestep = atof(argv[2]);
    double max_time = atof(argv[3]);

    double T_eq = 773.15;
    double T_const = 500 * timestep;
    double P_eq = 1E-4 * GPA_SCALE;
    double P_const = 500 * timestep;

    double lp = lp_0k; //initial

    for(double time = 0; time < max_time; time += timestep)
    {
        vv_one_step(forces, velocities, positions, NULL, &virial, timestep, N*lp, N_ATOMS, AL_MASS);
        kin_e = calculate_kinetic_energy(velocities, AL_MASS, N_ATOMS);

        eqb_temp_1_step(velocities, T_eq, T_const, timestep, kin_e, AL_MASS, N_ATOMS);
        eqb_press_1_step(velocities, virial, P_eq, P_const, timestep, kin_e, AL_MASS, N_ATOMS, &lp);

        fprintf(out,"%f\t%f\n",inst_temp_calc(kin_e, N_ATOMS), inst_press_calc(kin_e, virial, N*N*N*lp*lp*lp));
    }

    //production run
    for(double time = 0; time < max_time * 2; time += timestep)
    {

        vv_one_step(forces, velocities, positions, NULL, &virial, timestep, N*lp, N_ATOMS, AL_MASS);
        kin_e = calculate_kinetic_energy(velocities, AL_MASS, N_ATOMS);

        fprintf(out,"%f\t%f\n",inst_temp_calc(kin_e, N_ATOMS), inst_press_calc(kin_e, virial, N*N*N*lp*lp*lp));
    }

    destroy_2D_array(positions, N_ATOMS);
    destroy_2D_array(velocities, N_ATOMS);
    destroy_2D_array(forces, N_ATOMS);

}

void task2(int argc, char *argv[], gsl_rng *r)
{
    FILE *lp_0_k = fopen("op_text/lp_0_k.txt","r");
    FILE *out = fopen("op_text/energy_evol.txt","w");
    FILE *op_lattice = fopen("op_text/energy_evol_lattice.xyz","w");

    double **positions = create_2D_array(N_ATOMS, 3);
    double **velocities = create_2D_array(N_ATOMS, 3);
    double **forces = create_2D_array(N_ATOMS, 3);
    double pot_e = 0;
    double kin_e = 0;
    double lp_0k = 0;

    fscanf(lp_0_k,"%lf",&lp_0k);//read initial lp from previous task

    init_fcc(positions, N, lp_0k);

    int i;

    //randomize
    for(i = 0; i < N_ATOMS; i++)
    {
        positions[i][0] += ((gsl_rng_uniform(r) - 0.5) * 0.13 * lp_0k);
        positions[i][1] += ((gsl_rng_uniform(r) - 0.5) * 0.13 * lp_0k);
        positions[i][2] += ((gsl_rng_uniform(r) - 0.5) * 0.13 * lp_0k);

        forces[i][0] = 0;
        forces[i][1] = 0;
        forces[i][2] = 0;

        velocities[i][0] = 0;
        velocities[i][1] = 0;
        velocities[i][2] = 0;
    }

    //simulation settings
    double timestep = atof(argv[2]);
    double max_time = atof(argv[3]);

    for(double time = 0; time < max_time; time += timestep)
    {
        vv_one_step(forces, velocities, positions, &pot_e, NULL, timestep, N*lp_0k, N_ATOMS, AL_MASS);
        kin_e = calculate_kinetic_energy(velocities, AL_MASS, N_ATOMS);
        
        fprintf(out,"%f\t%f\n",kin_e,pot_e);
    }

    write_xyz(op_lattice, "Al", positions, velocities, lp_0k, N_ATOMS);

    destroy_2D_array(forces, N_ATOMS);
    destroy_2D_array(velocities, N_ATOMS);
    destroy_2D_array(positions, N_ATOMS);
}

void task1()
{
    FILE * pot_e_lp = fopen("op_text/pot_e_lp.txt","w");

    double **positions = create_2D_array(N_ATOMS, 3);

    double lp = 0;
    double lp_ll = 4;
    double lp_ul = 4.08;
    double lp_step = 0.01;

    double pot_e = 0;

    for(lp = lp_ll; lp < lp_ul; lp += lp_step)
    {
        init_fcc(positions , N , lp);
        pot_e = get_energy_AL(positions , N * lp , N_ATOMS);
        fprintf(pot_e_lp, "%f\t%f\n", lp, pot_e);
    }

    destroy_2D_array(positions, N_ATOMS);
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
        task1();
    else if(choice == 2)
        task2(argc, argv, r);
    else if(choice == 3)
        task3(argc, argv);
    else if(choice == 4)
        task4(argc, argv);

    gsl_rng_free(r);

    return 0;
}