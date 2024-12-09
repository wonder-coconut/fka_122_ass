#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "tools.h"

#define N_ATOMS 256
#define LP 4.236855
#define T 1000

void static_factor(int argc, char *argv[])
{
    FILE *pos = fopen("op_text/pos_evol_t4_temp.txt","r");

    FILE *sq = fopen("op_text/sq_t7.txt","w");
    FILE *qgrid = fopen("op_text/qgrid_t7.txt","w");

    int n_max = atoi(argv[2]);
    int q_dim = pow(n_max, 3);

    double L = 4 * LP;

    double **positions = create_2D_array(N_ATOMS,3);
    double **q_grid = create_2D_array(q_dim,3);    

    //generate grid
    int i,j,k,t,c,q;
    i = j = k = t = c = q = 0;
    for(i = 0; i < n_max; i++)
    {
        for(j = 0; j < n_max; j++)
        {
            for(k = 0; k < n_max; k++)
            {
                q_grid[c][0] = i*2*M_PI/L;
                q_grid[c][1] = j*2*M_PI/L;
                q_grid[c][2] = k*2*M_PI/L;
                fprintf(qgrid,"%f\t%f\t%f\n",q_grid[c][0],q_grid[c][1],q_grid[c][2]);
                c++;
            }
        }
    }
    double complex s_q = 0;
    double dot = 0;
    double complex exp_arg = 0;
    double complex temp = 0;

    c = 0;
    //every grid point in q_grid
    for(q = 0; q < q_dim; q++)
    {
        s_q = 0;
        //every timestep until t=1000
        for(t = 0; t < T; t++)
        {
            temp = 0;
            //every position at timestep t
            for(i = 0; i < N_ATOMS; i++)
            {
                fscanf(pos,"%lf\t%lf\t%lf\t",&positions[i][0],&positions[i][1],&positions[i][2]);
                dot = dot_product(q_grid[q],positions[i],3);
                exp_arg = I*dot;
                temp += cexp(exp_arg);//compute exponent
            }
            s_q += cpow(cabs(temp),2); //summation into sq
        }
        fprintf(sq,"%f\n",creal(s_q)/(N_ATOMS*T)); //ensemble average and normalization
    }

    destroy_2D_array(positions,N_ATOMS);
    destroy_2D_array(q_grid,N_ATOMS);
}