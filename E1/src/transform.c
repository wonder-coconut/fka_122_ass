#include "tools.h"
#include <math.h>
#include <stdlib.h>

/*
 * Implement the functions in this file.
 * Utilize the fact that the transformation can be written as
 * a matrix multiplication. To this end, you can use the following
 * functions from C4. No need to declare them, Yata already knows about them
 *

double **create_2D_array(unsigned int n, unsigned int m);

void destroy_2D_array(double **array, unsigned int n);

void matrix_vector_multiplication(double *result, double **A, double *b, unsigned int n, unsigned int m);

 **/

void gen_transform_mat(double **transform_mat, const unsigned int N)
{
    double coeff1 = M_PI / (N + 1);

    int i,j;

    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            transform_mat[i][j] = sin(coeff1 * (i + 1) * (j + 1));
}

/**
 * Transform to normal modes.
 *
 * Parameters
 * ----------
 *  positions - Vector of positions
 *  Q - Vector of Q-coordinates
 *  N - Number of particles (length of vectors)
 *
*/

void transform_to_normal_modes(double *positions, double **transform_mat, double *Q, const unsigned int N)
{
    double coeff2 = sqrt(2.0/ (N + 1));

    int i;
    
    matrix_vector_multiplication(Q, transform_mat, positions, N, N);

    for(i = 0; i < N; i++)
        Q[i] *= coeff2;
}



/**
 * Calculate the normal mode energies
 *
 * Parameters
 * ----------
 *  energies - Vector where energies will be written
 *  positions - Vector of positions
 *  velocities - Vector of positions
 *  N - Number of particles (length of vectors)
 *
*/
void calculate_normal_mode_energies(double *energies, double *positions, double **transform_mat, double *velocities, const unsigned int N)
{
    double *omega = (double *)malloc(sizeof(double) * N);
    double *P = (double *)malloc(sizeof(double) * N);
    double *Q = (double *)malloc(sizeof(double) * N);

    transform_to_normal_modes(positions,transform_mat,Q,N);
    transform_to_normal_modes(velocities,transform_mat,P,N);

    double coeff = M_PI/(2*(N+1));

    int i;
    for(i = 0; i < N; i++)
    {
        omega[i] = 2 * sin(coeff * (i+1.0));
        energies[i] = (P[i] * P[i] + omega[i]*omega[i] * Q[i]*Q[i]) / 2;
    }

    free(omega);
    free(P);
    free(Q);
}
