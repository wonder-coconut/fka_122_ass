#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#include "linalg.h"

/* compile using :
gcc C2_task1.c linalg.c -o C2 -I${loc}/include -L${loc}/lib -Wl,-rpath,${loc}/lib -lgsl -lgslcblas -lm

Execute using :
./C2

Input parameters: <array dimension>
*/
double scalar_driver(int);
double distance_driver(int);

int main(int argc, char *argv[])
{
    
    //sanity check
    if(argc == 1)
    {
        printf("please enter a value.\n");
        exit(1);
    }

    //dimension of random arrays input from arguments
    int dim = atoi(argv[1]);

    //sanity check
    if(dim <= 0)
    {
        printf("invalid input, must be a positive integer\n");
        exit(1);
    }

    //task1.1
    printf("\nScalar Product:\t%.5f\n-*-*-*-*-*-*-*-*-*-*-*-*-\n",scalar_driver(dim));

    //task1.2
    printf("\nEuclidean Distance:\t%f\n",distance_driver(dim));
}

double scalar_driver(int dim)
{
    //dynamic array allocation
    double *ar1 = (double*)malloc(dim * sizeof(double));
    double *ar2 = (double*)malloc(dim * sizeof(double));

    //rng to initialize arrays
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    time_t seed = time(NULL);

    gsl_rng_set(r, seed);

    printf("random arrays:\n");

    int i;
    for(i = 0; i < dim; i++)
    {
        ar1[i] = gsl_rng_uniform(r);
        ar2[i] = gsl_rng_uniform(r);
        printf("%.3f\t%.3f\n",ar1[i],ar2[i]);
    }

    //scalar product
    double res = scalar_product(ar1,ar2, dim);

    free(ar1);
    free(ar2);

    return res;
}

double distance_driver(int dim)
{
    //file creation for output
    FILE* op = fopen("matrix.csv","w");

    int ncol = 3;

    //dynamic 2D array allocation
    double *arr = (double *)malloc(dim * ncol * sizeof(double));

    //coordinate points pair allocation
    double *p1 = (double *)malloc(ncol * sizeof(double));
    double *p2 = (double *)malloc(ncol * sizeof(double));

    //rng to initialize arrays
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    time_t seed = time(NULL);

    gsl_rng_set(r, seed);

    int i,j;
    for(i = 0; i < dim; i++)
    {   
        printf("%d) ",i);
        for(j = 0; j < ncol; j++)
        {
            arr[i * ncol + j] = gsl_rng_uniform(r);
            printf("%.3f  ", arr[i * ncol + j]);

            //write to file
            if(j == ncol - 1)
                fprintf(op, "%.3f", arr[i * ncol + j]);
            else
                fprintf(op, "%.3f,", arr[i * ncol + j]);

        }
        printf("\n");
        if(i != dim - 1)
            fprintf(op,"\n");
    }

    //two random points in dim range
    int i1 = (int)(gsl_rng_uniform(r) * dim);
    int i2 = (int)(gsl_rng_uniform(r) * dim);

    printf("\nRandom points:\t%d, %d\n",i1,i2);
    
    for(j = 0; j < ncol; j++)
    {
        p1[j] = arr[i1 * ncol + j];
        p2[j] = arr[i2 * ncol + j];
    }

    if(op == NULL)
    {
        printf("file creation failed\n");
        exit(1);
    }

    fclose(op);
    free(arr);

    double res = coordinate_distance(p1,p2);

    free(p1);
    free(p2);

    return res;
}