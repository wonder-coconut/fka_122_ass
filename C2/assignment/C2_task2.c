#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

/* Compile:
gcc C2_task2.c -o C2 -I${loc}/include -L${loc}/lib -Wl,-rpath,${loc}/lib -lgsl -lgslcblas

Execute:
./C2

Input parameters: <array dimension>
*/

int main(int argc, char *argv[])
{    
    //output file
    FILE* op = fopen("rng_dist.txt","w");
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

    //dynamic array
    float *arr = (float*)malloc(dim * sizeof(float));

    //rng init
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    time_t seed = time(NULL);

    gsl_rng_set(r, seed);

    int i;
    for(i = 0; i < dim; i++)
    {
        arr[i] = gsl_rng_uniform(r);
        if(i == dim - 1)
            fprintf(op,"%f",arr[i]);
        else
            fprintf(op,"%f\n",arr[i]);
    }
    return 0;
}