#include <math.h> //add symbolic link lm
#include <stdio.h>

double scalar_product(double *a, double *b, int dim)
{
    double res = 0;

    int i;
    for(i = 0; i < dim; i++)
        res += a[i] * b[i];

    return res;
}

double coordinate_distance(double *a, double *b)
{
    //asuming 3D space
    int dim = 3;
    double distance = 0;

    int i;
    for(i = 0; i < dim; i++)
    {
        printf("%f   %f :\t",a[i],b[i]);
        printf("%f\n",(a[i] - b[i]) * (a[i] - b[i]));
        distance += (a[i] - b[i]) * (a[i] - b[i]);
    }
    distance = sqrt(distance);

    return distance;
}