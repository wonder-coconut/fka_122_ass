#include <stdio.h>

//importing math header for square root, note: add corresponding -lm flag in makefile(s)
#include <math.h>

#include "vector.h"

void
elementwise_addition(
		     double *res,
		     double *v1,
		     double *v2,
		     unsigned int len
	            )
{
    for(int i = 0; i < len; ++i){
	res[i] = v1[i] + v2[i];
    }
}

void
constant_multiplication(
			   double *res,
			   double *v1,
			   double a,
			   unsigned int len
	                  )
{
    for(int i = 0; i < len; ++i){
	res[i] = v1[i] * a;
    }
}

double
dot_product(
	    double *v1,
	    double *v2,
	    unsigned int len
	   )
{
    double result = 0;
    for(int i = 0; i < len; ++i){
	result += v1[i] * v2[i];
    }
    return result;
}

//added function to calculate the L2 norm of a vector
double norm(double *v1, unsigned int len)
{
	double res = 0.0;

	for(int i = 0; i < len ; i++)
		res += v1[i] * v1[i];
    
	return sqrt(res);
}