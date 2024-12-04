#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tools.h"

double variance(double *v, unsigned int len)
{
    double v_mean = average(v,len);

    int i;
    double temp = 0;

    for(i = 0; i < len; i++)
        temp += v[i] * v[i];

    double temp_avg = temp/len;
    return temp_avg - (v_mean * v_mean);
    
}
/*
 * The following functions are to be used:
*/
/*
double
average(
        double *v,
        unsigned int len
       );
double
variance(
         double *v,
         unsigned int len
        );
*/


/* ***************************************
*
* Autocorrelation
*
* Parameters
* ----------
*  data - the raw data from which the
*         autocorrelation should be
*         calculated
*  data_len     - the len of data
*  time_lag_ind - the lag (index) at 
*                 which the autocorrelation
*                 should be calculated
*
* Returns
* -------
*  The autocorrelation at a specific time
*  lag
*
* ***************************************/
double autocorrelation(
			   double *data,
			   int data_len,
			   int time_lag_ind)
{
    double ret = 0;
    
    double variance_denom = variance(data, data_len);

    double *data_lagged = (double *)malloc((data_len - time_lag_ind) * sizeof(double));
    double *temp1 = (double *)malloc((data_len - time_lag_ind) * sizeof(double));
    
    int i;

    for(i = 0; i < data_len - time_lag_ind; i++)
        data_lagged[i] = data[i + time_lag_ind];
    
    for(i = 0; i < data_len; i++)
    {
        if(i < data_len - time_lag_ind)
            temp1[i] = data_lagged[i] * data[i];  
    }

    double f_mean = average(data, data_len);

    double variance_num = average(temp1, data_len - time_lag_ind) - (f_mean * f_mean);
    
    ret = variance_num/variance_denom;
        
    free(data_lagged);
    free(temp1);
    return ret;
}

/* ***************************************
*
* Block average
*
* Parameters
* ----------
*  data       - the raw data from which the
*               autocorrelation should be
                calculated
*  data_len   - the length of data
*  block_size - the size of the block
*
* Returns
* -------
* The statistical inefficency for a given
* block size
*
* ***************************************/
double block_average(double *data,
	                int data_len,
	                int block_size
	               )
{
    int m_b = data_len / block_size;

    double *block_data = (double *)malloc(m_b * sizeof(double));

    int i,j;
    double temp;
    for(j = 0; j < m_b; j++)
    {
        temp = 0;
        for(i = 0; i < block_size; i++)
        {
            temp += data[i + j * block_size];
        }
        block_data[j] = temp/block_size;
    }

    double ret = 0;

    ret = block_size * (variance(block_data, m_b)) / variance(data, data_len);

    free(block_data);

    return ret;
}