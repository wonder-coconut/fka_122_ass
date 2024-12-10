#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tools.h"

double variance(double *v, unsigned int len)
{
    //square of standard deviation
    double sd = standard_deviation(v,len);
    return sd*sd;
}


double autocorrelation(double *data, int data_len, int time_lag_ind)
{
    //autocorrelation value var
    double ret = 0;
    
    //denominator calc
    double denom = variance(data, data_len);

    double *data_lagged = (double *)malloc((data_len - time_lag_ind) * sizeof(double));
    double *temp1 = (double *)malloc((data_len - time_lag_ind) * sizeof(double));
    
    //shifting data by time lag
    int i;
    for(i = 0; i < data_len - time_lag_ind; i++)
        data_lagged[i] = data[i + time_lag_ind];
    
    //data lagged * original data
    for(i = 0; i < data_len; i++)
    {
        if(i < data_len - time_lag_ind)
            temp1[i] = data_lagged[i] * data[i];
    }
    //mean
    double data_mean = average(data,data_len);
    //numerator
    double num = average(temp1, data_len - time_lag_ind) - (data_mean * data_mean);
    
    //autocorrelation
    ret = num/denom;
    
    free(data_lagged);
    free(temp1);
    
    return ret;
}

double block_average(double *data, int data_len, int block_size)
{
    //division of blocks
    int m_b = data_len / block_size;
    double *block_data = (double *)malloc(m_b * sizeof(double));

    int i,j;
    double temp;
    //iteration through all blocks
    for(j = 0; j < m_b; j++)
    {
        temp = 0;
        //average for each block
        for(i = 0; i < block_size; i++)
            temp += data[i + j * block_size];
        
        block_data[j] = temp/block_size;
    }

    double ret = 0;
    //statistical inefficiency
    ret = block_size * (variance(block_data, m_b)) / variance(data, data_len);

    free(block_data);
    return ret;
}