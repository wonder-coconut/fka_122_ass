#pragma once

#ifndef MCMC_STEP_H
#define MCMC_STEP_H

typedef struct{
    double energy;
    double probability;
    int accepted;
} mcmc_step;

#endif