#pragma once

#ifndef MCMC_STEP_H
#define MCMC_STEP_H

typedef struct{
    double energy;
    double ln_psi_grad;
    int accepted;
} mcmc_step;

#endif