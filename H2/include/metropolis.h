#pragma once
#include <gsl/gsl_rng.h>
#include "mcmc_step.h"

mcmc_step mcmc_displace_all(double *r1, double *r2, double d, double alpha, gsl_rng *r);