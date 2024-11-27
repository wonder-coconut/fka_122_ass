#pragma once
#include <gsl/gsl_rng.h>
#include "result_t.h"

result_t MCMC_step_displace_all(double *,double ,gsl_rng *);