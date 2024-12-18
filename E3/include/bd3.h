#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "result_t.h"

#pragma once

result_t BD3(double initial_position, double initial_velocity, double w0, double dt, double eta, double kB, double mass, double T, gsl_rng *k);