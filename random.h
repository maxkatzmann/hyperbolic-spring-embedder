#pragma once

#include <gsl/gsl_rng.h>

void initRNG(long long seed);

// Returns a random integer from [0, M-1].
int randint(int M);

// Returns a random double from [0,1)
double randdbl();

// Returns a random double from (0,1)
double randdblpos();

double randgaussian(double sigma);

gsl_rng* rng();
