#include "random.h"

#include <gsl/gsl_randist.h>

// TODO make singleton
gsl_rng* _rng;

void initRNG(long long seed) {
  const gsl_rng_type* T;
  gsl_rng_env_setup();

  T = gsl_rng_default;
  _rng = gsl_rng_alloc(T);
  gsl_rng_set(_rng, seed);
}

int randint(int M) {
  return (int)gsl_rng_uniform_int(_rng, (unsigned long int)M);
}

double randdbl() { return gsl_rng_uniform(_rng); }

double randdblpos() { return gsl_rng_uniform_pos(_rng); }

double randgaussian(double sigma) { return gsl_ran_gaussian(_rng, sigma); }

gsl_rng* rng() { return _rng; }
