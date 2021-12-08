#include "rng.hpp"
#include <cmath>
#include <stdlib.h>
RandGen* initializeRNG(rng_t j)
{
    RandGen* rng = new RandGen;
    rng->u = j ^ rng->v;
    rng->v = rng->u;
    rng->w = rng->v;
    return rng;
};
RandGen* initializeRNG(rng_t j, double mu, double sigma)
{
    RandGen* rng = new RandGen;
    rng->u = j ^ rng->v;
    rng->v = rng->u;
    rng->w = rng->v;
    rng->mu = mu;
    rng->sigma = sigma;
    return rng;
};