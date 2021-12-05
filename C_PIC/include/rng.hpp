#ifndef RANDNUMS
#define RANDNUMS 
//Taken from Numerical Recipes in C (3rd Edition)
typedef unsigned long long rng_t;
typedef struct
{
    rng_t u, v = 4101842887655102017LL, w = 1;
    double mu, sigma;
} RandGen;
RandGen* initializeRNG(rng_t j);
RandGen* initializeRNG(rng_t j, double mu, double sigma);
inline rng_t int64RNG(RandGen* rng)
{
    rng->u = rng->u * 2862933555777941757LL + 7046029254386353087LL;
    rng->v ^= rng->v >> 17; rng->v ^= rng->v << 31; rng->v = rng->v >> 8;
    rng->w = 4294957665U*(rng->w & 0xffffffff) + (rng->w >> 32);
    rng_t x = rng->u ^ (rng->u << 21); x ^= x >> 35; x ^= x << 4;
    return (x+rng->v) ^ rng->w;
};
inline double doubleRNG(RandGen* rng)
{
    return 5.42101086242752217e-20 * int64RNG(rng);
};
inline unsigned int int32RNG(RandGen* rng)
{
    return (unsigned int)int64RNG(rng);
};

//Based on Leva normal deviate quadratic curve method from Numerical Recipes in C 3rd edition
inline double normal_deviation(RandGen* rng)
{
    double u, v, x, y, q;
    do
    {
        u = doubleRNG(rng);
        v = 1.7156*(doubleRNG(rng)-0.5);
        x = u - 0.449871;
        y = fabs(v) + 0.386595;
        q = x*x + y*(0.19600*y-0.25472*x);
    }while(q > 0.27597 && (q > 0.27846 || v*v > -4*log(u)*u*u));
    return rng->mu + rng->sigma*v/u;
}
#endif
