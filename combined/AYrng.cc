#include "AYlinalg.hh"
#include <math.h>

extern "C"
{
  #include "AYaux.h"
}

AYrng::AYrng(): seed((uint64_t)1e15), carry(lcg_fwd(seed, jump))
{}
AYrng::AYrng(uint64_t seed_, uint64_t jump_): seed(seed_), jump(seed_), carry(lcg_fwd(seed, jump)) {}
AYrng::~AYrng() {}
void AYrng::rng_init(uint64_t seed_, uint64_t jump_)
{seed = seed_; jump = jump_; carry = lcg_fwd(seed, jump);}
double AYrng::rand_gen() { return 0.0;}

uniform::uniform(double low_, double high_): AYrng(), low(low_), high(high_)
{}
uniform::~uniform(){}
double uniform::rand_gen() {return knuth_random_uni(low, high, &carry);}

normal::normal(double mu_, double var_): AYrng(), mu(mu_), var(var_) {}
normal::~normal() {}
double normal::rand_gen()
{return boxmuller_knuth(mu, var, &carry);}
