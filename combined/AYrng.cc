#include "AYlinalg.hh"
#include <math.h>

extern "C"
{
  #include "AYaux.h"
}

AYrng::AYrng(): jump(1000), seed((uint64_t) rand()) {}
AYrng::~AYrng() {}
void AYrng::rng_init(uint64_t seed_, uint64_t jump_) {}
double AYrng::rand_gen() { return 0.0;}

uniform::uniform(double low_, double high_): AYrng(), low(low_), high(high_) {}
uniform::~uniform(){}
void uniform::rng_init(uint64_t seed_, uint64_t jump_)
{
  if (seed_ != 0){seed = seed_;}
  if (jump_ != 0){jump = jump_;}
  carry = lcg_fwd(seed, jump);
}
double uniform::rand_gen() {return knuth_random_uni(low, high, &carry);}

normal::normal(double mu_, double var_): AYrng(), mu(mu_), var(var_) {}
normal::~normal() {}
void normal::rng_init(uint64_t seed_, uint64_t jump_)
{
  if (seed_ != 0){seed = seed_;}
  if (jump_ != 0){jump = jump_;}
  carry = lcg_fwd(seed, jump);
}
double normal::rand_gen()
{
  return boxmuller_knuth(mu, var, &carry);
}
