#include "AYlinalg.hh"
#include <math.h>

extern "C"
{
  #include "AYaux.h"
}
AYrng::AYrng(): carry(lcg_fwd(seed, jump)) {}
AYrng::AYrng(uint64_t seed_) {seed+=seed_; carry=lcg_fwd(seed, jump);}
AYrng::AYrng(uint64_t seed_, uint64_t jump_): seed(seed_), jump(seed_), carry(lcg_fwd(seed, jump)) {}
AYrng::~AYrng() {if (gsl_gen!=NULL) gsl_rng_free(gsl_gen);}
void AYrng::rng_true_init(uint64_t seed_, uint64_t jump_)
{seed = seed_; jump = jump_; carry = lcg_fwd(seed, jump);}
void AYrng::rng_init(uint64_t seed_)
{seed+=seed_; jump=100; carry=lcg_fwd(seed, jump);}
double AYrng::dseed() {return ((double)seed)/(lcg_sze());}
double AYrng::dcarry() {return ((double)carry)/(lcg_sze());}
void AYrng::rng_init_gsl(uint64_t seed_)
{
  gsl_gen = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(gsl_gen, seed_);
}
double AYrng::rand_uni_gsl(double low_, double high_)
{
  if (gsl_gen==NULL) rng_init_gsl();
  return (gsl_rng_uniform(gsl_gen)-0.5)*(high_-low_) + ((high_+low_)/2.0);
}
double AYrng::rand_gau_gsl_old(double mu_, double var_)
{
  if (gsl_gen==NULL) rng_init_gsl();
  double uni1=gsl_rng_uniform(gsl_gen), uni2=gsl_rng_uniform(gsl_gen);
  return boxmuller_2uni(mu_, var_, uni1, uni2);
}


AYuniform::AYuniform(double low_, double high_): AYrng(), low(low_), high(high_) {}
AYuniform::AYuniform(double low_, double high_, uint64_t seed_): AYrng(seed_), low(low_), high(high_) {}
AYuniform::~AYuniform(){}
double AYuniform::rand_gen() {return knuth_random_uni(low, high, &carry);}

AYnormal::AYnormal(double mu_, double var_): AYrng(), mu(mu_), var(var_) {}
AYnormal::AYnormal(double mu_, double var_, uint64_t seed_): AYrng(seed_), mu(mu_), var(var_) {}
AYnormal::~AYnormal() {}
double AYnormal::rand_gen()
{return boxmuller_knuth(mu, var, &carry);}
