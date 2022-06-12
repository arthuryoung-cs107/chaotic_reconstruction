#ifndef MH_TOOLS_HH
#define MH_TOOLS_HH

#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

template <typename T> T ** Tmatrix(int M_, int N_)
{
  T * chunk = new T[M_*N_],
    ** rows = new T*[M_];
  for (int i = 0,j=0; i < M_; i++,j+=N_)
    rows[i] = chunk+j;
  return rows;
}

template <typename T> void free_Tmatrix(T ** Tmat_)
{
  delete [] Tmat_[0];
  delete [] Tmat_;
}


void fseek_SAFE(FILE *fp,long int offset,int origin);
void fread_SAFE(void *ptr,size_t size,size_t count,FILE *fp);

struct MH_rng
{
  MH_rng(int seed_=1): seed(seed_), gsl_gen(gsl_rng_alloc(gsl_rng_taus2)) {}
  ~MH_rng() {gsl_rng_free(gsl_gen);}

  const int seed;
  gsl_rng * const gsl_gen;

  inline double rand_uni(double low_=0.0,double high_=1.0) {return (gsl_rng_uniform(gsl_gen)-0.5)*(high_-low_) + 0.5*(high_+low_);}

  inline double rand_gau(double mu_=0.0, double sigma_=1.0) {return mu_+gsl_ran_gaussian(gsl_gen,sigma_);}
};

struct gaussian_likelihood
{
  gaussian_likelihood(double noise_tol_, double cl_): noise_tol(noise_tol_), cl(cl_), coeff(cl*noise_tol) {}
  ~gaussian_likelihood() {}

  const double  noise_tol,
                cl,
                coeff;

  inline double expected_r2(int * obs_, int nbeads_)
  {
    int net_obs=0;
    for (int i = 0; i < nbeads_; i++) net_obs+=obs_[i];
    return coeff*coeff*((double)net_obs);
  }
};

struct event_detector
{
  event_detector(int ncomp_, int nstates_, double alpha_tol_):  ncomp(ncomp_), nstates(nstates_), alpha_tol(alpha_tol_),
  evcount_comp_state(Tmatrix<int>(ncomp,nstates)),
  r2_state_comp(Tmatrix<double>(nstates,ncomp)), alpha_state_comp(Tmatrix<double>(nstates,ncomp)),r2_comp_regime(Tmatrix<double>(ncomp,ncomp)), INTr2_comp_history(Tmatrix<double>(ncomp,3)) {}
  ~event_detector() {free_Tmatrix<int>(evcount_comp_state);
  free_Tmatrix<double>(r2_state_comp);
  free_Tmatrix<double>(alpha_state_comp);
  free_Tmatrix<double>(r2_comp_regime);
  free_Tmatrix<double>(INTr2_comp_history);}


  const int ncomp, // nbeads
            nstates; // Frames

  int ** const evcount_comp_state;

  const double alpha_tol;

  double  ** const r2_state_comp,
          ** const alpha_state_comp,
          ** const r2_comp_regime,
          ** const INTr2_comp_history;

  inline double alpha_comp(double *a_, double t_p1, double t_m1)
  {return log(a_[0]/a_[2])/log(t_p1/t_m1);}
  inline double alpha_comp(double a_p1, double a_m1, double t_p1, double t_m1)
  {return log(a_p1/a_m1)/log(t_p1/t_m1);}
};

#endif
