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

struct system_struct
{
  system_struct(int ncomp_, int nstates_): ncomp(ncomp_), nstates(nstates_) {}
  ~system_struct() {}

  const int ncomp, // nbeads
            nstates; // Frames
};

struct event_detector: public virtual system_struct
{
  event_detector(int ncomp_, int nstates_, int dof_, double alpha_tol_): system_struct(ncomp_,nstates_), dof(dof_), ndof(ncomp*dof), alpha_tol(alpha_tol_),
  r2_state_comp(Tmatrix<double>(nstates,ncomp)), alpha_state_comp(Tmatrix<double>(nstates,ncomp)),
  INTr2_comp_history(Tmatrix<double>(ncomp,3)), r2_regime_comp(Tmatrix<double>(ncomp,ncomp)) {}
  ~event_detector()
  {
    free_Tmatrix<double>(r2_state_comp); free_Tmatrix<double>(alpha_state_comp);
    free_Tmatrix<double>(INTr2_comp_history); free_Tmatrix<double>(r2_regime_comp);
  }

  const int dof, // x,y
            ndof; // dof per state
        int stev_early,
            stev_late;

  const double alpha_tol;

  double  ** const r2_state_comp,
          ** const alpha_state_comp,
          ** const INTr2_comp_history,
          ** const r2_regime_comp;

  inline double alpha_comp(double *a_, double t_p1, double t_m1)
  {return log(a_[0]/a_[2])/log(t_p1/t_m1);}
  inline double alpha_comp(double a_p1, double a_m1, double t_p1, double t_m1)
  {return log(a_p1/a_m1)/log(t_p1/t_m1);}
};

struct event_block: public virtual system_struct
{
  event_block(int ncomp_, int nstates_): system_struct(ncomp_, nstates_),
  stev_comp(new int[ncomp]), stev_ordered(new int[ncomp]), comps_ordered(new int[ncomp]),
  nev_state_comp(Tmatrix<int>(nstates,ncomp)), nobs_state_comp(Tmatrix<int>(nstates,ncomp)),
  rho2_regime(new double[ncomp]),
  mur2_state_comp(Tmatrix<double>(nstates,ncomp)), stdr2_state_comp(Tmatrix<double>(nstates,ncomp)),
  mualpha_state_comp(Tmatrix<double>(nstates,ncomp)), stdalpha_state_comp(Tmatrix<double>(nstates,ncomp)) {}
  ~event_block()
  {
    delete [] stev_comp; delete [] stev_ordered; delete [] comps_ordered;
    free_Tmatrix<int>(nev_state_comp); free_Tmatrix<int>(nobs_state_comp);
    delete [] rho2_regime;
    free_Tmatrix<double>(mur2_state_comp); free_Tmatrix<double>(stdr2_state_comp);
    free_Tmatrix<double>(mualpha_state_comp); free_Tmatrix<double>(stdalpha_state_comp);
  }

  int stev_earliest,
      stev_latest;

  int * const stev_comp,
      * const stev_ordered,
      * const comps_ordered,
      ** const nev_state_comp,
      ** const nobs_state_comp;

  double  * const rho2_regime,
          ** const mur2_state_comp,
          ** const stdr2_state_comp,
          ** const mualpha_state_comp,
          ** const stdalpha_state_comp;


  virtual void synchronise_event_data(int stev_earliest_, int stev_latest_, int *stev_c_, int *stev_o_, int *comps_o_, double *rho2_r_);
  virtual void consolidate_event_data(int ntest_);
  virtual void clear_event_data();
  virtual void define_event_block(double sigma_scaled_);

  inline void report_event_data(int ** nev_s_c_, int ** nobs_s_c_, double ** r2_s_c_, double **alpha_s_c_)
  {
    for (int i = 0; i < ncomp*stev_latest; i++)
    {
      nev_s_c_[0][i] += nev_state_comp[0][i];
      nobs_s_c_[0][i] += nobs_state_comp[0][i];
      r2_s_c_[0][i] += mur2_state_comp[0][i];
      alpha_s_c_[0][i] += mualpha_state_comp[0][i];
    }
  }
  inline void set_state_counts(int &nstate_obs_, int &nstate_stable_, int &nstate_unstable_)
  {
    int nstate_stable_local=0;
    for (int i = 0; i < ncomp; i++) nstate_stable_local+=stev_comp[i];
    nstate_obs_=ncomp*stev_ordered[ncomp-1];
    nstate_stable_=nstate_stable_local;
    nstate_unstable_=nstate_obs_-nstate_stable_;
  }
  inline int earliest(int *frames_ordered_, int &early_index_, int index_start_)
  {
    int out = frames_ordered_[index_start_]; early_index_= index_start_;
    for (int i = index_start_+1; i < nbeads; i++)
      if (frames_ordered_[i])<out) out=frames_ordered_[early_index_=i];
    return out;
  }
  inline void earliest_recursive(int *frames_ordered_, int *beads_ordered_ int i_next)
  {
    int early_index, early_frame, index_temp, frame_temp;
    early_frame = earliest(frames_ordered_,early_index,i_next);
    index_temp=i_next; frame_temp=frames_ordered_[index_temp];
    beads_ordered_[index_temp]=early_index; frames_ordered_[index_temp]=early_frame;
    beads_ordered_[early_index]=index_temp; frames_ordered_[early_index]=frame_temp;
    if (i_next<nbeads-1) earliest_recursive(frames_ordered_,beads_ordered_,i_next+1)
  }
};

struct gaussian_likelihood
{
  gaussian_likelihood(double sigma_noise_, double cl_): sigma_noise(sigma_noise_), cl(cl_), sigma_scaled(cl*sigma_noise) {}
  ~gaussian_likelihood() {}

  const double  sigma_noise,
                cl,
                sigma_scaled;

  inline double compute_weight(double r_, double r_min_, double rho_)
  {
    // arguments provided unsquared so that we can avoid issues with machine precision
    return exp(0.5*((r_min_-r_)/rho_)*((r_min_+r_)/rho_));
  }
  inline double compute_prob(double r_, double rho_)
  {return (1.0/(rho_*sqrt(root2pi)))*exp(-0.5*(r_/rho_)*(r_/rho_));}  
  inline double expected_r2(int * obs_states_, int ncomp_, int ndof_=2)
  {
    int net_obs=0;
    for (int i = 0; i < ncomp_; i++) net_obs+=ndof_*obs_states_[i];
    return sigma_scaled*sigma_scaled*((double)net_obs);
  }

  private:
    const double root2pi=sqrt(2.0*M_PI);

};

double

#endif
