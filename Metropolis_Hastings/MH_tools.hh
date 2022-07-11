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
void print_row_vec(int *ivec_, int n_);
void print_row_vec(double *dvec_, int n_);

struct MH_rng
{
  MH_rng(int seed_=1): seed(seed_), gsl_gen(gsl_rng_alloc(gsl_rng_taus2)) {}
  ~MH_rng() {gsl_rng_free(gsl_gen);}

  const int seed;
  gsl_rng * const gsl_gen;

  inline double rand_uni(double low_=0.0,double high_=1.0) {return (gsl_rng_uniform(gsl_gen)-0.5)*(high_-low_) + 0.5*(high_+low_);}

  inline double rand_gau(double mu_=0.0, double sigma_=1.0) {return mu_+gsl_ran_gaussian(gsl_gen,sigma_);}
};

struct event_block
{
  event_block(int ncomp_, int nstates_, int dim_=2);
  ~event_block();

  bool stable_flag;

  const int ncomp, // nbeads, in a general sense
            nstates, // Frames, in a general sense
            dim, // dimension of each component (x,y)
            dof; // dof per state (i.e. sum x and y for every bead)

  int stev_earliest, // state index of earliest detected event (pertains to one system component)
      stev_latest, // state index of latest detected event (pertains to one system component)
      * const stev_comp, // state index of events for each component in the system
      * const stev_ordered, // state indices of events arranged chronologically
      * const comps_ordered, // system component indices arranged in order of their experienced events
      ** const nev_state_comp, // counts of events experienced at each sequence state for many trials
      ** const nobs_state_comp; // number of observations of system state for many trials of event detection

  double  rho2_objective,
          * const rho2stable_comp, // expected STABLE residual for each system component in the event block
          * const rho2unstable_comp,
          ** const mur2_state_comp, // workspaces for collecting statistics on state sequence
          ** const stdr2_state_comp, // ...
          ** const mualpha_state_comp, // ...
          ** const stdalpha_state_comp; // ...

  // debugging

  inline void print_event_block(double sigma_scaled_, const char indent_[]="     ")
  {
    printf("%s(event_block) stev_earliest: %d, stev_latest: %d, nst_full: %d, nst_stable: %d, nst_unstable: %d\n",
    indent_,
    stev_earliest,
    stev_latest,
    nst_full(),
    nst_stable(),
    nst_unstable());
    printf("%s(event_block) rho2_objective: %e, netrho2: %e, rho2stable: %e, rho2unstable: %e\n",
    indent_,
    rho2_objective,
    compute_netrho2(),
    compute_netrho2stable(),
    compute_netrho2unstable());
    printf("%s(event_block) stev_comp: ",indent_); print_row_vec(stev_comp, ncomp);
    printf("%s(event_block) stev_ordered: ",indent_); print_row_vec(stev_ordered, ncomp);
    printf("%s(event_block) rho2stable_comp: ",indent_); print_row_vec(rho2stable_comp, ncomp);
    printf("%s(event_block) rho2unstable_comp: ",indent_); print_row_vec(rho2unstable_comp, ncomp);
  }

  virtual void clear_event_data();
  virtual void consolidate_event_data();
  virtual void define_event_block(double sigma_scaled_);
  virtual void synchronise_event_data(int stev_e_, int stev_l_, int *stev_c_, int *stev_o_, int *comps_o_,double *rho2s_c_, double *rho2us_c_);
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
  inline int earliest(int *frames_ordered_, int &early_index_, int index_start_)
  {
    int out = frames_ordered_[index_start_]; early_index_= index_start_;
    for (int i = index_start_+1; i < ncomp; i++)
      if (frames_ordered_[i]<out) out=frames_ordered_[early_index_=i];
    return out;
  }
  inline void earliest_recursive(int *stev_ordered_, int *comps_ordered_, int i_next)
  {
    int early_index, // index position of next earliest event frame in stev_ordered_
        early_frame=earliest(stev_ordered_,early_index,i_next), // value of next earliest event frame
        early_bead=comps_ordered_[early_index]; // bead id associated with next earliest event frame

    stev_ordered_[early_index]=stev_ordered_[i_next];
    stev_ordered_[i_next]=early_frame;

    comps_ordered_[early_index]=comps_ordered_[i_next];
    comps_ordered_[i_next]=early_bead;
    if (i_next<ncomp-1) earliest_recursive(stev_ordered_,comps_ordered_,i_next+1);
  }

  inline bool check_stev_convergence()
  {
    bool conv_out=true;
    for (int i = 0; i < ncomp; i++) if (stev_comp[i]<(nstates-1)) conv_out=false;
    return conv_out;
  }

  inline int stev_last() {return stev_ordered[ncomp-1];}
  inline int nst_full() {return stev_ordered[ncomp-1]*ncomp;}
  inline int nst_stable() {int out=0; for (int i = 0; i < ncomp; i++) out+=stev_comp[i]; return out;}
  inline int nst_unstable() {return nst_full()-nst_stable();}

  inline double compute_netrho2()
    {double out=0.0; for (int i = 0; i < ncomp; i++) out+=rho2stable_comp[i]+rho2unstable_comp[i]; return out;}
  inline double compute_netrho2stable()
    {double out=0.0; for (int i = 0; i < ncomp; i++) out+=rho2stable_comp[i]; return out;}
  inline double compute_netrho2unstable()
    {double out=0.0; for (int i = 0; i < ncomp; i++) out+=rho2unstable_comp[i]; return out;}
};

struct event_detector: public event_block
{
  event_detector(int ncomp_, int nstates_, int dim_, double alpha_tol_): event_block(ncomp_,nstates_,dim_),
  alpha_tol(alpha_tol_),
  r2_state_comp(Tmatrix<double>(nstates,ncomp)), alpha_state_comp(Tmatrix<double>(nstates,ncomp)),
  INTr2_comp_history(Tmatrix<double>(ncomp,3)) {}

  ~event_detector()
  {
    free_Tmatrix<double>(r2_state_comp); free_Tmatrix<double>(alpha_state_comp);
    free_Tmatrix<double>(INTr2_comp_history);
  }

  int stev_early,
      stev_late;

  const double alpha_tol;

  double  ** const r2_state_comp,
          ** const alpha_state_comp,
          ** const INTr2_comp_history;

  inline void update_integral_history(double INT_now_, int ibead_)
  {
    INTr2_comp_history[ibead_][2]=INTr2_comp_history[ibead_][1];
    INTr2_comp_history[ibead_][1]=INTr2_comp_history[ibead_][0];
    INTr2_comp_history[ibead_][0]+=INT_now_;
  }

  inline double alpha_comp(double *a_, double t_p1, double t_m1)
    {return log(a_[0]/a_[2])/log(t_p1/t_m1);}
  inline double alpha_comp(double a_p1, double a_m1, double t_p1, double t_m1)
    {return log(a_p1/a_m1)/log(t_p1/t_m1);}
};

struct gaussian_likelihood
{
  gaussian_likelihood(double sigma_noise_, double cl_): sigma_noise(sigma_noise_), noise_scale(cl_), sigma_scaled(noise_scale*sigma_noise) {}
  ~gaussian_likelihood() {}

  const double  sigma_noise,
                noise_scale,
                sigma_scaled;

  inline void print_gaussian_likelihood(const char indent_[])
    {printf("%s(gaussian_likelihood) sigma_noise = %e, noise_scale = %e, sigma_scaled = %e\n", indent_,sigma_noise,noise_scale,sigma_scaled);}

  inline double compute_weight(double r_, double r_min_, double rho_)
    {// arguments provided unsquared so that we can avoid issues with machine precision
    return exp(0.5*((r_min_-r_)/rho_)*((r_min_+r_)/rho_));}

  inline double compute_prob(double r_, double rho_)
    {return (1.0/(rho_*root2pi))*exp(-0.5*(r_/rho_)*(r_/rho_));}

  inline double comp_rho2(int * obs_states_, int ncomp_, int dim_=2)
  {
    int net_stobs=0; // number of observed system STATES
    for (int i = 0; i < ncomp_; i++) net_stobs+=obs_states_[i];
    return (sigma_scaled*((double)(dim_*net_stobs)))*sigma_scaled;
  }

  private:
    const double root2pi=sqrt(2.0*M_PI);
};

#endif
