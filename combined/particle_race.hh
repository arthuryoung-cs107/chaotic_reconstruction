#ifndef RACE_HH
#define RACE_HH

#include <cstdio>
#include <gsl/gsl_rng.h>
#include "swirl.hh"
#include "RYdat2AYdat.hh"
#include "AYlinalg.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

extern "C"
{
  #include "AYaux.h"
}

struct record
{
  int frscore;
  int l2score;
  int global_index;
  bool success;

  double * params;

  record() {}
  record(int i_): global_index(i_) {}
  ~record()

  inline void check_success(int frscore_, double l2score_, int frmin_, double l2min_)
  {
    frscore = frscore_; l2score = l2score_;
    if (frscore>frmin_) success = false; // worse than worst leader run, by frame score
    else if (frscore<frmin_) success = true; // better than best leader run, by frame score
    else if (l2score<l2min_) success = true; // better than best leader run, by l2 error
    else success = false; // worse than worst leader run by l2 error
  }
  inline bool isbetter(record * rcomp_)
  {
    if (frscore>rcomp_->frscore) return false;
    else if (frscore<rcomp_->frscore) return true;
    else if (l2score<rcomp_->l2score) return true;
    else return false;
  }

  inline void take_vals(record * rtake_, int len_ )
  {
    frscore = rtake_->frscore; l2score = rtake_->l2score;
    for (int i = 0; i < len_; i++) params[i] = rtake_->params[i];
  }

  inline bool isworse(record * rcomp_)
  {return !(isbetter(rcomp_));}
};

struct referee
{
    /** The maximum simulation timestep to use. */
    const double dt_sim;
    /** the perturbation variance for resampled particles */
    const double gau_var;
    // the strength of the exponential weight function for the resampling of particles
    const double lambda;


    /** number of particles being tested in each generation */
    const int npool;
    /** The number of best performing particles we store */
    const int nlead;
    /** The length of the parameter vector corresponding to a particle */
    const int param_len;

    // structure of records with which we will be using to compare particles
    record **leaders, **pool, **leader_board, **pool_leaders;

    // memory chunks to store parameters associated with particles
    double **pool_params, **lead_params, *sample_weights;

    bool alloc_flag;

    referee(int n_leaders_, int npool_, int param_len_, double dt_sim_, double gau_var_, double lambda_): nlead(n_leaders_), npool(npool_), param_len(param_len_), dt_sim(dt_sim_), gau_var(gau_var_), lambda(lambda_) {}

    referee(referee &ref_): nlead(ref_.nlead), npool(ref_.npool), param_len(ref_.param_len), dt_sim(ref_.dt_sim), gau_var(ref_.gau_var), lambda(ref_.lambda) {}
    ~referee()
    {
      if (alloc_flag)
      {
        delete leader_board;
        for (int i = 0; i < 2*nlead; i++) delete leaders[i];
        delete leaders;
        free_AYdmatrix(lead_params);
      }
    }
    void alloc_records()
    {
      leaders = new record*[nlead+npool];
      leader_board = new record*[nlead+npool];
      for (int i = 0; i < nlead+npool; i++) leaders[i] = new record(i);
      pool = leaders + nlead;
      pool_leaders = leader_board + nlead;

      params = AYdmatrix(nlead+npool, param_len);
      pool_params = params + nlead;

      alloc_flag = true;
    }
};

class runner : public swirl
{
    public:
      const int thread_id, param_len, Frames;
      double *pvals, *x0, t0, ctheta0, pos_err_acc, tol;
      int frame;

      runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, double tol_);
      ~runner();

    private:
      void init_ics(double *x0_, double t0_, double ctheta0_);
      void reset_sim(double *ptest_);
      void run_race(int frames, double *ts_, double *xs_, double *d_ang_);
      bool is_lost(double *f_);
};

class race : public referee {
    public:
        /** The minimum parameters. */
        swirl_param sp_min;
        /** The maximum parameters. */
        swirl_param sp_max;
        /** The simulation time in seconds. */
        double t_phys;
        /** The total number of beads. */
        int n;
        /** The total number of snapshots. */
        int nsnap;
        /** The total number of particles. */
        int npar;
        /** The time points of the snapshots. */
        double* ts;
        /** The bead positions at the snapshots. */
        double* xs;
        /** The dish angle data at the snapshots. */
        double* d_ang;

        /** A count of the recorded best performing particles */
        int leader_count;
        /** A count of the generations of tested particles */
        int gen_count;
        /** frame score associated w/ poorest particle on leaderboard  */
        int frscore_min;
        /** error score associated w/ poorest particle on leaderboard  */
        int l2score_min;

        /** the data index we take to be the initial conditions of the experiment */
        const int ic_index;

        ODR_struct odr;

        race(referee &rparam,swirl_param &sp_min_,swirl_param &sp_max_,swirl_param &sp_rnd_,wall_list &wl_,double t_phys_,ODR_struct &odr_,int ic_index_=0);
        ~race();
        void init_race(int npar_);
    private:
        /** The number of threads. */
        const int nt;
        /** The array AYuniform generators. */
        AYuniform * uni;
        /** A reference to the list of walls for the swirling simulation. */
        wall_list &wl;
        /** The array of GSL random number generators. */
        gsl_rng** const rng;
        /** The array of proximity grids. */
        proximity_grid** const pg;

        /** Array of swirl simulations used by each thread to test particles */
        runner ** runners;
        /** count of particles in current pool beating worst leader particles */
        int pool_success_count;

        bool check_pool_results();
        void resample_pool();
        int collect_pool_leaders();

#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

int find_worst(record ** r, int ncap)
{
  int worst_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[worst_index]->isbetter(r[i]))
      worst_index = i;
  return worst_index;
}

#endif
