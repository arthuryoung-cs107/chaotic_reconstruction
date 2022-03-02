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
  int index;
  int rank;

  double * params;

  record();
  ~record();

  int check_success(int frscore_, double l2score_, int frmin_, double l2min_)
  {
    frscore = frscore_; l2score = l2score_;
    if (frscore>frmin_) return 0; // worse then worst leader run, by frame score
    else if (frscore<frmin_) return 1; // better than best leader run, by frame score
    else if (l2score<l2min_) return 1; // better than best leader run, by l2 error
    else return 0; // worse then worst leader run by l2 error
  }
};

struct referee
{
    /** The maximum simulation timestep to use. */
    const double dt_sim;
    /** The length scale of the Gaussian (in pixels) for updating the particle
     * weights. */
    const double gau_scale;
    /** The coefficient that appears in the Gaussian for updating the particle
     * weights. */
    const double gau_coeff;

    /** number of particles being tested in each generation */
    const int npool;
    /** The number of best performing particles we store */
    const int nlead;
    /** The length of the parameter vector corresponding to a particle */
    const int param_len;

    // structure of records with which we will be using to compare particles
    record *leaders, *pool;

    // memory chunks to store parameters associated with particles
    double **pool_params, **lead_params;

    // integer vectors indicating the relative ranking of particles
    int *lead_ranking, *pool_ranking;

    bool alloc_flag;

    referee(int n_leaders_, int npool_, int param_len_, double dt_sim_, double gau_scale_): nlead(n_leaders_), npool(npool_), param_len(param_len_), dt_sim(dt_sim_), gau_scale(gau_scale_), gau_coeff(0.5/(gau_scale*gau_scale)) {}

    referee(referee &ref_): nlead(ref_.nlead), npool(ref_.npool), param_len(ref_.param_len), dt_sim(ref_.dt_sim), gau_scale(ref_.gau_scale), gau_coeff(0.5/(gau_scale*gau_scale)) {}
    ~referee()
    {
      if (alloc_flag)
      {
        delete lead_ranking;
        delete pool_ranking;
        delete leaders;
        delete pool_leaders;
        free_AYdmatrix(lead_params);
        free_AYdmatrix(pool_params);
      }
    }
    void alloc_records()
    {
      lead_ranking = new int[nlead];
      pool_ranking = new int[nlead];

      leaders = new record[nlead];
      pool_leaders = new record[nlead];

      lead_params = AYdmatrix(nlead, param_len);
      pool_params = AYdmatrix(npool, param_len);

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

#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

#endif
