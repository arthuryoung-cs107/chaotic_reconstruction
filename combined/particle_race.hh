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
    const int nleaders;
    /** The length of the parameter vector corresponding to a particle */
    const int param_len;

    double **pool=NULL, **leaders = NULL, *l2_score=NULL, *l2_pscore=NULL;
    int *frame_score=NULL, *frame_pscore=NULL;

    referee(int npool_, int n_leaders_, int param_len_, double dt_sim_, double gau_scale_): npool(npool_), nleaders(n_leaders_), param_len(param_len_), dt_sim(dt_sim_), gau_scale(gau_scale_), gau_coeff(0.5/(gau_scale*gau_scale)) {}
    referee(referee &ref_): npool(ref_.npool), nleaders(ref_.nleaders), param_len(ref_.param_len), dt_sim(ref_.dt_sim), gau_scale(ref_.gau_scale), gau_coeff(0.5/(gau_scale*gau_scale)) {}
    ~referee()
    {
      if (pool!=NULL) free_AYdmatrix(pool);
      if (leaders!=NULL) free_AYdmatrix(leaders);
      if (l2_score!=NULL) delete l2_score;
      if (l2_pscore!=NULL) delete l2_pscore;
      if (frame_score!=NULL) delete frame_score;
      if (frame_pscore!=NULL) delete frame_pscore;
    }
    void alloc_records()
    {
      pool = AYdmatrix(npool, param_len);
      leaders = AYdmatrix(nleaders, param_len);
      l2_score = new double[nleaders];
      frame_score = new int[nleaders];
      l2_pscore = new double[npool];
      frame_pscore = new int[npool];
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

        int *sorting_pool=NULL;

        bool check_pool_results();
        void resample_pool();

#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

#endif
