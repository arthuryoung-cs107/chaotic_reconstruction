#ifndef PARTICLE_RACE_HH
#define PARTICLE_RACE_HH

#include <cstdio>
#include <gsl/gsl_rng.h>
#include "swirl.hh"
#include "RYdat2AYdat.hh"
#include "AYlinalg.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

struct record
{
  int frscore=0;
  double l2score=0.0;
  int global_index;
  bool success;

  double * params;

  record() {}
  record(int i_): global_index(i_) {}
  ~record() {}

  inline bool check_success(int frscore_, double l2score_, int frmin_, double l2min_)
  {
    frscore = frscore_; l2score = l2score_;
    if (frscore<frmin_) success = false; // worse than worst leader run, by frame score
    else if (frscore>frmin_) success = true; // better than best leader run, by frame score
    else if (l2score<l2min_) success = true; // better than best leader run, by l2 error
    else success = false; // worse than worst leader run by l2 error
    return success;
  }
  inline bool isbetter(record * rcomp_)
  {
    if (frscore<rcomp_->frscore) return false;
    else if (frscore>rcomp_->frscore) return true;
    else if (l2score<rcomp_->l2score) return true;
    else return false;
  }

  void take_vals(record * rtake_, int len_ )
  {
    frscore = rtake_->frscore; l2score = rtake_->l2score;
    for (int i = 0; i < len_; i++) params[i] = rtake_->params[i];
  }

  inline bool isworse(record * rcomp_)
  {return !(isbetter(rcomp_));}

  inline double phi(int F_)
  {return ((double) frscore)/((double) F_);}

  inline double w(int F_, double lambda_)
  {return exp(lambda_*((double) frscore)/((double) F_));}

  inline double var(int F_, double var_)
  {return var_*exp(-2.0*sqrt(var_)*((double) frscore)/((double)F_));}

  void print_record()
  {
    printf("record ID:%d, fr:%d, l2:%e, params: ", global_index, frscore, l2score);
    for (int i = 0; i < 12; i++) printf("%e ", params[i]);
  }
};

class runner : public swirl
{
    public:
      const int thread_id, param_len, Frames;
      double *pvals, *x0, t0, ctheta0, pos_err_acc, tol, t_phys, t0_raw;
      int frame;

      runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, double tol_);
      ~runner();

      void init_ics(double t_phys_ , double *x0_, double t0_raw_, double ctheta0_);
      void reset_sim(double *ptest_);
      void run_race(double dt_sim_, double *ts_, double *xs_, double *d_ang_);

      void print_raw_ics();
      void print_params();
      void print_current_pos();
    private:
      bool is_lost(double *f_);
};

struct referee
{
    /** The maximum simulation timestep to use. */
    const double dt_sim;
    /** the perturbation variance for resampled particles */
    const double gau_var;
    /* the strength of the exponential weight function for the resampling of particles */
    const double lambda;
    /* the clearance we leave for uniform resampling when we HAVE NOT filled our leaderboard */
    const double rs_fill_factor; // maybe 0.5
    /* the clearance we leave for uniform resampling when we HAVE filled our leaderboard */
    const double rs_full_factor; // maybe 0.99?

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

    bool alloc_flag=false;

    referee(int nlead_=100, int npool_=1000, int param_len_=12, double dt_sim_=0.002, double gau_var_=0.01, double lambda_coeff_=1.0, double rs_fill_factor_=0.5, double rs_full_factor_=0.99): nlead(nlead_), npool(npool_), param_len(param_len_), dt_sim(dt_sim_), gau_var(gau_var_), lambda(lambda_coeff_*log(((double) 1e16-1)/((double) nlead-1))), rs_fill_factor(rs_fill_factor_), rs_full_factor(rs_full_factor_) {}

    referee(referee &ref_): nlead(ref_.nlead), npool(ref_.npool), param_len(ref_.param_len), dt_sim(ref_.dt_sim), gau_var(ref_.gau_var), lambda(ref_.lambda), rs_fill_factor(ref_.rs_fill_factor), rs_full_factor(ref_.rs_full_factor) {}

    ~referee();

    void alloc_records();
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
        int Frames;
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
        double l2score_min;

        /** the data index we take to be the initial conditions of the experiment */
        const int ic_index;

        ODR_struct * odr;

        race(referee &rparam,swirl_param &sp_min_,swirl_param &sp_max_,wall_list &wl_,double t_phys_,ODR_struct * odr_,int ic_index_=0);
        ~race();
        void init_race();
        void start_race(int gen_max_, bool verbose_=true);

        void make_best_swirl(char * name_);
          void make_best_swirl(const char *name_)
            {make_best_swirl((char*)name_);}
    private:
        /** The number of threads. */
        const int nt;
        /** The array AYuniform generators. */
        AYuniform * uni;
        /** A reference to the list of walls for the swirling simulation. */
        wall_list &wl;
        /** The array of random number generators. */
        AYrng ** rng;
        /** The array of proximity grids. */
        proximity_grid** const pg;

        /** Array of swirl simulations used by each thread to test particles */
        runner ** runners;
        /** count of particles in current pool beating worst leader particles */
        int pool_success_count;
        int pool_candidates;

        bool check_pool_results();
        void resample_pool();
        int collect_pool_leaders();

        void print_pool_leaders();
        void print_leaders();
        void print_leader_candidates();
        void print_best();
        void print_worst_leader();
        void print_reference_positions(int len_=10);
        void print_pool_params();
        void print_pool_records();
#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

int find_worst(record ** r, int ncap);
int find_best(record ** r, int ncap);

#endif
