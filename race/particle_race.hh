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

const int record_int_len=7;
const int record_double_len=1;

struct record
{
  bool success;
  int len;
  int Frames;

  const int global_index; // index position of this record
  int frscore; // frame score
  int gen; // generation in which this particle was generated
  int parent_gen; // generation this particle comes from
  int parent_count; // number of particle ancestors
  int parent_global_index; // index position of parent particle
  int dup_count; // number of times this particle has been duplicated

  double l2score;
  double z;

  double gau_h;
  double gau_lambda;

  double * params;

  record(int i_): global_index(i_) {}
  ~record() {}

  void init(double * params_, int len_, int Frames_, double gau_h_, double gau_lambda_);
  void reset_record(int gen_=0);
  void resample(int gen_, double * dmin_, double *dmax_, AYrng * r_);
  void duplicate(record *parent_, int gen_, double *dmin_, double *dmax_, AYrng * r_);
  void take_vals(record * rtake_);

  inline void init_leader(double * params_, int len_, int Frames_, double gau_h_, double gau_lambda_)
  {init(params_, len_, Frames_, gau_h_, gau_lambda_); reset_record();}
  inline void init_pool(double * params_, int len_, int Frames_, double gau_h_, double gau_lambda_, double * dmin_, double *dmax_, AYrng * r_)
  {init(params_, len_, Frames_, gau_h_, gau_lambda_); resample(0,dmin_,dmax_,r_);}

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
  inline bool isworse(record * rcomp_)
  {return !(isbetter(rcomp_));}

  inline double phi(int F_)
  {return ((double) frscore)/((double) F_);}

  inline double w(int F_, double lambda_)
  {return exp(lambda_*((double) frscore)/((double) F_));}

  inline double w(double lambda_z_)
  {return lambda_z_*exp(-lambda_z_*z);}

  inline double z_eval(int frscore_min_)
  {return z;}

  inline double var()
  {return gau_h*exp(-2.0*gau_lambda*((double) frscore)/((double)Frames));}

  void print_record();
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
      int run_race(double dt_sim_, double *ts_, double *xs_, double *d_ang_, record * rec_, double frscore_min_, double l2score_min_);

      void print_raw_ics();
      void print_params();
      void print_current_pos();
      void print_reset_pos() {reset_sim(pvals); print_current_pos();}
    private:
      bool is_lost(double *f_);
};

struct referee
{
    /** The maximum simulation timestep to use. */
    const double dt_sim;
    /** the max perturbation variance for duplicated particles */
    const double gau_var_high;
    /** the lowest perturbation variance for duplicated particles */
    const double gau_var_low;
    /** the commensurate exponential factor for the above two */
    const double gau_lambda;
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

    referee(int nlead_, int npool_, int param_len_, double dt_sim_, double gau_var_high_, double gau_var_low_, double lambda_coeff_, double rs_fill_factor_, double rs_full_factor_): nlead(nlead_), npool(npool_), param_len(param_len_), dt_sim(dt_sim_), gau_var_high(gau_var_high_), gau_var_low(gau_var_low_), gau_lambda(log(gau_var_high_/gau_var_low_)), lambda((lambda_coeff_)*log(((double) 1e16-1)/((double) nlead_-1))), rs_fill_factor(rs_fill_factor_), rs_full_factor(rs_full_factor_) {}

    referee(referee &ref_): nlead(ref_.nlead), npool(ref_.npool), param_len(ref_.param_len), dt_sim(ref_.dt_sim), gau_var_high(ref_.gau_var_high), gau_var_low(ref_.gau_var_low), gau_lambda(ref_.gau_lambda), lambda(ref_.lambda), rs_fill_factor(ref_.rs_fill_factor), rs_full_factor(ref_.rs_full_factor) {}

    ~referee();

    void alloc_records();

    void print_referee_params();
};

class reporter : public ODR_struct
{
  public:
    bool staged_flag = false;
    int nlead, npool, len;

    int * dup_vec;

    double *sample_weights;

    record ** leaders;

    reporter(): ODR_struct() {}
    ~reporter() {}
    void init_race(char * proc_loc_, char * rydat_dir_, char * file_name_);
      void init_race(const char *proc_loc_, const char *rydat_dir_, const char *file_name_)
        {init_race((char*)proc_loc_,(char*)rydat_dir_,(char*)file_name_);}
    void write_gen_diagnostics(int gen_count_, int leader_count_, int worst_leader_, int best_leader_);
    void close_diagnostics(int gen_count_, int leader_count_, int worst_leader_, int best_leader_);
    ODR_struct * spawn_swirlODR(char *name_);

};

class race : public referee
{
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
        /** frame score associated w/ best particle on leaderboard  */
        int frscore_best;
        /** error score associated w/ best particle on leaderboard  */
        double l2score_best;

        /** the data index we take to be the initial conditions of the experiment */
        const int ic_index;

        reporter * odr;

        race(referee &rparam,swirl_param &sp_min_,swirl_param &sp_max_,wall_list &wl_,double t_phys_,reporter * odr_,int ic_index_=0);
        ~race();
        void init_race();
        void start_race(int gen_max_, bool verbose_=true);

        void make_best_swirl(char * name_);
          void make_best_swirl(const char *name_)
            {make_best_swirl((char*)name_);}
    private:
        /** The number of threads. */
        const int nt;
        /** A reference to the list of walls for the swirling simulation. */
        wall_list &wl;
        /** The array of random number generators. */
        AYrng ** rng;
        /** The array of proximity grids. */
        proximity_grid** const pg;

        int *dup_vec, **dup_mat;

        /** Array of swirl simulations used by each thread to test particles */
        runner ** runners;
        /** count of particles in current pool beating worst leader particles */
        int pool_success_count;
        int pool_candidates;
        int best_leader, worst_leader;

        bool debugging_flag=true;
        bool z_weight_flag=true;


        void stage_diagnostics();

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
