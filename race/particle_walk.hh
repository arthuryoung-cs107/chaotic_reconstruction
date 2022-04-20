#ifndef PARTICLE_WALK_HH
#define PARTICLE_WALK_HH

#include "particle_race.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

const int grade_int_len=8;
const int grade_double_len=1;

struct grade
{
  bool success;
  int len,Frames;

  const int global_index;
        int school; // refers to the class A particle lineage
        int frscore;
        int gen; // generation in which this particle was generated
        int parent_gen; // generation this particle comes from
        int parent_count; // number of particle ancestors
        int parent_global_index; // index position of parent particle
        int dup_count; // number of times this particle has been duplicated

  double  l2score;
  double  z;

  double * const params;

  grade(int i_, double * params_): global_index(i_), params(params_) {}
  ~grade() {}

  void reset_grade(int gen_);
  void resample(int gen_, double * dmin_, double *dmax_, AYrng * r_);
  void duplicate(grade *parent_, int gen_, double *dmin_, double *dmax_, AYrng * r_);
  void take_vals(grade * rtake_);

  inline void init(int len_, int Frames_)
  {len=len_; Frames=Frames_;}
  inline void init_leader(double * params_, int len_, int Frames_)
  {init(len_, Frames_); reset_grade(0);}
  inline void init_pool(double * params_, int len_, int Frames_, double * dmin_, double *dmax_, AYrng * r_)
  {init(len_, Frames_); resample(0,dmin_,dmax_,r_);}

  inline bool check_success(int frscore_, double l2score_, int frmin_, double l2min_)
  {
    frscore = frscore_; l2score = l2score_;
    return success=(l2score<l2min_)?true:false; // conditioning on average frame error
  }
  inline bool isbetter(grade * rcomp_)
  {return l2score<rcomp_->l2score;}

  inline bool isworse(grade * rcomp_)
  {return !(isbetter(rcomp_));}

  inline double w(double lambda_z_)
  {return lambda_z_*exp(-lambda_z_*z);}

  inline double z_eval(int frscore_min_)
  {return (z=(1.0/((double)(frscore-frscore_min_+1)))*(1.0/((double)frscore))*l2score);}

  inline double var()
  {return gau_h*exp(-2.0*gau_lambda*((double) frscore)/((double)Frames));}

};

class walker: public swirl
{
    public:
      const int thread_id;
      const int param_len;
      const int Frames;
      const int nlead;
      const int npool;

      const double tol;
      const double t_phys;
      const double dt_sim;
      const double t_wheels; // rycroft's default is 0.012

      int frame;

      int *lead_dup_count, // length = nlead
          *frame_kill_count; // length = Frames

      double  pos_err_acc,
              *frame_err_data; // [pos_err|max_err|pos_err_dead|pos_err_alive]

      walker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, double tol_, double t_phys_, double dt_sim_, double t_wheels_, double *ts_, double *xs_, double *d_ang_, int ic_index_, int nlead_, int npool_);
      ~walker();

      int start_walking(grade * gra_, int frscore_min_, double l2score_min_);

      void reset_diagnostics()
      {n_test=0; for (int i = 0; i < 4*Frames; i++) frame_err_data[i] = 0.0; for (int i = 0; i < Frames; i++) frame_kill_count[i] = 0;}
      void consolidate_diagnostics();
      void update_diagnostics(int *frame_kill_count_, double *mean_frame_err_data_);

    private:
      bool dead;
      int frame_kill,
          n_test;

      double  *pvals,
              *ts, *xs, *d_ang,
              t0, *x0, ctheta0, t0_raw;

      void compute_error(double *f_, double dur_);
};

struct guide
{
  /** The maximum simulation timestep to use. */
  const double dt_sim;
        double t_wheels;

  const int nlead;
  const int npool;
  const int nA;
  const int param_len;

  grade **grades,
        **leaders, **pool,
        **classA, **classB,
        **leader_board, **candidates;

  // memory chunks to store parameters associated with particles
  double  **param_chunk;

  int *lead_dup_count, // length = nlead
      *frame_kill_count; // length = Frames

  double  *sample_weights, // length = nlead
          *mean_frame_err_data; // [pos_err|max_err|pos_err_dead|pos_err_alive]

  bool alloc_flag=false;

  guide(int nlead_, int npool_, int nA_, int param_len_, double dt_sim_, double t_wheels_): nlead(nlead_), npool(npool_), nA(nA_), param_len(param_len_), dt_sim(dt_sim_), t_wheels(t_wheels_) {}
  guide(guide &gui_): nparent(gui_.nparent), nlead(gui_.nlead), npool(gui_.npool), param_len(gui_.param_len), dt_sim(gui_.dt_sim), t_wheels(gui_.t_wheels) {}

  ~guide();

  void alloc_grades(int nt_, int Frames_);
};

class pedestrian : public ped_struct
{
  public:
    bool staged_flag = false;
    int nlead, npool, nA, len, walk_id;



    double * sample_weights;

    grade ** leaders;

    pedestrian() : ped_struct() {}
    ~pedestrian() {}
    void init_walk(char * proc_loc_, char * rydat_dir_, char * file_name_, int walk_id_);
      void init_walk(const char *proc_loc_, const char *rydat_dir_, const char *file_name_, int walk_id_)
        {init_walk((char*)proc_loc_,(char*)rydat_dir_,(char*)file_name_, walk_id_);}

    void write_gen_diagnostics(int gen_count_, int leader_count_, int worst_leader_, int best_leader_, double t_wheels_);
    void close_diagnostics(int gen_count_, int leader_count_, int worst_leader_, int best_leader_);
    ODR_struct * spawn_swirlODR(char *name_);
};

class walk : public guide
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
    /** The time points of the snapshots. */
    double* ts;
    /** The bead positions at the snapshots. */
    double* xs;
    /** The dish angle data at the snapshots. */
    double* d_ang;

    /** the data index we take to be the initial conditions of the experiment */
    const int ic_index;

    pedestrian * ped;

    walk(guide &gui_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_,double t_phys_,pedestrian * ped_,int ic_index_=0);
    ~walk();

    void init_walk();
    void start_walk(int gen_max_, bool verbose_=true);

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

      walker ** walkers;

      void stage_diagnostics();

      void train_classA(bool verbose_);


#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

int find_worst_grade(record ** r, int ncap);
int find_best_grade(record ** r, int ncap);

#endif
