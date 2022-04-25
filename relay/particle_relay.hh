#ifndef PARTICLE_RELAY_HH
#define PARTICLE_RELAY_HH

#ifdef _OPENMP
#include "omp.h"
#endif

const int record_int_len=7;
const int record_double_len=2;

struct record
{
  bool success;
  int len,Frames;

  const int beads;

  const int global_index;
        int gen; // generation in which this particle was generated
        int parent_gen; // generation this particle comes from
        int parent_count; // number of particle ancestors
        int parent_global_index; // index position of parent particle
        int dup_count; // number of times this particle has been duplicated

  double  residual,
          event_residual,
          weight;

  int * const event_positions;

  double * const params;
  double * const residual_data;
  double * const alpha_data;

  record(int i_, double * params_): global_index(i_), params(params_) {}
  ~record() {}

  void reset_record(int gen_);
  void resample(int gen_, double * dmin_, double *dmax_, AYrng * r_);
  void duplicate(record *parent_, int gen_, double *dmin_, double *dmax_, AYrng * r_, double * var_);
  void take_vals(record * rtake_);

  inline void init(int len_, int Frames_)
  {len=len_; Frames=Frames_;}
  inline void init_leader(int len_, int Frames_)
  {init(len_, Frames_); reset_record(0);}
  inline void init_pool(int len_, int Frames_, double * dmin_, double *dmax_, AYrng * r_)
  {init(len_, Frames_); resample(0,dmin_,dmax_,r_);}
  inline bool check_success(double residual_, double residual_worst_)
  {residual = residual_; return success=residual<residual_worst_;}
  inline bool isbetter(record * rcomp_)
  {return residual<rcomp_->residual;}
  inline bool isworse(record * rcomp_)
  {return !(isbetter(rcomp_));}
  inline double w(double c_, double b_)
  {return weight = exp(c_-0.5*(residual/b_)*(residual/b_));}
};

class runner: public swirl
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
      const double alpha_tol=100.0;

      int frame;

      int *lead_dup_count, // length = nlead
          *frame_kill_count; // length = Frames

      double  pos_res_acc,
              *frame_res_data, // [pos_res|max_err|pos_res_dead|pos_res_alive]
              *param_mean;

      runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, double tol_, double t_phys_, double dt_sim_, double t_wheels_, double *ts_, double *xs_, double *d_ang_, int ic_index_, int nlead_, int npool_);
      ~runner();

      void reset_sim(double *ptest_);

      int run_relay(record * gra_, int frscore_min_, double residual_min_);
      void detect_events(record* rec_, int start_, int end_);

      void reset_diagnostics();
      void consolidate_diagnostics();
      void update_diagnostics(int *frame_kill_count_, double *gen_frame_res_data_);

    private:
      bool  dead;

      int   n_test,
            *kill_frames,
            **event_frame_count;

      double  *pvals,
              *ts, *xs, *d_ang,
              **pos_res, **alpha_INTpos_res,
              **INTpos_res;

      void compute_error(double *f_, double dur_);


      inline double advance_runner()
      {
        double dur=(ts[frame]-ts[frame-1])/t_phys, ctheta=d_ang[frame-1], comega=d_ang[frame]-d_ang[frame-1];
        if(comega>M_PI) comega-=2*M_PI; else if(comega<-M_PI) comega+=2*M_PI; comega/=dur;
        advance(dur, ctheta, comega, dt_sim);
        return dur;
      }
      inline void clear_event_data()
      {for (int i = 0; i < n*Frames; i++) event_frame_count[0][i] = 0;}
      inline double alpha_comp(double *a_, double t_m1, double t_p1)
      {
        double alpha_val = log(a_[0]/a_[2])/log(t_p1/t_m1);
        a_[2] = a_[1]; a_[1] = a_[0]; // shift the integral history along
        return alpha_val;
      }
};

struct referee
{
  /** The maximum simulation timestep to use. */
  const double dt_sim;
  const double noise_tol; // of same order as standard deviation
  const double max_weight_ceiling = 1e13;

  const int nlead;
  const int npool;
  const int nA;
  const int param_len;

  record  **records,
          **leaders, **pool,
          **leader_board, **candidates;

  // memory chunks to store parameters associated with particles
  double  **param_chunk;

  int *lead_dup_count,
      *global_event_frame_count;

  double  *sample_weights, // length = nlead
          *gen_frame_res_data, // [pos_res|max_err|pos_res_dead|pos_res_alive]
          *lead_par_w_mean,
          *lead_par_w_var;

  bool alloc_flag=false;

  referee(int nlead_, int npool_, int nA_, int param_len_, double dt_sim_, double t_wheels_): nlead(nlead_), npool(npool_), nA(nA_), param_len(param_len_), dt_sim(dt_sim_), t_wheels(t_wheels_) {}
  referee(referee &ref_): nlead(ref_.nlead), npool(ref_.npool), nA(ref_.nA), param_len(ref_.param_len), dt_sim(ref_.dt_sim), t_wheels(ref_.t_wheels) {}

  ~referee();

  void alloc_records(int nt_, int Frames_);
};

class reporter : public ODR_struct
{
  public:
    bool staged_flag = false;
    int relay_id;
    int nlead, npool, nA, param_len;

    int *lead_dup_count,
        *frame_kill_count;

    double  *sample_weights,
            *gen_frame_res_data,
            *gen_param_mean,
            *gen_param_var;

    record ** leaders;

    reporter() : ODR_struct() {}
    ~reporter() {}
    void init_relay(char * proc_loc_, char * rydat_dir_, char * file_name_, int relay_id_);
      void init_relay(const char *proc_loc_, const char *rydat_dir_, const char *file_name_, int relay_id_)
        {init_relay((char*)proc_loc_,(char*)rydat_dir_,(char*)file_name_, relay_id_);}

    void write_gen_diagnostics(int gen_count_, int leader_count_, int worst_leader_, int best_leader_, int dup_count_, int resample_count_, int dup_unique_, int repl_count_, double t_wheels_, double min_res_);
    void close_diagnostics(int gen_count_, int worst_leader_, int best_leader_, double t_wheels_, double min_res_);
    ODR_struct * spawn_swirlODR(char *name_);
};

class relay : public referee
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

    int leader_count,
        pool_success_count,
        pool_candidates,
        gen_count;

    double  residual_worst,
            residual_best,
            gau_scale_sqrt,
            max_weight_factor;

    reporter * rep;

    relay(referee &ref_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_,double t_phys_,reporter * ped_,int ic_index_=0);
    ~relay();

    void init_relay();
    void start_relay(int gen_max_, bool verbose_=true);

    void make_best_swirl(char * name_);
      void make_best_swirl(const char *name_)
        {make_best_swirl((char*)name_);}

    private:
      /** The number of threads. */
      const int nt;

      int best_leader,
          worst_leader,
          dup_count,
          resample_count,
          dup_unique,
          repl_count,
          latest_event,
          event_observations;

      bool debugging_flag=true;

      /** A reference to the list of walls for the swirling simulation. */
      wall_list &wl;
      /** The array of random number generators. */
      AYrng ** rng;
      /** The array of proximity grids. */
      proximity_grid** const pg;

      runner ** runners;

      void stage_diagnostics();

      void learn_first_leg(int gen_max_, bool verbose_);

      bool check_pool_results();
      int collect_candidates();

      double compute_leader_statistics();
      void resample_pool();

      inline void clear_global_event_data()
      {for (int i = 0; i < n*Frames; i++) global_event_frame_count[0][i] = 0;}

#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

int find_worst_record(record ** r, int ncap);
int find_best_record(record ** r, int ncap);

#endif
