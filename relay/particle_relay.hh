#ifndef PARTICLE_RELAY_HH
#define PARTICLE_RELAY_HH

#include "swirl.hh"
#include "RYdat2AYdat.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

const int record_int_len=6;
const int record_double_len=7;

const int record_int_chunk_count = 1;
const int record_double_chunk_count = 3;

const double log_P = log(1e14);
const double par_acc_P = 1e14;
const double log_C = log(1.0/sqrt(2.0*M_PI));

struct record
{
  bool success;

  const int beads;
  const int Frames;
  const int len;
  const int global_index;
        int gen; // generation in which this particle was generated
        int parent_gen; // generation this particle comes from
        int parent_count; // number of particle ancestors
        int parent_global_index; // index position of parent particle
        int dup_count; // number of times this particle has been duplicated

  double residual;
  double root_residual;
  double net_residual;
  double net_smooth_residual;
  double weight;
  double px;
  double zeta;

  int * int_params;
  double * double_params;

  int * const int_chunk;
  double * const double_chunk;

  int * const event_positions;

  double * const smooth_residual;
  double * const stiff_residual;
  double * const alpha_data;

  double * const params;

  bool smooth_training=true;

  record(int beads_, int Frames_, int len_, int i_, int * int_chunk_, double * params_, double * double_chunk_): beads(beads_), Frames(Frames_), len(len_), global_index(i_), int_chunk(int_chunk_), double_chunk(double_chunk_), event_positions(int_chunk_), smooth_residual(double_chunk_), stiff_residual(double_chunk_+beads_), alpha_data(double_chunk_+2*beads_), params(params_)
  {int_params=&(gen); double_params=&(residual);}

  ~record() {}

  void reset_record(int gen_, int p_gen_=-1, int p_count_=0, int p_gi_=-1);
  void resample(int gen_, double * dmin_, double *dmax_, AYrng * r_);
  void duplicate(record *parent_, int gen_, double *dmin_, double *dmax_, AYrng * r_, double * var_);
  void take_vals(record * rtake_);

  inline void init_leader()
  {reset_record(0);}
  inline void init_pool(double * dmin_, double *dmax_, AYrng * r_)
  {resample(0,dmin_,dmax_,r_);}
  inline bool check_success(double net_residual_, double net_smooth_residual_, double residual_worst_)
  {
    net_residual=net_residual_; net_smooth_residual=net_smooth_residual_;
    residual = (smooth_training)?net_smooth_residual:net_residual;
    root_residual=sqrt(residual); return success=(residual<residual_worst_);
  }
  inline bool isbetter(record * rcomp_)
  {return residual<rcomp_->residual;}
  inline bool isworse(record * rcomp_)
  {return !(isbetter(rcomp_));}

  // being as careful as possible for floating point precision
  inline double w(double lambda_, double beta_, double tau_)
  {
    double itau=1.0/tau_;
    zeta = (0.5*(root_residual*itau))*(root_residual*itau);
    px = exp(log_C+(log(itau)-zeta));
    return weight = lambda_*exp(log_P+(beta_-zeta));
  }

  inline void init_training(bool smooth_training_)
  {
    residual=(smooth_training=smooth_training_)?residual=net_smooth_residual:net_residual;
    root_residual=sqrt(residual);
  }
};

struct events
{
  const int beads;
  const int Frames;

  int ** global_event_frame_count;
  int * const event_frames;

  int * const bead_order;
  int * const event_sorted;

  events(int Frames_, int beads_, int ** global_event_frame_count_, int *event_frames_in_): beads(beads_), Frames(Frames_),
  global_event_frame_count(global_event_frame_count_), event_frames(event_frames_in_),
  bead_order(new int[beads_]), event_sorted((new int[beads_]))
  {}
  ~events() {delete [] bead_order; delete [] event_sorted;}

  void define_relay_event_block(int event_block_id_, int * obs_vec, double * tau_vec, int * early_late_events, double tau_coeff);

  inline int earliest(int *index, int start_index_)
  {
    int out = event_sorted[start_index_]; *index=start_index_;
    for (int i = start_index_+1; i < beads; i++) if (event_sorted[i]<out) out=event_sorted[(*index)=i];
    return out;
  }
  inline void earliest_recursive(int i_next)
  {
    int i_it;
    int early_it = earliest(&i_it,i_next), early_temp = event_sorted[i_next], i_temp = bead_order[i_next];
    event_sorted[i_next]=event_sorted[i_it]; bead_order[i_next]=bead_order[i_it];
    event_sorted[i_it]=early_temp; bead_order[i_it]=i_temp;
    if (i_next<beads-1) earliest_recursive(i_next+1);
  }
};

class runner: public swirl
{
    public:
      const int thread_id;
      const int param_len;
      const int Frames;
      const int nlead;
      const int npool;

      const double t_phys;
      const double dt_sim;
      const double alpha_tol;

      int *lead_dup_count,
          **event_frame_count;

      double  *pvals, *ts, *xs, *d_ang, *comega_s,
              *param_acc,
              **pos_res, **alpha_INTpos_res,
              **INTpos_res,
              *sim_pos,
              *bead_event_res_mat;

      runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, int nlead_, int npool_, double t_phys_, double dt_sim_, double alpha_tol_, double *ts_, double *xs_, double *d_ang_, double *comega_s);
      ~runner();

      void reset_sim(double *ptest_, double t0_, double ctheta0_, double comega0_, double *x0_);

      int run_relay(record * rec_, int start_, int earliest_, int latest_, int * event_frames_ordered_, int * bead_order, double residual_worst_);
      int start_detection(int start_, double * params_, double *t_history_, int * event_frames, double * smooth_residual, double * net_residual_);
      void detect_events(record* rec_, int start_, int end_);

      inline void clear_event_data()
      {for (int i = 0; i < n*Frames; i++) event_frame_count[0][i] = 0;}

      // debugging stuff
      double *TEST_pos, *TEST_pos_res, *TEST_alpha_INTpos_res, *TEST_INTpos_res;
      int start_test_detection(int start_, double * params_, double *t_history_, int * event_frames, double * smooth_residual, double * net_residual_);
      void test_detection(record * rec_, int start_, int end_, int i_);

      void init_test_run(int Frame_end_);
      void free_test_buffers();

    private:

      inline double alpha_comp(double *a_, double t_p1, double t_m1)
      {return log(a_[0]/a_[2])/log(t_p1/t_m1);}
      inline double alpha_comp(double a_p1, double a_m1, double t_p1, double t_m1)
      {return log(a_p1/a_m1)/log(t_p1/t_m1);}

};

struct referee
{
  const int nlead;
  const int npool;
  const int param_len;

  /** The maximum simulation timestep to use. */
  const double dt_sim;
  const double noise_tol; // projected standard deviation of noised data
  const double alpha_tol; // threshold for event detection
  const double rs_full_factor; // ratio of duplicated particles to resampled particles per generation
  const double data_scale;

  record  **records,
          **leaders, **pool,
          **leader_board, **candidates;

  events * ev;

  int **record_int_chunk,
      **global_event_frame_count,
      *lead_dup_count,
      *event_frames;

      // conservative estimate for pre-collision frame

  double  **param_chunk,
          **record_double_chunk,
          *sample_weights,
          *lead_par_w_mean,
          *lead_par_w_var;

  bool alloc_flag=false;

  referee(int nlead_, int npool_, int param_len_, double dt_sim_, double noise_tol_, double alpha_tol_, double rs_full_factor_, double data_scale_): nlead(nlead_), npool(npool_), param_len(param_len_), dt_sim(dt_sim_), noise_tol(noise_tol_), alpha_tol(alpha_tol_), rs_full_factor(rs_full_factor_), data_scale(data_scale_) {}
  referee(referee &ref_): nlead(ref_.nlead), npool(ref_.npool), param_len(ref_.param_len), dt_sim(ref_.dt_sim), noise_tol(ref_.noise_tol), alpha_tol(ref_.alpha_tol) ,rs_full_factor(ref_.rs_full_factor), data_scale(ref_.data_scale) {}

  ~referee();

  void alloc_records(int nt_, int Frames_, int beads_);
};



class reporter : public ODR_struct
{
  public:
    bool staged_flag = false;
    int relay_id;

    int nlead,
        npool,
        param_len,
        beads;

    int *lead_dup_count,
        *global_event_frame_count,
        *event_frames,
        *gen_int_vec,
        *event_int_vec;

    double  *sample_weights,
            *lead_par_w_mean,
            *lead_par_w_var,
            *gen_double_vec,
            *event_double_vec;

    record ** leaders;
    record ** pool;

    reporter() : ODR_struct() {}
    ~reporter() {}
    void init_relay(char * proc_loc_, char * rydat_dir_, char * file_name_, int relay_id_, bool noise_data_, double noise_tol_in_);
      void init_relay(const char *proc_loc_, const char *rydat_dir_, const char *file_name_, int relay_id_, bool noise_data_, double noise_tol_in_)
        {init_relay((char*)proc_loc_,(char*)rydat_dir_,(char*)file_name_, relay_id_, noise_data_, noise_tol_in_);}


    void load_relay(double *ts_, double *xs_, double *d_ang_);
    void write_startup_diagnostics(int *header_, int *int_params_, double *double_params_);
    void write_event_diagnostics(int event_);
    void write_gen_diagnostics(int gen_count_, int leader_count_);
    void close_diagnostics(int gen_count_);

  private:

    bool noise_data;
    double noise_tol_in;
};

const int gen_int_len=10;
const int gen_double_len=9;
const int event_int_len=5;
const int event_double_len=5;

class relay : public referee
{
  public:
    /** The minimum parameters. */
    swirl_param sp_min;
    /** The maximum parameters. */
    swirl_param sp_max;
    /** The simulation time in seconds. */
    const double t_phys;
    /** The total number of beads. */
    const int n;
    /** The total number of snapshots. */
    const int Frames;
    /** The time points of the snapshots. */
    double* ts;
    /** The bead positions at the snapshots. */
    double* xs;
    /** The dish angle data at the snapshots. */
    double* d_ang;
    /** The dish average rotational speed. */
    double* comega_s;

    int gen_count;
    int leader_count;
    int pool_success_count;
    int pool_candidates;
    int best_leader;
    int worst_leader;
    int repl_count;
    int dup_count;
    int dup_unique;
    int resample_count;

    int earliest_event;
    int latest_event;
    int net_observations;
    int smooth_obs;
    int stiff_obs;

    double residual_best;
    double residual_worst;
    double root_res_best;
    double root_res_worst;
    double pos_res_global;
    double lambda;
    double beta;
    double w_sum;
    double w_max;
    double t_wheel;

    double tau;
    double tau_sqr;
    double tau_full;
    double tau_smooth;
    double tau_stiff;

    reporter * rep;

    relay(referee &ref_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_,double t_phys_,reporter * rep_);
    ~relay();

    void init_relay();
    void start_relay(int gen_max_, bool verbose_=true);

    void make_best_swirl(char * name_);
      void make_best_swirl(const char *name_)
        {make_best_swirl((char*)name_);}

    protected:
      /** The number of threads. */
      const int nt;

      bool  debugging_flag=true;

      /** A reference to the list of walls for the swirling simulation. */
      wall_list &wl;
      /** The array of random number generators. */
      AYrng ** rng;
      /** The array of proximity grids. */
      proximity_grid** const pg;

      runner ** runners;

      double * param_acc_factor; // for overflow resolution

      void stage_diagnostics(int gen_max_);

      void find_events(int min_frame_, int latest_frame_);
      int train_event_block(int event_block, int gen_max_, double tol_leeway_=0.0);

      void check_gen0();

      bool check_pool_results(double tol_leeway_=0.0);
      int collect_candidates();

      void compute_leader_statistics();
      void resample_pool();

      inline void clear_global_event_data()
      {for (int i = 0; i < n*Frames; i++) global_event_frame_count[0][i] = 0;}

#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

struct medic
{
  const int Frame_end;
  const int beads;

  runner * rt;

  char* test_directory;
  size_t buf_end;

  medic(int Frame_end_, runner * rt_, char * test_directory_): Frame_end(Frame_end_), beads(rt_->n), rt(rt_), test_directory(new char[strlen(test_directory_)+100])
  {strcpy(test_directory, test_directory_); buf_end = strlen(test_directory);}
  ~medic()
  {delete [] test_directory;}

  void write_test_run(record * rec_, int i_);
};

class doctor : public relay
{
  public:

    char * test_buffer=NULL;

    doctor(referee &ref_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_,double t_phys_,reporter * rep_): relay(ref_, sp_min_, sp_max_, wl_,t_phys_,rep_) {}
    ~doctor() {if (test_buffer!=NULL) delete [] test_buffer;}

    void init_test(int test_id_, int test_relay_id_);
    void test_run(int Frame_end_);


    double *TEST_ref_pos=NULL;
};

int find_worst_record(record ** r, int ncap);
int find_best_record(record ** r, int ncap);

#endif
