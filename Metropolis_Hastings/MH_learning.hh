#ifndef MH_LEARNING_HH
#define MH_LEARNING_HH

#include "swirl.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

struct MH_io
{
  MH_io(char * proc_loc_, char * test_dir_, char * data_name_, int id_, bool noise_data_=false, double noise_sigma_=0.0);
  ~MH_io();

  char  * const proc_loc,
        * const test_dir,
        * const data_name,
        * const obuf,
        * const ibuf;

  const size_t fullbuf_len;
  size_t  obuf_len,
          ibuf_len;

  const int id;

  int Frames,
      nbeads;

  void load_reference(double *ts_, double *xs_, double *d_ang_, double * comega_s_, double t_phys_);

  private:
    const bool noise_data;
    const double noise_sigma;

    void read_fisml(char * ibuf_);
};

const int record_int_len=5;
const int record_double_len=2;
struct record
{
  record(int rid_, int nbeads_, int Frames_, int len_, int * int_chunk_, double * double_chunk_, double * u_): rid(rid_),
  nbeads(nbeads_), Frames(Frames_), len(len_),
  iparams(&gen), dparams(&residual),
  u(u_) {}
  ~record() {}

  bool success;

  const int rid,
            nbeads,
            Frames,
            len;

  int gen, // generation in which this particle was generated
      parent_gen, // generation this particle comes from
      parent_count, // number of particle ancestors
      parent_rid, // index position of parent particle
      dup_count; // number of times this particle has been duplicated

  int * iparams;

  double  residual,
          w;

  double  * const dparams,
          * const u;
};

class thread_worker: public swirl
{
    public:

      thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int thread_id_, int ichunk_len_, int dchunk_len_, int param_len_, int nbeads_, int Frames_, int nlead_, int npool_, double dt_sim_, double t_phys_, double *ts_, double *xs_, double *d_ang_, double *comega_s_);
      ~thread_worker();

    protected:

      const int thread_id,
                ichunk_len,
                dchunk_len,
                param_len,
                nbeads,
                Frames,
                nlead,
                npool;

      int * const ichunk;

      const double  dt_sim,
                    t_phys;

      double  * const dchunk,
              * const ts,
              * const xs,
              * const d_ang,
              * const comega_s,
              * pvals;
};

struct MH_params
{
  MH_params(MH_io *io_, int nlead_, int npool_, int param_len_, double dt_sim_, double t_phys_, double sigma_): io(io_),
  nbeads(io_->nbeads), Frames(io_->Frames),
  nlead(nlead_), npool(npool_), param_len(param_len_),
  dt_sim(dt_sim_), t_phys(t_phys_), sigma(sigma_) {}
  MH_params(MH_params &par_): io(par_.io),
  nbeads(par_.nbeads), Frames(par_.Frames),
  nlead(par_.nlead), npool(par_.npool), param_len(par_.param_len),
  dt_sim(par_.dt_sim), t_phys(par_.t_phys), sigma(par_.sigma) {}

  ~MH_params() {}

  const int nbeads, // number of beads
            Frames, // number of observation frames
            nlead, // number of leaders we store for review
            npool, // number of parameter sets we evaluate at a time
            param_len; // length of parameter vector


  const double  dt_sim, // The maximum simulation timestep to use
                t_phys, // simulation time scale
                sigma; // expected standard deviation of noise on reference data

  MH_io * const io;
};

struct MH_train_inputs
{
  MH_train_inputs(MH_params *par_, swirl_param *sp_min_, swirl_param *sp_max_, wall_list *wl_, double t_wheels_=0.0): par(par_) , sp_min(sp_min_), sp_max(sp_max_), wl(wl_), t_wheels(t_wheels_) {}
  ~MH_train_inputs() {}

  MH_params * const par;
  swirl_param * const sp_min,
              * const sp_max;
  wall_list * const wl;
  const double t_wheels;
};

class MH_trainer : public MH_params
{
  public:

    MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_, double t_wheels_=0.0);
    MH_trainer(MH_train_inputs &mhti): MH_trainer(*(mhti.par), *(mhti.sp_min), *(mhti.sp_max), *(mhti.wl), mhti.t_wheels) {}
    ~MH_trainer();

    swirl_param sp_min, // lower boundary of parameter space U
                sp_max; // upper boundary of parameter space U

    int leader_count, // current number of leaders
        nsuccess, // number of parameters better than current worst leader in recent trial
        ncandidates, // number of candidates to compare to current leaders
        bleader_rid, // index of current best leader
        wleader_rid, // index of current worst leader
        nreplace, // number of leader replacements following evaluation of recent trial
        ndup, // total number of leader duplication and perturbations from recent resampling
        ndup_unique, // total number of unique duplications from recent resampling
        nredraw; // total number of trial particles drawn from proposal distribution in recent resampling

    double  rho2, // current expected residual
            bres, // current best leader residual
            wres, // current worst leader residual
            t_wheels; // current training wheels

    double * const  ts, // wall time of observed data
           * const  xs, // 2D observed position data
           * const  d_ang, // observed angular position of dish
           * const  comega_s, // observed average rotational speed of dish
           ** const u_chunk; // space for parameters

    protected:

      const int nt; // number of worker threads

      wall_list &wl; // a reference to the list of walls for the swirling simulation.
      proximity_grid ** const pg; // array of proximity grids.
      AYrng ** rng; // random number generators

#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

class MH_doctor : public MH_trainer
{
  public:

    MH_doctor(MH_train_inputs &mhti);
    ~MH_doctor();

  private:

};

int find_worst_record(record ** r, int ncap);
int find_best_record(record ** r, int ncap);

#endif
