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

  const int id,
            fullbuf_len;

  int Frames,
      nbeads;

  size_t  obuf_len,
          ibuf_len;

  private:
    const bool noise_data;
    const double noise_sigma;

    void read_fisml(char * ibuf_, int * Frames_);
};

const int record_int_len=5;
const int record_double_len=2;
struct record
{
  record(int rid_, int nbeads_, int Frames_, int len_, int * int_chunk_, double * double_chunk_, double * params_): rid(rid_),
  nbeads(nbeads_), Frames(Frames_), len(len_),
  int_params(&gen), double_params(&residual),
  params(params_) {}
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

  double  residual,
          w;
};

class thread_worker: public swirl
{
    public:

      thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int nbeads_, int thread_id_, int param_len_, int Frames_, int nlead_, int npool_, double dt_sim_, double t_phys_, double *ts_, double *xs_, double *d_ang_, double *comega_s_);
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
              * const comega_s;
};

struct MH_params
{
  MH_params(int nlead_, int npool_, int param_len_, double dt_sim_, double t_phys_, double sigma_): nlead(nlead_), npool(npool_), param_len(param_len_), dt_sim(dt_sim_), t_phys(t_phys_), sigma(sigma_) {}
  MH_params(MH_params &par_): nlead(par_.nlead), npool(par_.npool), param_len(par_.param_len), dt_sim(par_.dt_sim), t_phys(par_.t_phys), sigma(par_.sigma) {}

  ~MH_params();

  const int nlead, // number of leaders we store for review
            npool, // number of parameter sets we evaluate at a time
            param_len; // length of parameter vector

  const double  dt_sim, // The maximum simulation timestep to use
                t_phys, // simulation time scale
                sigma; // expected standard deviation of noise on reference data
};

class MH_trainer : public MH_params
{
  public:

    MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_, double t_wheels_=0.0);

    ~MH_trainer();

    swirl_param sp_min, // lower boundary of parameter space U
                sp_max; // upper boundary of parameter space U

    const int nbeads, // the total number of beads
              Frames; // the total number of data frames

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
            wres; // current worst leader residual

    double * const  ts, // wall time of observed data
           * const  xs, // 2D observed position data
           * const  d_ang, // observed angular position of dish
           * const  comega_s; // observed average rotational speed of dish

    int ** record_int_chunk;

    double  **record_double_chunk,
            **param_chunk;

    record  ** const records,
            ** const leaders, ** const pool,
            ** leader_board, ** candidates;

    protected:

      const int nt; // number of worker threads

      wall_list &wl; // a reference to the list of walls for the swirling simulation.
      proximity_grid** const pg; // array of proximity grids.
      AYrng ** rng; // random number generators

      thread_worker ** const workers; // team of threads for solving and sampling

#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

int find_worst_record(record ** r, int ncap);
int find_best_record(record ** r, int ncap);

#endif
