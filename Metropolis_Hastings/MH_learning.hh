#ifndef MH_LEARNING_HH
#define MH_LEARNING_HH

#include "MH_auxiliary.hh"
#include "MH_tools.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

struct record: public record_struct
{
  record(record_struct &rec_, int rid_, double * u_): record_struct(rec_), rid(rid_), u(u_) {}
  ~record() {}

  bool success;

  const int rid;

  double * const u;

  virtual int isworse(record * r_) = 0;
  virtual int isbetter(record * r_) = 0;
};

class thread_worker: public swirl, public thread_worker_struct
{
    public:

      thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, thread_worker_struct &tws, int thread_id_);
      ~thread_worker();

    protected:

      const int thread_id;

      int * const ichunk;

      double  * const u,
              * const dchunk;

};

class MH_trainer : public MH_params
{
  public:

    MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_);
    MH_trainer(MH_train_struct &mhti): MH_trainer(*(mhti.par), *(mhti.sp_min), *(mhti.sp_max), *(mhti.wl)) {}
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
            wres; // current worst leader residual

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
      inline int get_nt() {return omp_get_max_threads();}
#else
      inline int thread_num() {return 0;}
      inline int get_nt() {return 1;}
#endif
};

const int basic_record_ilen=5;
const int basic_record_dlen=2;

struct basic_record: public record
{
  basic_record(record_struct &rec_, int rid_, double * u_, int * int_chunk_, double * double_chunk_): record(rec_, rid_, u_), iparams(&gen), dparams(&r2) {}
  ~basic_record();

  int gen, // generation in which this particle was generated
      parent_gen, // generation this particle comes from
      parent_count, // number of particle ancestors
      parent_rid, // index position of parent particle
      dup_count; // number of times this particle has been duplicated

  int * const iparams;

  double  r2,
          w;

  double * const dparams;
};

class basic_thread_worker: public thread_worker, public event_detector
{
    public:

      basic_thread_worker(thread_work_struct &tws_, int thread_id_, double alpha_tol_): thread_worker(tws_, thread_id_), event_detector(alpha_tol_) {}
      ~basic_thread_worker() {}

    protected:

      int * const lead_dup_count, // lead_dup_count = ichunk+0
          ** const event_frame_count; // event_frame_count[0] = ichunk + nlead

      double  ** const r2_bead,
              ** const INTr2_bead,
              ** const alpha_bead,
              * const p_sim,
              * const r2_bead_mat;

      virtual int get_ichunk_len()  {return nlead + // lead_dup_count
                                            (nbeads*Frames); // event_frame_count
                                    }

      virtual int get_dchunk_len()  {return (nbeads*Frames) + // r2_bead
                                            (nbeads*3) + // INTr2_bead
                                            (nbeads*Frames) + // alpha_bead
                                            (2*nbeads*Frames) + // p_sim
                                            (nbeads*nbeads); // r2_bead_mat
                                    }
};

class basic_MH_trainer: public MH_trainer, public gaussian_likelihood
{
    public:

      basic_MH_trainer(MH_train_struct &mhts_, double t_wheels0_): MH_trainer(mhts_), gaussian_likelihood(mhts_.par.sigma, mhts_.sp_min.cl_im), t_wheels(t_wheels0_) {}
      ~basic_MH_trainer() {}

      double t_wheels; // current drift fraction
};
int find_worst_record(record ** r_, int ncap_);
int find_best_record(record ** r_, int ncap_);

#endif
