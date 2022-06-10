#ifndef MH_LEARNING_HH
#define MH_LEARNING_HH

#include "MH_auxiliary.hh"
#include "MH_tools.hh"

int find_worst_record(record ** r_, int ncap_);
int find_best_record(record ** r_, int ncap_);

struct record: public record_struct
{
  record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record_struct(rec_), rid(rid_), ichunk(ichunk_), dchunk(dchunk_), u(u_) {}
  ~record() {}

  const int rid;

  int * const ichunk;

  double  * const dchunk,
          * const u;

  inline void draw_ranuni(MH_rng * ran_, double * umin_, double * umax_)
  {for (int i = 0; i < ulen; i++) u[i] = umin_[i]+(umax_[i]-umin_[i])*ran_->rand_uni();}

  virtual int isworse(record * r_) = 0;
  virtual int isbetter(record * r_) = 0;

};

class thread_worker: public swirl, public thread_worker_struct
{
    public:

      thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_): swirl(sp_, pg_, wl_, tws.nbeads), thread_worker_struct(tws_),
      thread_id(thread_id_),
      u(&Kn), p(new double[2*nbeads*Frames]) {}
      ~thread_worker() {delete p;}

    protected:

      const int thread_id;

      double  * const u,
              * const p;
};


const int MH_train_ilen=9;
const int MH_train_dlen=3;
class MH_trainer : public MH_params
{
  public:

    MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_, int ichunk_width_, int dchunk_width_);
    MH_trainer(MH_train_struct &mhts_, int ichunk_width_, int dchunk_width_)
    : MH_trainer(*(mhts_.par), *(mhts_.sp_min), *(mhts_.sp_max), *(mhts_.wl), ichunk_width_, dchunk_width_) {}
    ~MH_trainer();

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

    protected:

      const int nt, // number of worker threads
                ichunk_width,
                dchunk_width;

      int ** const ichunk; // space for records to store integer parameters

      double * const  ts, // wall time of observed data
             * const  xs, // 2D observed position data
             * const  d_ang, // observed angular position of dish
             * const  comega_s, // observed average rotational speed of dish
             ** const uchunk, // space for parameters
             ** const dchunk; // space for records to store double parameters

      swirl_param sp_min, // lower boundary of parameter space U
                  sp_max; // upper boundary of parameter space U
      wall_list &wl; // a reference to the list of walls for the swirling simulation.
      proximity_grid ** const pg; // array of proximity grids.
      MH_rng ** rng; // random number generators
};

const int basic_rec_ilen=5;
const int basic_rec_dlen=2;
struct basic_record: public record
{
  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record(rs_, rid_, ichunk_, dchunk_, u_), basic_rec_ints(&gen), basic_rec_dubs(&r2) {}
  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double * umin_, double * umax_): basic_record(rs_, rid_, ichunk_, dchunk_, u_) {draw_ranuni(ran_,umin_,umax_);}
  ~basic_record();

  bool success;

  int gen, // generation in which this particle was generated
      parent_gen, // generation this particle comes from
      parent_count, // number of particle ancestors
      parent_rid, // index position of parent particle
      dup_count; // number of times this particle has been duplicated

  int * const basic_rec_ints;

  double  r2,
          w;

  double * const basic_rec_dubs;
};

class basic_thread_worker: public thread_worker, public event_detector
{
    public:

      basic_thread_worker(thread_work_struct &tws_, int thread_id_, double alpha_tol_): thread_worker(tws_, thread_id_), event_detector(nbeads, Frames, alpha_tol_) {}
      ~basic_thread_worker() {}
};

class basic_MH_trainer: public MH_trainer, public gaussian_likelihood
{
    public:

      basic_MH_trainer(MH_train_struct &mhts_, int ichunk_width_, int dchunk_width_, double t_wheels0_=-1.0): MH_trainer(mhts_, ichunk_width_, dchunk_width_), gaussian_likelihood(mhts_.par.sigma, mhts_.sp_min.cl_im),
      t_wheels(t_wheels0_), apply_training_wheels(t_wheels>0.0)
      umin(&(sp_min.Kn)), umax(&(sp_max.Kn)) {}
      ~basic_MH_trainer() {}

    protected:

      bool apply_training_wheels;

      double t_wheels; // current drift fraction

      double  * const umin,
              * const umax;
};

#endif
