#ifndef MH_LEARNING_HH
#define MH_LEARNING_HH

#include "MH_auxiliary.hh"
#include "MH_tools.hh"

struct record: public record_struct
{
  record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record_struct(rs_), rid(rid_), ichunk(ichunk_), dchunk(dchunk_), u(u_) {}
  ~record() {}

  const int rid;

  int * const ichunk;

  double  * const dchunk,
          * const u;

  inline void draw_ranuni(MH_rng * ran_, double * umin_, double * umax_)
  {for (int i = 0; i < ulen; i++) u[i] = umin_[i]+(umax_[i]-umin_[i])*ran_->rand_uni();}

  virtual int isworse(record * r_) {printf("(record) WARNING: using uninitialized record comparison\n"); return 0;}
  virtual int isbetter(record * r_) {printf("(record) WARNING: using uninitialized record comparison\n"); return 0;}

  virtual int ilen_full() = 0;
  virtual int dlen_full() = 0;

  virtual void write_ints(FILE * file_) = 0;
  virtual void write_dubs(FILE * file_) = 0;
  inline void write_chunks(FILE * file_)
  {
    fwrite(ichunk, sizeof(int), ichunk_len, file_);
    fwrite(dchunk, sizeof(double), dchunk_len, file_);
    fwrite(u, sizeof(double), ulen, file_);
  }
  inline void write_record_data(FILE * file_)
  {
    write_ints(file_);
    write_dubs(file_);
    write_chunks(file_);
  }
};

class thread_worker: public swirl, public thread_worker_struct
{
    public:

      thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_): swirl(sp_, pg_, wl_, tws_.nbeads), thread_worker_struct(tws_),
      thread_id(thread_id_),
      u(&Kn), psim(new double[2*nbeads*Frames]) {}
      ~thread_worker() {delete psim;}

    protected:

      void reset_sim(double *utest_, double t0_, double ctheta0_, double comega0_, double *p0_);

      const int thread_id;

      double  * const u,
              * const psim;
};


const int MHT_ilen=10;
const int MHT_dlen=3;
class MH_trainer : public MH_params
{
  public:

    MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_, int ichunk_width_, int dchunk_width_);
    MH_trainer(MH_train_struct &mhts_, int ichunk_width_, int dchunk_width_)
    : MH_trainer(*(mhts_.par), *(mhts_.sp_min), *(mhts_.sp_max), *(mhts_.wl), ichunk_width_, dchunk_width_) {}
    ~MH_trainer();

    virtual void run(bool verbose_) = 0;
    virtual void stage_diagnostics() = 0;
    virtual void close_diagnostics() = 0;

  protected:

    const int nt, // number of worker threads
              ichunk_width,
              dchunk_width;

    int leader_count, // current number of leaders
        gen_count, // how many pools have been drawn
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

    int * const MHT_ints,
        ** const ichunk; // space for records to store integer parameters

    double * const  ts, // wall time of observed data
           * const  xs, // 2D observed position data
           * const  d_ang, // observed angular position of dish
           * const  comega_s, // observed average rotational speed of dish
           * const  MHT_dubs,
           ** const uchunk, // space for parameters
           ** const dchunk; // space for records to store double parameters

    swirl_param sp_min, // lower boundary of parameter space U
                sp_max; // upper boundary of parameter space U
    wall_list &wl; // a reference to the list of walls for the swirling simulation.
    proximity_grid ** const pg; // array of proximity grids.
    MH_rng ** rng; // random number generators

    virtual void write_it_ints(FILE * file_) {fwrite(MHT_ints, sizeof(int), MHT_ilen, file_);}
    virtual void write_it_dubs(FILE * file_) {fwrite(MHT_dubs, sizeof(double), MHT_dlen, file_);}
    virtual void initialize_run()
    {
      for (int i = 0; i < MHT_ilen; i++) MHT_ints[i]=0;
      for (int i = 0; i < MHT_dlen; i++) MHT_dubs[i]=0.0;
    }
    virtual int it_ilen_full() {return basic_MH_trainer::it_ilen_full() + genetic_train_it_ilen;}
    virtual int it_dlen_full() {return basic_MH_trainer::it_dlen_full() + genetic_train_it_dlen;}
};

const int basic_rec_ilen=7;
const int basic_rec_dlen=2;
struct basic_record: public record
{
  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record(rs_, rid_, ichunk_, dchunk_, u_), basic_rec_ints(&gen), basic_rec_dubs(&r2) {init_record();}
  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double * umin_, double * umax_): basic_record(rs_, rid_, ichunk_, dchunk_, u_) {draw_ranuni(ran_,umin_,umax_);}
  ~basic_record() {}

  bool success;

  int gen, // 0: generation in which this particle was generated
      Class, // 1: leader Class this particle belongs to, if any
      dup_count // 2: number of times this particle has been duplicated
      parent_count, // 3: number of particle ancestors
      parent_gen, // 4: generation this particle comes from
      parent_Class, // 5: leader Class this particle comes from
      parent_rid; // 6: index position of parent particle

  double  r2,
          w;

  int * const basic_rec_ints;

  double * const basic_rec_dubs;

  inline void init_record()
  {
    gen=dup_count=parent_count=0; r2=w=0.0;
    parent_gen=Class=parent_rid=-1;
  }

  virtual int isworse(basic_record * r_) {return r2>r_->r2;}
  virtual int isbetter(basic_record * r_) {return r2<r_->r2;}

  virtual void write_ints(FILE * file_) {fwrite(basic_rec_ints, sizeof(int), basic_rec_ilen, file_);}
  virtual void write_dubs(FILE * file_) {fwrite(basic_rec_dubs, sizeof(double), basic_rec_dlen, file_);}
  virtual int ilen_full() {return basic_rec_ilen;}
  virtual int dlen_full() {return basic_rec_dlen;}
};

const int basic_tw_ilen=3;
const int basic_tw_dlen=3;
class basic_thread_worker: public thread_worker, public event_detector
{
    public:

      basic_thread_worker(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_): thread_worker(sp_, pg_, wl_, tws_, thread_id_), event_detector(nbeads, Frames, 2, alpha_tol_),
      basic_tw_ints(&nf_obs), basic_tw_dubs(&net_r2) {}
      ~basic_thread_worker() {}

    protected:

      int nf_obs,
          nf_stable,
          nf_unstable;

      double  net_r2,
              net_r2_stable,
              net_r2_unstable;

      int * const basic_tw_ints;

      double  * const basic_tw_dubs;
};

class basic_MH_trainer: public MH_trainer, public gaussian_likelihood
{
    public:

      basic_MH_trainer(MH_train_struct &mhts_, int ichunk_width_, int dchunk_width_, double t_wheels0_=-1.0): MH_trainer(mhts_, ichunk_width_, dchunk_width_), gaussian_likelihood(sigma, mhts_.sp_min->cl_im),
      apply_training_wheels(t_wheels0_>0.0), t_wheels0(t_wheels0_),
      umin(&(sp_min.Kn)), umax(&(sp_max.Kn)) {}
      ~basic_MH_trainer() {}

    protected:

      bool apply_training_wheels;

      const double t_wheels0; // initial drift fraction

      double t_wheels; // current drift fraction

      double  * const umin,
              * const umax;

      virtual int it_ilen_full() {return MH_trainer::it_ilen_full();}
      virtual int it_dlen_full() {return MH_trainer::it_dlen_full();}

      virtual void initialize_run()
      {
        MH_trainer::initialize_run();
        t_wheels=t_wheels0;
      }
};

int find_worst_record(record ** r_, int ncap_);
int find_best_record(record ** r_, int ncap_);

#endif
