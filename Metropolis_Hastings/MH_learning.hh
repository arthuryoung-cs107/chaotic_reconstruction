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

  virtual int take_record(record * rec_) = 0;

  virtual int isworse(record * r_) {printf("(record::isworse) WARNING: using uninitialized record comparison\n"); return 0;}
  virtual int isbetter(record * r_) {printf("(record::isbetter) WARNING: using uninitialized record comparison\n"); return 0;}

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


const int MHT_it_ilen=10;
const int MHT_it_dlen=3;
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
            br2, // current best leader residual
            wr2; // current worst leader residual

    int * const MHT_it_ints,
        ** const ichunk; // space for records to store integer parameters

    double * const  ts, // wall time of observed data
           * const  xs, // 2D observed position data
           * const  d_ang, // observed angular position of dish
           * const  comega_s, // observed average rotational speed of dish
           * const  MHT_it_dubs,
           ** const uchunk, // space for parameters
           ** const dchunk; // space for records to store double parameters

    swirl_param sp_min, // lower boundary of parameter space U
                sp_max; // upper boundary of parameter space U
    wall_list &wl; // a reference to the list of walls for the swirling simulation.
    proximity_grid ** const pg; // array of proximity grids.
    MH_rng ** rng; // random number generators

    inline void initialize_MHT_run()
    {
      for (int i = 0; i < MHT_it_ilen; i++) MHT_it_ints[i]=0;
      for (int i = 0; i < MHT_dlen; i++) MHT_it_dubs[i]=0.0;
    }
    inline int find_worst_record(record ** r_, int ncap_)
    {
      int worst_index = 0;
      for (int i = 1; i < ncap_; i++)
        if (r_[worst_index]->isbetter(r_[i]))
          worst_index = i;

      return worst_index;
    }
    inline int find_best_record(record ** r_, int ncap_)
    {
      int best_index = 0;
      for (int i = 1; i < ncap_; i++)
        if (r_[best_index]->isworse(r_[i]))
          best_index = i;
      return best_index;
    }
    inline void pick_nworst_records(record ** rin_, record ** rout_, int n_, int ncap_)
    {
      for (int i = 0; i < n_; i++) rout_[i]=rin_[i];
      int i_best_worst = find_best_record(rout_,n_);
      for (int i = n_; i < ncap_; i++)
        if (rout_[i_best_worst]->isbetter(rin_[i]))
        {
          rout_[i_best_worst] = rin_[i];
          i_best_worst=find_best_record(rout_,n_);
        }
    }
    inline void pick_nbest_record(record ** rin_, record ** rout_, int n_, int ncap_)
    {
      for (int i = 0; i < n_; i++) rout_[i]=rin_[i];
      int i_worst_best = find_worst_record(rout_,n_);
      for (int i = n_; i < ncap_; i++)
        if (rout_[i_worst_best]->isworse(rin_[i]))
        {
          rout_[i_worst_best] = rin_[i];
          i_worst_best=find_worst_record(rout_,n_);
        }
    }
    inline int take_records(record ** rin_, record ** rout_, int * repl_list_, int ncap_)
    {
      int nrepl=0;
      for (int i = 0; i < ncap_; i++) if (rout_[i]->take_record(rin_[i])) repl_list[nrepl++]=i;
      return nrepl;
    }
    inline void write_MHT_it_ints(FILE * file_) {fwrite(MHT_it_ints, sizeof(int), MHT_it_ilen, file_);}
    inline void write_MHT_it_dubs(FILE * file_) {fwrite(MHT_it_dubs, sizeof(double), MHT_dlen, file_);}
};

const int basic_rec_ilen=7;
const int basic_rec_dlen=2;
struct basic_record: public record
{
  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record(rs_, rid_, ichunk_, dchunk_, u_), basic_rec_ints(&gen), basic_rec_dubs(&r2) {init_basic_record();}
  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double * umin_, double * umax_): basic_record(rs_, rid_, ichunk_, dchunk_, u_) {draw_ranuni(ran_,umin_,umax_);}
  ~basic_record() {}

  bool success;

  int gen, // 0: generation in which this particle was generated
      Class, // 1: leader Class this particle belongs to, if any
      dup_count, // 2: number of times this particle has been duplicated
      parent_count, // 3: number of particle ancestors
      parent_rid, // 4: index position of parent particle
      parent_gen, // 5: generation this particle comes from
      parent_Class; // 6: leader Class this particle comes from

  double  r2,
          w;

  int * const basic_rec_ints;

  double * const basic_rec_dubs;

  virtual int isworse(basic_record * r_) {return r2>r_->r2;}
  virtual int isbetter(basic_record * r_) {return r2<r_->r2;}

  inline void init_basic_record(int gen_=0,int Class_=-1,int dup_count_=0,int parent_count_=0,int parent_rid_=-1,int parent_gen_=-1,int parent_Class_=-1,double r2_=0.0,double w_=0.0)
  {
    gen=gen_; Class=Class_; dup_count=dup_count_;
    parent_count=parent_count_; parent_rid=parent_rid_
    parent_gen=parent_gen_; parent_Class=parent_Class_;
    r2=r2_; w=w_;
  }
  inline int take_basic_record(basic_record * rec_)
  {
    if (rid!=rec_->rid) // successful replacement
    {
      memcpy(basic_rec_ints, rec_->basic_rec_ints, basic_rec_ilen*sizeof(int));
      memcpy(basic_rec_dubs, rec_->basic_rec_dubs, basic_rec_dlen*sizeof(double));
      return 1;
    }
    else return 0;
  }
  inline int basic_rec_ilen_full() {return basic_rec_ilen;}
  inline int basic_rec_dlen_full() {return basic_rec_dlen;}
  inline void write_basic_rec_ints(FILE * file_) {fwrite(basic_rec_ints, sizeof(int), basic_rec_ilen, file_);}
  inline void write_basic_rec_dubs(FILE * file_) {fwrite(basic_rec_dubs, sizeof(double), basic_rec_dlen, file_);}
};

const int basic_tw_ilen=3;
const int basic_tw_dlen=3;
class basic_thread_worker: public thread_worker, public event_detector
{
    public:

      basic_thread_worker(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_): thread_worker(sp_, pg_, wl_, tws_, thread_id_), event_detector(nbeads, Frames, 2, alpha_tol_),
        basic_tw_ints(&nf_obs), basic_tw_dubs(&net_r2),
        int_wkspc(new int[nlead]) {}
      ~basic_thread_worker() {delete [] int_wkspc;}

      int * const int_wkspc;

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
      umin(&(sp_min.Kn)), umax(&(sp_max.Kn)),
      ndup_leaders(new int[nlead]), irepl_leaders(new int[nlead]),
      w_leaders(new double[nlead]) {}

      ~basic_MH_trainer() {delete [] ndup_leaders; delete [] irepl_leaders; delete [] w_leaders;}

    protected:

      bool apply_training_wheels;

      const double t_wheels0; // initial drift fraction

      double t_wheels; // current drift fraction

      int * const ndup_leaders,
          * const irepl_leaders;

      double  * const umin,
              * const umax,
              * const w_leaders;

      virtual double compute_weights(double r2_min_, basic_record ** recs_, int n_);
      virtual void respawn_pool(double w_sum_, basic_thread_worker **tws_, basic_record ** pool_, basic_record ** leaders_);

      virtual void duplicate_u(basic_record * rec_pool_, basic_record * rec_lead_, MH_rng * rng_t_);
      virtual void redraw_u(basic_record * rec_pool_, MH_rng * rng_t_)
      {rec_pool_->draw_ranuni(rng_t_,umin,umax); rec_pool_->init_basic_record(gen_count);}
      virtual int basic_train_it_ilen_full() {return MHT_it_ilen;}
      virtual int basic_train_it_dlen_full() {return MHT_it_dlen;}
      inline void initialize_basic_trainer_run() {MH_trainer::initialize_MHT_run(); t_wheels=t_wheels0;}
};

#endif
