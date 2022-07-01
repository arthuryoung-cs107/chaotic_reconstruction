#ifndef MH_LEARNING_HH
#define MH_LEARNING_HH

#include "MH_auxiliary.hh"
#include "MH_tools.hh"

struct record: public record_struct
{
  record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_);
  ~record() {}

  const int rid;

  int * const ichunk;

  double  * dub_compare_bad,
          * const dchunk,
          * const u;

  virtual void write_ints(FILE * file_) {}
  virtual void write_dubs(FILE * file_) {}

  inline int isworse(record * rec_) {return (*dub_compare_bad)>(*(rec_->dub_compare_bad));}
  inline int isbetter(record * rec_) {return (*dub_compare_bad)<(*(rec_->dub_compare_bad));}
  inline void draw_ranuni(MH_rng * ran_, double * umin_, double * umax_)
  {for (int i = 0; i < ulen; i++) u[i] = umin_[i]+(umax_[i]-umin_[i])*ran_->rand_uni();}
  inline int take_record_chunks(record *rec_)
  {
    if (rid!=rec_->rid)
    {
      memcpy(ichunk,rec_->ichunk,ichunk_len*sizeof(int));
      memcpy(dchunk,rec_->dchunk,dchunk_len*sizeof(double));
      memcpy(u,rec_->u,ulen*sizeof(double));
      return 1;
    }
    else return 0;
  }
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

      thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_);

      ~thread_worker() {delete psim;}

    protected:
      const int thread_id;

      double  * const u = &(Kn),
              * const psim;

      void reset_sim(double *utest_, double t0_, double ctheta0_, double comega0_, double *p0_);

      inline double * advance_sim(int f_local_,double *t_history_)
      {
        advance((ts[f_local_]-ts[f_local_-1])/t_phys,d_ang[f_local_-1], comega_s[f_local_],dt_sim);
        t_history_[2]=t_history_[1]; t_history_[1]=t_history_[0]; t_history_[0]=ts[f_local_];
        return xs+(f_local_*ndof);
      }
      inline double compute_residual(double xs_, double ys_, double xr_, double yr_)
      {
        double  x_now=(xs_-cx)*cl_im+cx_im, y_now=(ys_-cy)*cl_im+cy_im,
                xerr=x_now-xr_, yerr=y_now-yr_;
        return xerr*xerr+yerr*yerr;
      }
      inline double compute_residual(double xs_, double ys_, double &x_now_, double &y_now_, double xr_, double yr_)
      {
        x_now_=(xs_-cx)*cl_im+cx_im; y_now_=(ys_-cy)*cl_im+cy_im;
        double xerr=x_now_-xr_, yerr=y_now_-yr_;
        return xerr*xerr+yerr*yerr;
      }
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

    swirl_param sp_min, // lower boundary of parameter space U
                sp_max; // upper boundary of parameter space U

    const int nt, // number of worker threads
              ichunk_width,
              dchunk_width;

    int leader_count, // 0: current number of leaders
        gen_count, // 1: how many pools have been drawn
        nsuccess, // 2: number of parameters better than current worst leader in recent trial
        ncandidates, // 3: number of candidates to compare to current leaders
        bleader_rid, // 4: index of current best leader
        wleader_rid, // 5: index of current worst leader
        nreplace, // 6: number of leader replacements following evaluation of recent trial
        ndup, // 7: total number of leader duplication and perturbations from recent resampling
        ndup_unique, // 8: total number of unique duplications from recent resampling
        nredraw, // 9: total number of trial particles drawn from proposal distribution in recent resampling
        * const MHT_it_ints = &leader_count,
        ** const ichunk; // space for records to store integer parameters

    double  rho2, // 0: current expected residual
            br2, // 1: current best leader residual
            wr2, // 2: current worst leader residual
            * const MHT_it_dubs = &rho2,
            * const umin = &(sp_min.Kn), // lower bound on u parameters
            * const umax = &(sp_max.Kn), // upper bound on u parameters
            * const ts, // wall time of observed data
            * const xs, // 2D observed position data
            * const d_ang, // observed angular position of dish
            * const comega_s, // observed average rotational speed of dish
            ** const uchunk, // space for parameters
            ** const dchunk; // space for records to store double parameters

    wall_list &wl; // a reference to the list of walls for the swirling simulation.
    proximity_grid ** const pg; // array of proximity grids.
    MH_rng ** rng; // random number generators

    inline void redraw_u_uni(record * rec_pool_, MH_rng * rng_t_)
    {rec_pool_->draw_ranuni(rng_t_,umin,umax);}
    inline void initialize_MHT_run()
    {
      for (int i = 0; i < MHT_it_ilen; i++) MHT_it_ints[i]=0;
      for (int i = 0; i < MHT_it_dlen; i++) MHT_it_dubs[i]=0.0;
    }
    inline void write_MHT_it_ints(FILE * file_) {fwrite(MHT_it_ints, sizeof(int), MHT_it_ilen, file_);}
    inline void write_MHT_it_dubs(FILE * file_) {fwrite(MHT_it_dubs, sizeof(double), MHT_it_dlen, file_);}
    inline double max(double a_,double b_ ) {return (a_>b_)?a_:b_;}
    inline double min(double a_,double b_ ) {return (a_<b_)?a_:b_;}
};

const int basic_rec_ilen=7;
const int basic_rec_dlen=2;
struct basic_record: public record
{
  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record(rs_, rid_, ichunk_, dchunk_, u_) {init_basic_record();}

  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double * umin_, double * umax_): basic_record(rs_, rid_, ichunk_, dchunk_, u_) {draw_ranuni(ran_,umin_,umax_);}

  ~basic_record() {}

  bool success;

  int gen, // 0: generation in which this particle was generated
      Class, // 1: leader Class this particle belongs to, if any
      dup_count, // 2: number of times this particle has been duplicated
      parent_count, // 3: number of particle ancestors
      parent_rid, // 4: index position of parent particle
      parent_gen, // 5: generation this particle comes from
      parent_Class, // 6: leader Class this particle comes from
      * const basic_rec_ints = &gen;

  double  r2,
          w,
          * const basic_rec_dubs = &r2;

  inline double get_r2() {return *dub_compare_bad;}

  inline void init_basic_record(int gen_=0,int Class_=-1,int dup_count_=0,int parent_count_=0,int parent_rid_=-1,int parent_gen_=-1,int parent_Class_=-1,double r2_=0.0,double w_=0.0)
  {
    gen=gen_; Class=Class_; dup_count=dup_count_;
    parent_count=parent_count_; parent_rid=parent_rid_;
    parent_gen=parent_gen_; parent_Class=parent_Class_;
    r2=r2_; w=w_;
  }

  inline int take_basic_record(basic_record * rec_)
  {
    if (take_record_chunks(rec_)) // successful replacement
    {
      memcpy(basic_rec_ints, rec_->basic_rec_ints, basic_rec_ilen*sizeof(int));
      memcpy(basic_rec_dubs, rec_->basic_rec_dubs, basic_rec_dlen*sizeof(double));
      return 1;
    }
    else return 0;
  }

  inline int basic_rec_ilen_full() {return basic_rec_ilen;}
  inline int basic_rec_dlen_full() {return basic_rec_dlen;}

  inline void write_basic_rec_ints(FILE * file_)
  {fwrite(basic_rec_ints, sizeof(int), basic_rec_ilen, file_);}

  inline void write_basic_rec_dubs(FILE * file_)
  {fwrite(basic_rec_dubs, sizeof(double), basic_rec_dlen, file_);}

};

const int basic_tw_ilen=4;
const int basic_tw_dlen=4;
class basic_thread_worker: public thread_worker, public event_detector
{
    public:

      basic_thread_worker(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_): thread_worker(sp_, pg_, wl_, tws_, thread_id_), event_detector(nbeads, Frames, 2, alpha_tol_),
      int_wkspc(new int[npool]), dub_wkspc(new double[ulen]) {}

      ~basic_thread_worker() {delete [] int_wkspc; delete [] dub_wkspc;}

      int * const int_wkspc;
      double * const dub_wkspc;

      inline void clear_basic_tw_training_data() {memset(int_wkspc,0,npool*sizeof(int));}

    protected:

      int nf_obs,
          nf_stable,
          nf_regime,
          nf_unstable,
          * const basic_tw_ints = &nf_obs;

      double  net_r2,
              net_r2_stable,
              net_r2_regime,
              net_r2_unstable,
              * const basic_tw_dubs = &net_r2;

      inline void update_integral_history(double INT_now_, int ibead_)
      {
        INTr2_comp_history[ibead_][2]=INTr2_comp_history[ibead_][1];
        INTr2_comp_history[ibead_][1]=INTr2_comp_history[ibead_][0];
        INTr2_comp_history[ibead_][0]+=INT_now_;
      }
};

class basic_MH_trainer: public MH_trainer, public gaussian_likelihood
{
    public:

      basic_MH_trainer(MH_train_struct &mhts_, int ichunk_width_, int dchunk_width_, double t_wheels0_=-1.0);

      ~basic_MH_trainer()
      {
        delete [] ndup_leaders; delete [] irepl_leaders; delete [] isuccess_pool;
        delete [] w_leaders;
        delete [] u_mean; delete [] u_var;
        delete [] u_wmean; delete [] u_wvar;
      }

    protected:

      bool apply_training_wheels;

      const double t_wheels0; // initial drift fraction

      double t_wheels; // current drift fraction

      int * const ndup_leaders,
          * const irepl_leaders,
          * const isuccess_pool;

      double  * const w_leaders,
              * const u_mean,
              * const u_var,
              * const u_wmean,
              * const u_wvar;

      virtual void duplicate_u(basic_record *rec_child_, basic_record *rec_parent_, MH_rng *rng_t_);

      virtual void redraw_u(basic_record * rec_pool_, MH_rng * rng_t_)
      {redraw_u_uni(rec_pool_, rng_t_); rec_pool_->init_basic_record(gen_count);}
      inline void clear_basic_trainer_training_data() {memset(isuccess_pool,0,npool*sizeof(int));nsuccess=0;}
      inline void initialize_basic_trainer_run() {initialize_MHT_run(); t_wheels=t_wheels0;}
      inline int basic_train_it_ilen_full() {return MHT_it_ilen;}
      inline int basic_train_it_dlen_full() {return MHT_it_dlen;}
      inline void write_basic_train_it_ints(FILE * file_) {write_MHT_it_ints(file_);}
      inline void write_basic_train_it_dubs(FILE * file_) {write_MHT_it_dubs(file_);}
      inline void write_ustats(FILE *file_)
      {
        fwrite(u_mean,sizeof(double),ulen,file_);
        fwrite(u_wmean,sizeof(double),ulen,file_);
        fwrite(u_var,sizeof(double),ulen,file_);
        fwrite(u_wvar,sizeof(double),ulen,file_);
      }
};

#endif
