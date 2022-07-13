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
  inline void print_record(const char indent_[]="  ")
  {printf("%s(record): ichunk_len=%d, dchunk_len=%d\n", indent_, ichunk_len,dchunk_len);}
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
        return xs+(f_local_*2*nbeads);
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
      inline double compute_diff_residual(double xdiff_, double ydiff_, double xr_, double yr_)
      {
        double  x_now=xdiff_*cl_im+cx_im, y_now=ydiff_*cl_im+cy_im,
                xerr=x_now-xr_, yerr=y_now-yr_;
        return xerr*xerr+yerr*yerr;
      }
      inline double compute_diff_residual(double xdiff_, double ydiff_, double &xnow_, double &ynow_, double xr_, double yr_)
      {
        xnow_=xdiff_*cl_im+cx_im; ynow_=ydiff_*cl_im+cy_im;
        double xerr=xnow_-xr_, yerr=ynow_-yr_;
        return xerr*xerr+yerr*yerr;
      }
};


const int MHT_it_ilen=11;
const int MHT_it_dlen=3;
class MH_trainer : public MH_params
{
  public:

    MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_, int ichunk_width_, int dchunk_width_);
    MH_trainer(MH_train_struct &mhts_, int ichunk_width_, int dchunk_width_)
    : MH_trainer(*(mhts_.par), *(mhts_.sp_min), *(mhts_.sp_max), *(mhts_.wl), ichunk_width_, dchunk_width_) {}
    ~MH_trainer();

    virtual void run(bool verbose_) = 0;

  protected:

    swirl_param sp_min, // lower boundary of parameter space U
                sp_max; // upper boundary of parameter space U

    const int nt, // number of worker threads
              ichunk_width,
              dchunk_width;

    int leader_count, // 0: current number of leaders
        gen_count,    // 1: how many pools have been drawn
        nsuccess,     // 2: number of parameters better than current worst leader in recent trial
        ncandidates,  // 3: number of candidates to compare to current leaders
        bleader_rid,  // 4: index of current best leader
        wleader_rid,  // 5: index of current worst leader
        nreplace,     // 6: number of leader replacements following evaluation of recent trial
        ndup,         // 7: total number of leader duplication and perturbations from recent resampling
        ndup_unique,  // 8: total number of unique duplications from recent resampling
        nredraw,      // 9: total number of trial particles drawn from proposal distribution in recent resampling
        nreload,      // 10: total number of leader particles reloaded into pool
        * const MHT_it_ints = &leader_count,
        ** const ichunk; // space for records to store integer parameters

    double  rho2, // 0: current expected residual
            br2,  // 1: current best leader residual
            wr2,  // 2: current worst leader residual
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

    // debugging

    inline void print_MHT(const char indent_[])
    {
      printf("%s(MH_trainer) leader_count: %d\n", indent_,leader_count);
      printf("%s(MH_trainer) gen_count: %d\n", indent_,gen_count);
      printf("%s(MH_trainer) nsuccess: %d\n", indent_,nsuccess);
      printf("%s(MH_trainer) ncandidates: %d\n", indent_,ncandidates);
      printf("%s(MH_trainer) bleader_rid: %d\n", indent_,bleader_rid);
      printf("%s(MH_trainer) wleader_rid: %d\n", indent_,wleader_rid);
      printf("%s(MH_trainer) nreplace: %d\n", indent_,nreplace);
      printf("%s(MH_trainer) ndup: %d\n", indent_,ndup);
      printf("%s(MH_trainer) ndup_unique: %d\n", indent_,ndup_unique);
      printf("%s(MH_trainer) nredraw: %d\n", indent_,nredraw);

      printf("%s(MH_trainer) rho2: %e\n",indent_,rho2);
      printf("%s(MH_trainer) br2: %e\n",indent_,br2);
      printf("%s(MH_trainer) wr2: %e\n",indent_,wr2);
    }

    // run
    inline void initialize_MHT_run()
    {
      for (int i = 0; i < MHT_it_ilen; i++) MHT_it_ints[i]=0;
      for (int i = 0; i < MHT_it_dlen; i++) MHT_it_dubs[i]=0.0;
    }

    // sampling
    inline void redraw_u_uni(record * rec_pool_, MH_rng * rng_t_)
    {rec_pool_->draw_ranuni(rng_t_,umin,umax);}

    // io
    virtual void stage_diagnostics() = 0;
    virtual void close_diagnostics() = 0;
    inline void write_MHT_it_ints(FILE * file_) {fwrite(MHT_it_ints, sizeof(int), MHT_it_ilen, file_);}
    inline void write_MHT_it_dubs(FILE * file_) {fwrite(MHT_it_dubs, sizeof(double), MHT_it_dlen, file_);}

    // aux
    inline double max(double a_,double b_ ) {return (a_>b_)?a_:b_;}
    inline double min(double a_,double b_ ) {return (a_<b_)?a_:b_;}
};

const int basic_rec_ilen=7;
const int basic_rec_dlen=1;
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

  double  w,
          * const basic_rec_dubs = &w;

  // sampling

  inline void init_basic_record(int gen_=0,int Class_=-1,int dup_count_=0,int parent_count_=0,int parent_rid_=-1,int parent_gen_=-1,int parent_Class_=-1,double w_=0.0)
  {
    gen=gen_; Class=Class_; dup_count=dup_count_;
    parent_count=parent_count_; parent_rid=parent_rid_;
    parent_gen=parent_gen_; parent_Class=parent_Class_;
    w=w_;
  }

  // debugging

  inline void print_basic_record(const char indent_[]=" ", bool print_all_=true)
  {
    printf("%s(basic_record) int data:\n", indent_);
    printf("%s gen: %d\n",indent_,gen);
    printf("%s Class: %d\n",indent_,Class);
    printf("%s dup_count: %d\n",indent_,dup_count);
    printf("%s parent_count: %d\n",indent_,parent_count);
    printf("%s parent_rid: %d\n",indent_,parent_rid);
    printf("%s parent_gen: %d\n",indent_,parent_gen);
    printf("%s parent_Class: %d\n",indent_,parent_Class);

    printf("%s(basic_record) double data:\n", indent_);
    printf("%s w: %e\n",indent_,w);
    printf("%s r2: %e\n",indent_,get_r2());
    if (print_all_) print_record();
  }

  // training

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

  inline double get_r2() {return *dub_compare_bad;}

  // io
  inline int basic_rec_ilen_full() {return basic_rec_ilen;}
  inline int basic_rec_dlen_full() {return basic_rec_dlen;}

  inline void write_basic_rec_ints(FILE * file_)
    {fwrite(basic_rec_ints, sizeof(int), basic_rec_ilen, file_);}

  inline void write_basic_rec_dubs(FILE * file_)
    {fwrite(basic_rec_dubs, sizeof(double), basic_rec_dlen, file_);}
};

class basic_thread_worker: public thread_worker, public event_detector
{
    public:

      basic_thread_worker(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_, int iwkspc_len_, int dwkspc_len_): thread_worker(sp_, pg_, wl_, tws_, thread_id_), event_detector(nbeads, Frames, 2, alpha_tol_),
      iwkspc_len(iwkspc_len_), dwkspc_len(dwkspc_len_),
      int_wkspc(new int[iwkspc_len]), dub_wkspc(new double[dwkspc_len]) {}

      ~basic_thread_worker() {delete [] int_wkspc; delete [] dub_wkspc;}

      const int iwkspc_len,
                dwkspc_len;

      int * const int_wkspc;

      double  netr2,
              netr2_stable,
              netr2_unstable,
              *r2_objective,
              * const dub_wkspc;

      inline void set_net_objective()
        {r2_objective=&netr2; rho2_objective=compute_netrho2();}
      inline void set_stable_objective()
        {r2_objective=&netr2_stable; rho2_objective=compute_netrho2stable();}
      inline void set_unstable_objective()
        {r2_objective=&netr2_unstable; rho2_objective=compute_netrho2unstable();}

      inline void clear_basic_tw_training_data()
      {
        memset(int_wkspc,0,iwkspc_len*sizeof(int));
        for (int i = 0; i < dwkspc_len; i++) dub_wkspc[i]=0.0;
      }

      inline bool check_success(double r2_threshold_) {return (*r2_objective)<r2_threshold_;}

    protected:

      // debugging
      inline void print_basic_tw(const char indent_[], double sigma_scaled_)
        {print_event_block(sigma_scaled_,indent_);}
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

      double  t_wheels, // current drift fraction
              r2_scale;

      int * const ndup_leaders,
          * const irepl_leaders,
          * const isuccess_pool;

      double  * const w_leaders,
              * const u_mean,
              * const u_var,
              * const u_wmean,
              * const u_wvar;

      // debugging
      inline void print_basic_MHT(const char indent_[]="     ")
        {print_MHT(indent_); print_gaussian_likelihood(indent_);}

      // run
      inline void initialize_basic_MHT_run() {initialize_MHT_run(); t_wheels=t_wheels0;}

      // sampling
      virtual void redraw_u(basic_record * rec_pool_, MH_rng * rng_t_)
        {redraw_u_uni(rec_pool_, rng_t_); rec_pool_->init_basic_record(gen_count);}

      void duplicate_u_basic(basic_record *child_, basic_record *parent_, MH_rng *rng_t_, double fac_low_=0.0);

      // training
      inline void clear_basic_MHT_training_data() {memset(isuccess_pool,0,npool*sizeof(int));nsuccess=0;}

      // io
      inline int basic_MHT_it_ilen_full() {return MHT_it_ilen;}
      inline int basic_MHT_it_dlen_full() {return MHT_it_dlen;}
      inline void write_basic_MHT_it_ints(FILE * file_) {write_MHT_it_ints(file_);}
      inline void write_basic_MHT_it_dubs(FILE * file_) {write_MHT_it_dubs(file_);}
      inline void write_ustats(FILE *file_)
      {
        fwrite(u_mean,sizeof(double),ulen,file_);
        fwrite(u_wmean,sizeof(double),ulen,file_);
        fwrite(u_var,sizeof(double),ulen,file_);
        fwrite(u_wvar,sizeof(double),ulen,file_);
      }
};

#endif
