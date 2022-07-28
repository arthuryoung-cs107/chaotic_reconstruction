#ifndef MH_LEARNING_HH
#define MH_LEARNING_HH

#include "MH_auxiliary.hh"
#include "MH_tools.hh"

struct record: public record_struct
{
  record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record_struct(rs_),
  rid(rid_),ichunk(ichunk_),dchunk(dchunk_), u(u_) {}
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

  // debugging
  void print_record(const char indent_[]="  ");

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
    void print_MHT(const char indent_[]);

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
    inline int max(int a_,int b_ ) {return (a_>b_)?a_:b_;}
    inline int min(int a_,int b_ ) {return (a_<b_)?a_:b_;}
};

#endif
