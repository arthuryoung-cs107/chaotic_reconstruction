#ifndef MH_AUXILIARY_HH
#define MH_AUXILIARY_HH

#include "swirl.hh"

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

struct MH_params
{
  MH_params(MH_io *io_, int param_len_, int nlead_, int npool_, double dt_sim_, double t_phys_, double sigma_): io(io_),
  param_len(param_len_), nbeads(io_->nbeads), Frames(io_->Frames), nlead(nlead_), npool(npool_), dt_sim(dt_sim_), t_phys(t_phys_), sigma(sigma_) {}
  MH_params(MH_params &par_): MH_params(par_.io,par_.param_len,par_.nlead,par_.npool,par_.dt_sim,par_.t_phys,par_.sigma) {}

  ~MH_params() {}

  const int param_len, // length of parameter vector
            nbeads, // number of beads
            Frames, // number of observation frames
            nlead, // number of leaders we store for review
            npool; // number of parameter sets we evaluate at a time


  const double  dt_sim, // The maximum simulation timestep to use
                t_phys, // simulation time scale
                sigma; // expected standard deviation of noise on reference data

  MH_io * const io;
};

struct record_struct
{
  record_struct(int param_len_, int nbeads_, int Frames_): param_len(param_len_), nbeads(nbeads_), Frames(Frames_) {}
  record_struct(record_struct &rec_): record_struct(rec_.param_len, rec_.nbeads, rec_.Frames) {}
  ~record_struct();

  const int param_len,
            nbeads,
            Frames;
};

struct thread_worker_struct
{
  thread_worker_struct(int param_len_, int nbeads_, int Frames_, int nlead_, int npool_, int ichunk_len_, int dchunk_len_, double dt_sim_, double t_phys_, double *ts_, double *xs_, double *d_ang_, double *comega_s_): param_len(param_len_), nbeads(nbeads_), Frames(Frames_), nlead(nlead_), npool(npool_), ichunk_len(ichunk_len_), dchunk_len(dchunk_len_), dt_sim(dt_sim_), t_phys(t_phys_), ts(ts_), xs(xs_), d_ang(d_ang_), comega_s(comega_s_) {}
  thread_worker_struct(thread_worker_struct &tws_): thread_worker_struct(tws_.param_len, tws_.nbeads, tws_.Frames, tws_.nlead, tws_.npool, tws_.ichunk_len, tws_.dchunk_len, tws_.dt_sim, tws_.t_phys, tws_.ts, tws_.xs, tws_.d_ang, tws_.comega_s) {}
  thread_worker_struct(int * ipars_, double * dpars_,double * ts_, double * xs_, double * d_ang_, double * comega_s_): thread_worker_struct(ipars_[0], ipars_[1], ipars_[2], ipars_[3], ipars_[4], get_ichunk_len(), get_dchunk_len(), ts_, xs_, d_ang_, comega_s_) {}

  ~thread_worker_struct() {}

  const int param_len,
            nbeads,
            Frames,
            nlead,
            npool,
            ichunk_len,
            dchunk_len;

  const double  dt_sim,
                t_phys;

  double  * const ts,
          * const xs,
          * const d_ang,
          * const comega_s;

  protected:

    virtual int get_ichunk_len() = 0;
    virtual int get_dchunk_len() = 0;
};

struct MH_train_struct
{
  MH_train_struct(MH_params *par_, swirl_param *sp_min_, swirl_param *sp_max_, wall_list *wl_): par(par_), sp_min(sp_min_), sp_max(sp_max_), wl(wl_) {}
  MH_train_struct(MH_train_struct &mhts_): MH_train_struct(mhts_.par, mhts_.sp_min, mhts_.sp_min, mhts_.sp_max, mhts_.wl) {}
  ~MH_train_struct() {}

  MH_params * const par;
  swirl_param * const sp_min,
              * const sp_max;
  wall_list * const wl;

  inline int get_record_len() {return par->nlead+par->npool;}
  inline int get_par_nbeads() {return par->io->nbeads;}
  inline int get_io_obuf_len() {return par->io->obuf_len;}
  inline int get_par_nbeads() {return par->nbeads;}
};

#endif
