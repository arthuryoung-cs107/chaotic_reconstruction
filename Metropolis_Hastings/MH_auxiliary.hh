#ifndef MH_AUXILIARY_HH
#define MH_AUXILIARY_HH

#include "swirl.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

const int full_ulen=14;
const int special_u_count=5;
const double special_u[][full_ulen] = {
/* 0  1    2      3     4     5     6   7    8   9   10    11    12    13
 {rad,mass,Kn    ,gb   ,gf   ,gw   ,mb ,mf  ,mw ,ds ,cx_im,cy_im,cl_im,wallsca} */
  0.5,1.0, 1000.0,40.0 ,40.0 ,40.0 ,0.5,0.25,0.5,1.8,203.0,178.0,27.6 ,1.0, // 0: default from Rycroft
  0.5,1.0, 500.0 ,5.0  ,5.0  ,5.0  ,0.1,0.1 ,0.1,1.8,203.0,178.0,27.6 ,1.0, // 1: min from Rycroft
  0.5,1.0, 5000.0,120.0,120.0,120.0,1.0,1.0 ,1.0,1.8,203.0,178.0,27.6 ,1.0, // 2: max from Rycroft
  0.5,1.0,1.193724226206541e3,0.047740229534684e3,0.050033742603846e3,0.025898751172936e3,0.000773880324000e3,0.000131158234738e3,0.000439930156805e3,1.8,203.0,178.0,27.6,1.0, // 3: perturbed from stat test
  0.5,1.0,968.783014762039e+000,41.2628475974082e+000,37.6024843085871e+000,40.1025717087601e+000,786.658543983309e-003,248.062481196825e-003,508.720872282600e-003,1.8,203.0,178.0,27.6,1.0 // 4: solution from nbeads=3, par_id=0, relay_id=5 relay
};

int set_special_u(int id_, double *vec_);
int set_special_u(const char *id_, double*vec_);

const int swirl_system_struct_const_ilen=3;
struct swirl_system_struct
{
  swirl_system_struct(int ulen_, int nbeads_, int Frames_): ulen(ulen_), nbeads(nbeads_), Frames(Frames_) {}
  swirl_system_struct(swirl_system_struct &sss_): swirl_system_struct(sss_.ulen, sss_.nbeads, sss_.Frames) {}
  ~swirl_system_struct() {}

  const int ulen, // length of parameter vector
            nbeads, // number of beads
            Frames; // number of observation frames
};

struct MH_io
{
  MH_io(char * proc_loc_, char * test_dir_, char * data_name_, int id_, bool noise_data_=false, double noise_sigma_=0.0);
  ~MH_io();

  const size_t fullbuf_len;
  size_t  obuf_len,
          ibuf_len;

  const int id;

  int Frames,
      nbeads;

  char  * const proc_loc,
        * const test_dir,
        * const data_name,
        * const obuf,
        * const ibuf;

  void load_reference(double *ts_, double *xs_, double *d_ang_, double * comega_s_, double t_phys_);

  private:
    const bool noise_data;
    const double noise_sigma;

    void read_fisml(char * ibuf_);
};

struct record_struct: public virtual swirl_system_struct
{
  record_struct(int ulen_, int nbeads_, int Frames_, int ichunk_len_, int dchunk_len_): swirl_system_struct(ulen_, nbeads_, Frames_), ichunk_len(ichunk_len_), dchunk_len(dchunk_len_) {}

  record_struct(record_struct &rs_): record_struct(rs_.ulen, rs_.nbeads, rs_.Frames, rs_.ichunk_len, rs_.dchunk_len) {}

  ~record_struct() {}

  const int ichunk_len,
            dchunk_len;
};

struct thread_worker_struct: public virtual swirl_system_struct
{
  thread_worker_struct(int ulen_, int nbeads_, int Frames_, int nlead_, int npool_, double dt_sim_, double t_phys_, double *ts_, double *xs_, double *d_ang_, double *comega_s_): swirl_system_struct(ulen_, nbeads_, Frames_),
  nlead(nlead_), npool(npool_), dt_sim(dt_sim_), t_phys(t_phys_), ts(ts_), xs(xs_), d_ang(d_ang_), comega_s(comega_s_) {}

  // thread_worker_struct(thread_worker_struct &tws_): swirl_system_struct(tws_.ulen, tws_.nbeads, tws_.Frames),
  // nlead(tws_.nlead), npool(tws_.npool), dt_sim(tws_.dt_sim), t_phys(tws_.t_phys), ts(tws_.ts), xs(tws_.xs), d_ang(tws_.d_ang), comega_s(tws_.comega_s) {}

  thread_worker_struct(thread_worker_struct &tws_): thread_worker_struct(tws_.ulen, tws_.nbeads, tws_.Frames, tws_.nlead, tws_.npool, tws_.dt_sim, tws_.t_phys, tws_.ts, tws_.xs, tws_.d_ang, tws_.comega_s) {}

  ~thread_worker_struct() {}

  const int nlead,
            npool;

  const double  dt_sim,
                t_phys;

  double  * const ts,
          * const xs,
          * const d_ang,
          * const comega_s;
};

struct MH_params: public virtual swirl_system_struct
{
  MH_params(MH_io *io_, int ulen_, int nlead_, int npool_, double dt_sim_, double t_phys_, double sigma_): swirl_system_struct(ulen_, io_->nbeads, io_->Frames),
  io(io_),
  nlead(nlead_), npool(npool_), dt_sim(dt_sim_), t_phys(t_phys_), sigma(sigma_) {}

  MH_params(MH_params &par_): MH_params(par_.io,par_.ulen,par_.nlead,par_.npool,par_.dt_sim,par_.t_phys,par_.sigma) {}

  ~MH_params() {}

  const int nlead, // number of leaders we store for review
            npool; // number of parameter sets we evaluate at a time

  const double  dt_sim, // The maximum simulation timestep to use
                t_phys, // simulation time scale
                sigma; // expected standard deviation of noise on reference data

  MH_io * const io;

  protected:

    inline void write_MH_params(FILE * file_)
    {
      int hlen=2, ilen=5, dlen=3,
          header[] = {hlen, ilen, dlen},
          ints[] = {ulen, nbeads, Frames, nlead, npool};
      double dubs[] = {dt_sim, t_phys, sigma};
      fwrite(header, sizeof(int), hlen+1, file_);
      fwrite(ints, sizeof(int), ilen, file_);
      fwrite(dubs, sizeof(double), dlen, file_);
    }

    #ifdef _OPENMP
      inline int thread_num() {return omp_get_thread_num();}
      inline int get_nt() {return omp_get_max_threads();}
    #else
      inline int thread_num() {return 0;}
      inline int get_nt() {return 1;}
    #endif
};

struct MH_train_struct
{
  MH_train_struct(MH_params *par_, swirl_param *sp_min_, swirl_param *sp_max_, wall_list *wl_): par(par_), sp_min(sp_min_), sp_max(sp_max_), wl(wl_) {}
  MH_train_struct(MH_train_struct &mhts_): MH_train_struct(mhts_.par, mhts_.sp_min, mhts_.sp_max, mhts_.wl) {}
  ~MH_train_struct() {}

  MH_params * const par;
  swirl_param * const sp_min,
              * const sp_max;
  wall_list * const wl;

  inline int get_par_nbeads() {return par->nbeads;}
};

#endif
