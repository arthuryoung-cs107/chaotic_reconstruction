// #include <sys/stat.h>
// #include <sstream>
// #include <string>
// #include <sstream>
// #include <fstream>

#include "MH_learning.hh"

extern "C"
{
  #include "AYaux.h"
}

MH_io::MH_io(char * proc_loc_, char * test_dir_, char * data_name_, int id_, bool noise_data_, double noise_sigma_): proc_loc(proc_loc_), test_dir(test_dir_), data_name(data_name_),
id(id_), fullbuf_len(strlen(proc_loc)+strlen(test_dir)+strlen(noise_data)),
noise_data(noise_data_), noise_sigma(noise_sigma_),
obuf(new char[fullbuf_len+100]), ibuf(new char[fullbuf_len+1])
{
  sprintf(obuf, "%s%s", proc_loc, test_dir); obuf_end = strlen(obuf);
  sprintf(ibuf, "%s%s%s", proc_loc, test_dir, data_name); ibuf_end = strlen(ibuf);

  FILE *refh_file = fopen();
}

thread_worker::thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int thread_id_, int ichunk_len_, int dchunk_len_, int param_len_, int nbeads_, int Frames_, int nlead_, int npool_, double dt_sim_, double t_phys_, double *ts_, double *xs_, double *d_ang_, double *comega_s_): swirl(sp_, pg_, wl_, nbeads_),
thread_id(thread_id_), ichunk_len(ichunk_len_), dchunk_len(dchunk_len_),
param_len(param_len_), nbeads(nbeads_), Frames(Frames_), nlead(nlead_), npool(npool_),
dt_sim(dt_sim_), t_phys(t_phys_),
ts(ts_), xs(xs_), d_ang(d_ang_), comega_s(comega_s_),
ichunk(new int[ichunk_len]), dchunk(new double[dchunk_len])
{pvals = &Kn;}

thread_worker::~thread_worker()
{delete [] ichunk; delete [] dchunk;}
