#include "MH_learning.hh"

extern "C"
{
  #include "AYaux.h"
}

thread_worker::thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, thread_worker_struct &tws, int thread_id_): swirl(sp_, pg_, wl_, tws.nbeads), thread_worker_struct(tws), thread_id(thread_id_),
ichunk(new int[tws.ichunk_len]), dchunk(new double[tws.dchunk_len])
{u = &Kn;}

thread_worker::~thread_worker()
{delete [] ichunk; delete [] dchunk;}

MH_trainer::MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_): MH_params(par_),
sp_min(sp_min_), sp_max(sp_max_), wl(wl_),
ts(new double[par_.Frames]), xs(new double[2*par_.nbeads*par_.Frames]), d_ang(new double[par_.Frames]), comega_s(new double[par_.Frames]),
nt(get_nt()),
pg(new proximity_grid*[nt]), rng(new AYrng*[nt]),
u_chunk(AYdmatrix(par_.nlead+par_.npool, par_.param_len))
{
  io->load_reference(ts, xs, d_ang, comega_s, t_phys);

  // Set up each runner's personal data
#pragma omp parallel
  {
    int t=thread_num();
    pg[t]=new proximity_grid();
    rng[t]=new AYrng();
    rng[t]->rng_init_gsl(t+1);
  }
}

MH_trainer::~MH_trainer()
{
  delete [] ts;
  delete [] xs;
  delete [] d_ang;
  delete [] comega_s;
  free_AYdmatrix(u_chunk);

  for (int i = 0; i < nt; i++)
  {
    delete pg[i];
    delete rng[i];
  }
  delete [] pg;
  delete [] rng;
}

basic_thread_worker::basic_thread_worker(thread_work_strct &tws_, int thread_id_, double alpha_tol_): thread_worker(tws_, thread_id_), event_detector(alpha_tol_), r2_bead(new double*[Frames_]), INTr2_bead()
{

}

int find_worst_record(record ** r_, int ncap_)
{
  int worst_index = 0;
  for (int i = 1; i < ncap_; i++)
    if (r_[worst_index]->isbetter(r_[i]))
      worst_index = i;

  return worst_index;
}

int find_best_record(record ** r_, int ncap_)
{
  int best_index = 0;
  for (int i = 1; i < ncap_; i++)
    if (r_[best_index]->isworse(r_[i]))
      best_index = i;
  return best_index;
}
