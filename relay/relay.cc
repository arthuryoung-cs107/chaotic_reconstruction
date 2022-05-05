#include "particle_relay.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

extern "C"
{
  #include "AYaux.h"
}

relay::relay(referee &ref_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_,double t_phys_,reporter * rep_): referee(ref_), sp_min(sp_min_), sp_max(sp_max_), wl(wl_), t_phys(t_phys_), rep(rep_), n(rep_->P), Frames(rep_->Frames),
ts(new double[Frames]), xs(new double[2*n*Frames]), d_ang(new double[Frames]), comega_s(new double[Frames]),
#ifdef _OPENMP
nt(omp_get_max_threads()), // each thread is a runner
#else
nt(1), // only one runner
#endif
pg(new proximity_grid*[nt]), rng(new AYrng*[nt]), runners(new runner*[nt])
{
  alloc_records(nt, Frames, n);
  rep->load_relay(ts, xs, d_ang);

  // initialize the dish velocity data. Will save time on calculations later
  comega_s[0]=0.0;
  for (int i = 1; i < Frames; i++)
  {
    double comega=d_ang[i]-d_ang[i-1];
    if(comega>M_PI) comega-=2*M_PI; else if(comega<-M_PI) comega+=2*M_PI;
    comega_s[i] = t_phys*comega/(ts[i]-ts[i-1]);
  }

  param_acc_factor = new double[param_len];
  double *dmax=&sp_max.Kn, d_tn = (double)(omp_get_max_threads());
  for (int i = 0; i < param_len; i++) param_acc_factor[i] = ((par_acc_P*0.1)/(((double)nlead)*dmax[i]))*(d_tn);

  // Set up each runner's personal data
#pragma omp parallel
  {
    int t=thread_num();
    pg[t]=new proximity_grid();
    rng[t]=new AYrng();
    rng[t]->rng_init_gsl(t+1);
    runners[t]=new runner(sp_min, pg[t], wl, n, t, param_len, Frames, nlead, npool, t_phys, dt_sim, alpha_tol, ts, xs, d_ang, comega_s);
  }
}

relay::~relay()
{
  delete [] ts;
  delete [] xs;
  delete [] d_ang;
  delete [] comega_s;
  for (int i = 0; i < nt; i++)
  {
    delete pg[i];
    delete rng[i];
    delete runners[i];
  }
  delete [] pg;
  delete [] rng;
  delete [] runners;

  delete param_acc_factor;
}

void relay::stage_diagnostics(int gen_max_)
{

  int startup_header[] = {10, 6};

  int startup_int_params[] = {
  gen_max_,
  rep->nlead=nlead,
  rep->npool=npool,
  rep->param_len=param_len,
  rep->beads=n,
  rep->Frames,
  record_int_len,
  record_double_len,
  record_int_chunk_count,
  record_double_chunk_count
  };

  double startup_double_params[] = {
  rep->dt_sim = dt_sim,
  rep->noise_tol = noise_tol,
  rep->alpha_tol = alpha_tol,
  rep->rs_full_factor=rs_full_factor,
  rep->data_scale=data_scale,
  rep->t_phys = t_phys
  };

  rep->lead_dup_count = lead_dup_count;
  rep->global_event_frame_count = *global_event_frame_count;
  rep->event_frames = event_frames;
  rep->gen_int_vec = &(gen_count);
  rep->postevent_int_vec = &(event_observations);

  rep->sample_weights = sample_weights;
  rep->lead_par_w_mean = lead_par_w_mean;
  rep->lead_par_w_var = lead_par_w_var;
  rep->gen_double_vec = &(residual_best);
  rep->postevent_double_vec = &(tau);

  rep->leaders = leaders;
  rep->pool = pool;

  rep->staged_flag = true;
  rep->write_startup_diagnostics(startup_header, startup_int_params, startup_double_params);
}

void relay::init_relay()
{
  // initialize the pool of testing particles
  leader_count=gen_count=0;
  residual_worst = DBL_MAX;
#pragma omp parallel
  {
    AYrng *r = rng[thread_num()];
    double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
#pragma omp for
    for (int i = 0; i < npool; i++)
      pool[i]->init_pool(dmin, dmax, r);
  }
  for (int i = 0; i < nlead; i++)
    leaders[i]->init_leader();
}

void relay::start_relay(int gen_max_, bool verbose_)
{
  if (debugging_flag) stage_diagnostics(gen_max_);
  bool relay_underway=true;
  bool first2finish = true;

  printf("searching for first events\n");
  // finding first events
#pragma omp parallel
  {
    runner *rt = runners[thread_num()];
    rt->clear_event_data();
#pragma omp for nowait
    for (int i = 0; i < npool; i++)
    {
      rt->detect_events(pool[i], 0, Frames);
    }
#pragma omp critical
    {
      if (first2finish)
      {
        for (int i = 0; i < n*Frames; i++) global_event_frame_count[0][i] = rt->event_frame_count[0][i];
        first2finish=false;
      }
      else for (int i = 0; i < n*Frames; i++) global_event_frame_count[0][i] += rt->event_frame_count[0][i];
    }
  }

  ev->define_relay_leg(0);
  if (debugging_flag) rep->write_event_diagnostics(0);

  int gen_max_smooth = gen_max_/2, gen_max_stiff = gen_max_-gen_max_smooth;

  if (debugging_flag) rep->write_postevent_diagnostics(0);]
  train_leg(gen_max_smooth, gen_max_stiff);
  if (debugging_flag) rep->close_diagnostics(gen_count);
}
