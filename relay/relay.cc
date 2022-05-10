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
    runners[t]=new runner(sp_min, pg[t], wl, n, t, param_len, Frames, nlead, npool, t_phys, dt_sim, alpha_tol, ts, xs, d_ang, comega_s, t_wheels);
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
  dt_sim,
  noise_tol,
  alpha_tol,
  rs_full_factor,
  data_scale,
  t_phys
  };

  // integer pointers
  rep->lead_dup_count = lead_dup_count;
  rep->global_event_frame_count = *global_event_frame_count;
  rep->event_frames = event_frames;
  rep->gen_int_vec = &(gen_count);
  rep->event_int_vec = &(earliest_event);

  // double pointers
  rep->sample_weights = sample_weights;
  rep->lead_par_w_mean = lead_par_w_mean;
  rep->lead_par_w_var = lead_par_w_var;
  rep->gen_double_vec = &(residual_best);
  rep->event_double_vec = &(tau);

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
  // bool relay_underway=true;
  // int event_block=0;
  // int start_frame=0;
  // int gen_smooth_max = gen_max_/2;
  //
  // if (debugging_flag) stage_diagnostics(gen_max_);
  //
  // #pragma omp parallel for
  // for (int i = 0; i < nlead+npool; i++)
  // {
  //   leaders[i]->init_training(true);
  // }
  //
  // find_events(start_frame,Frames);
  // ev->define_relay_event_block(event_block, &net_observations, &tau_full, &earliest_event, noise_tol*data_scale);
  // tau=tau_smooth; tau_sqr=tau*tau;
  // if (debugging_flag) rep->write_event_diagnostics(event_block);
  //
  // // smooth training
  // start_frame = train_event_block(event_block, gen_smooth_max,0.005);
  // resample_pool();
  //
  // printf("(event block %d) Completed smooth training. Best residual: %e, for tau^2=%e. tau/r= %e. Beginning stiff training: \n", event_block, residual_best, tau_sqr, tau/root_res_best);
  // #pragma omp parallel for
  // for (int i = 0; i < nlead+npool; i++)
  // {
  //   leaders[i]->init_training(false);
  // }
  // // reevaluate leaders with their stiff performance.
  // assess_leaders();
  //
  // tau=tau_full; tau_sqr=tau*tau;
  //
  // // full block training
  // int block_end = train_event_block(event_block, gen_max_,0.005);
  // resample_pool();
  // printf("(event block %d) Completed stiff training. Best residual: %e, for tau^2=%e. tau/r= %e.\n", event_block, residual_best, tau_sqr, tau/root_res_best);
  //
  // for (int i = 0; i < nlead; i++)
  // {
  //   pool[i]->take_vals(leaders[i]);
  //   leaders[i]->residual = DBL_MAX;
  // }
  // tau=noise_tol*data_scale*sqrt((double)(2*n*Frames)); tau_sqr=tau*tau;
  // residual_worst=DBL_MAX;
  //
  // printf("Reloading leaders. Training on full data.\n");
  //
  // block_end = train_event_block(event_block, gen_max_,0.01,true);
  //
  // if (debugging_flag) rep->close_diagnostics(gen_count);
}

void relay::start_block_relay(int gen_max_, bool verbose_)
{
  bool relay_underway=true;
  int event_block=0;
  int start_frame=0;

  for (int i = 0; i < n; i++) event_frames[i] = 0;

  if (debugging_flag) stage_diagnostics(gen_max_);

  #pragma omp parallel for
  for (int i = 0; i < nlead+npool; i++)
  {
    leaders[i]->init_training(true);
  }
  find_events(0,Frames);
  ev->define_relay_event_block(event_block, &net_observations, &tau_full, &earliest_event, noise_tol*data_scale);
  tau=tau_smooth; tau_sqr=tau*tau;
  if (debugging_flag) rep->write_event_diagnostics(event_block, npool, pool);
  do
  {
    // smooth training
    start_frame = train_event_block(event_block, 50,0.05);
    printf("\n(event block %d) Completed smooth training. Best residual: %e, for tau^2=%e. tau/r= %e. Beginning stiff training: \n", event_block, residual_best, tau_sqr, tau/root_res_best);
    reload_leaders(false);

    tau=tau_full; tau_sqr=tau*tau;

    // full block training
    int block_end = train_event_block(event_block, 50,0.1);
    printf("\n(event block %d) Completed stiff training. Best residual: %e, for tau^2=%e. tau/r= %e.", event_block, residual_best, tau_sqr, tau/root_res_best);
    reload_leaders(true);

    event_block++;
    find_events(0,Frames, true);
    bool same_frames = ev->define_relay_event_block(event_block, &net_observations, &tau_full, &earliest_event, noise_tol*data_scale);

    if (gen_count==gen_max_) relay_underway=false;
    else if (ev->check_last()) relay_underway=false;
    else
    {
      tau=tau_smooth; tau_sqr=tau*tau;
      if (debugging_flag) rep->write_event_diagnostics(event_block, nlead, leaders);
    }
  } while(relay_underway);
  if (debugging_flag) rep->close_diagnostics(gen_count);
}

void relay::assess_leaders()
{
  worst_leader = find_worst_record(leaders, leader_count);
  best_leader = find_best_record(leaders, leader_count);
  record * wl_rec = leaders[worst_leader], * bl_rec = leaders[best_leader];
  residual_best = bl_rec->residual; root_res_best = bl_rec->root_residual;
  residual_worst = wl_rec->residual; root_res_worst = wl_rec->root_residual;
}

void relay::reload_leaders(bool smooth_training_)
{
  #pragma omp parallel for
  for (int i = 0; i < nlead+npool; i++)
  {
    leaders[i]->init_training(smooth_training_);
  }
  if (smooth_training_) assess_leaders(); // reevaluate leaders on their stiff performance
  else
  {
    for (int i = 0; i < nlead; i++)
    {
      pool[i]->take_vals(leaders[i]);
      leaders[i]->residual=DBL_MAX;
    }
    residual_worst=DBL_MAX;
  }
}
