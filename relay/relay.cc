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
  rep->load_filter(ts, xs, d_ang);

  data_scale = sp_min_.cl_im;

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
  rep->nlead=nlead;
  rep->npool=npool;
  rep->param_len=param_len;
  rep->beads=n;

  rep->dt_sim = dt_sim;
  rep->noise_tol = noise_tol;
  rep->alpha_tol = alpha_tol;
  rep->t_phys = t_phys;

  rep->lead_dup_count = lead_dup_count;
  rep->global_event_frame_count = *global_event_frame_count;
  rep->event_end = event_end;
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
  rep->write_startup_diagnostics(gen_max_);
}

void relay::init_relay()
{
  // initialize the pool of testing particles
  leader_count=gen_count=0;
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
  learn_first_leg(gen_max_, verbose_);
}

void relay::learn_first_leg(int gen_max_, bool verbose_)
{
  bool training_underway=true;
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
  if (debugging_flag) rep->write_event_diagnostics(0);
  check_gen0();
  if (debugging_flag) rep->write_postevent_diagnostics(0);

  residual_worst = DBL_MAX;
  do // train off of the first events
  {
    int success_local=0;
#pragma omp parallel
    {
      runner *rt = runners[thread_num()];
#pragma omp for reduction(+:success_local) nowait
      for (int i = 0; i < npool; i++)
      {
        success_local += rt->run_relay(pool[i], 0, event_end, earliest_event, latest_event, residual_worst);
      }
    }
    gen_count++;
    printf("(gen %d): %d candidates. ", gen_count, success_local);
    if (check_pool_results()) training_underway=false; // we win
    else if (gen_count == gen_max_) training_underway=false; // we give up
    else resample_pool(); // we try again
    if (debugging_flag) rep->write_gen_diagnostics(gen_count, leader_count);
  } while (training_underway);
  if (debugging_flag) rep->close_diagnostics(gen_count);
  printf("\n");
}

void relay::check_gen0()
{
  printf("(gen 0): First events identified. Event frames -");
  latest_event = event_observations = 0;
  earliest_event = Frames;
  for (int i = 0; i < n; i++) for (int j = 0; j < Frames; j++)
    if (global_event_frame_count[i][j])
    {
      event_end[i] = j-1;
      event_observations+=2*event_end[i];
      if (event_end[i]>latest_event) latest_event = event_end[i];
      if (event_end[i]<earliest_event) earliest_event = event_end[i];
      printf(" %d", j-1);
      break;
    }
  // event_observations = n*2*earliest_event;
  printf(". Earliest: %d, latest: %d ", earliest_event, latest_event);

  tau = (noise_tol*data_scale)*sqrt((double)(event_observations));
  tau_sqr = tau*tau;

  if (gen0_resample_flag)
  {
    int nresamp = npool-nlead;
    for (int i = 0; i < nlead; i++) leader_board[i] = pool[i];
      int worst_best = find_worst_record(pool, nlead);
      for (int i = nlead; i < npool; i++)
      {
        if (leader_board[worst_best]->isworse(pool[i]))
        {
          candidates[i-nlead] = leader_board[worst_best];
          leader_board[worst_best] = pool[i];
          worst_best = find_worst_record(leader_board, nlead);
        }
        else candidates[i-nlead] = pool[i];
      }

      int best_surviving = find_best_record(leader_board, nlead),
          best_killed = find_best_record(candidates, nresamp),
          worst_killed = find_worst_record(candidates, nresamp),
          i_bs = leader_board[best_surviving]->global_index,
          i_ws = leader_board[worst_best]->global_index,
          i_bk = candidates[best_killed]->global_index,
          i_wk = candidates[worst_killed]->global_index;

      double  r_bs = leader_board[best_surviving]->residual,
              r_ws = leader_board[worst_best]->residual,
              r_bk = candidates[best_killed]->residual,
              r_wk = candidates[worst_killed]->residual;

    printf(". Resampling %d weakest particles - \n", nresamp);
    #pragma omp parallel
      {
        AYrng *r = rng[thread_num()];
        double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
    #pragma omp for
        for (int i = 0; i < nresamp; i++)
          candidates[i]->resample(0, dmin, dmax, r);
      }
      printf("Done. Residuals %e to %e (%d, %d) kept. Residuals %e to %e (%d, %d) resampled. Beginning training:\n", r_bs, r_ws, i_bs, i_ws, r_bk, r_wk, i_bk, i_wk);
  }
  else
  {
    int w_gen0 = find_worst_record(pool, npool),
        b_gen0 = find_best_record(pool, npool);
    printf(". Residual range: [%e,%e]. Skipping gen0 resampling. Beginning training\n", pool[b_gen0]->residual, pool[w_gen0]->residual);
  }
}

int relay::collect_candidates()
{
  pool_success_count = 0;
  for (int i = 0; i < npool; i++)
    if (pool[i]->success)
      candidates[pool_success_count++] = pool[i];
  if (pool_success_count > nlead) // if we have lots of good candidates
  {
    // we seek to pick the nlead best grades for comparison.
    int worst_best = find_worst_record(candidates, nlead);
    for (int i = nlead; i < pool_success_count; i++)
      if (candidates[worst_best]->isworse(candidates[i]))
      {
        candidates[worst_best] = candidates[i];
        worst_best = find_worst_record(candidates, nlead);
      }
    return nlead;
  }
  else return pool_success_count;
}

bool relay::check_pool_results()
{
  int pool_candidates_local = pool_candidates = collect_candidates();
  // if we now have a full leader roster
  repl_count=0;
  if (leader_count + pool_candidates_local >= nlead)
  {
    // fill up remainder of leaders
    for (int i = 0, gap = nlead-leader_count; i < gap; i++)
    {
      leader_board[leader_count] = leaders[leader_count];
      leaders[leader_count++]->take_vals(candidates[--pool_candidates_local]);
    }

    int worst_best = find_worst_record(leader_board, nlead);

    // consider the pool candidates, which are positioned adjacent to the current leaders on the leaderboard
    for (int i = nlead; i < nlead + pool_candidates_local; i++)
      if (leader_board[worst_best]->isworse(leader_board[i]))
      {
        leader_board[worst_best] = leader_board[i];
        worst_best = find_worst_record(leader_board, nlead);
      }

    worst_leader = worst_best; // this will be the worst leader
    residual_worst = leader_board[worst_best]->residual;

    double resacc_local = 0.0;
    for (int i = 0; i < nlead; i++)
    {
      resacc_local+=leader_board[i]->residual;
      if (leaders[i]->global_index != leader_board[i]->global_index)
      {
        repl_count++;
        leaders[i]->take_vals(leader_board[i]);
        leader_board[i] = leaders[i];
      }
    }
    pos_res_global=resacc_local/((double)(leader_count-1));
  }
  // otherwise, we can just fill in the leaderboard
  else for (int i = 0; i < pool_candidates; i++)
  {
    repl_count++;
    leader_board[leader_count] = leaders[leader_count];
    leaders[leader_count++]->take_vals(candidates[i]);
  }

  best_leader = find_best_record(leaders, leader_count);
  record * wl_rec = leaders[worst_leader], * bl_rec = leaders[best_leader];
  residual_best = bl_rec->residual;
  printf("Best/Worst: (ID, gen, parents, residual) = (%d/%d %d/%d %d/%d %e/%e), %d replacements. ", bl_rec->global_index, wl_rec->global_index, bl_rec->gen, wl_rec->gen, bl_rec->parent_count, wl_rec->parent_count, residual_best, residual_worst, repl_count);
  if (sqrt(residual_worst)<tau) return true;
  return false;
}
