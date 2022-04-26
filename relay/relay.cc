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

  // initialize the dish velocity data. Will save time on calculations later
  comega_s[0]=0.0;
  for (int i = 1; i < Frames; i++)
  {
    double comega=d_ang[i]-d_ang[i-1];
    if(comega>M_PI) comega-=2*M_PI; else if(comega<-M_PI) comega+=2*M_PI; comega/=dur;
    comega_s[i] = t_phys*comega/(ts[i]-ts[i-1]);
  }

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
}

void relay::stage_diagnostics(int gen_max_)
{
  rep_->nlead=nlead;
  rep_->npool=npool;
  rep_->param_len=param_len;
  rep_->beads=n;
  rep_->Frames=Frames;

  rep_->dt_sim = dt_sim;
  rep_->noise_tol = noise_tol;
  rep_->alpha_tol = alpha_tol;
  rep_->max_weight_ceiling = max_weight_ceiling;
  rep_->t_phys = t_phys;

  rep_->lead_dup_count = lead_dup_count;
  rep_->global_event_frame_count = *global_event_frame_count;
  rep_->event_end = event_end;
  rep_->gen_int_vec = &(gen_count);
  rep_->postevent_int_vec = &(event_observations);

  rep_->sample_weights = sample_weights;
  rep_->lead_par_w_mean = lead_par_w_mean;
  rep_->lead_par_w_var = lead_par_w_var;
  rep_->gen_double_vec = &(residual_best);
  rep_->postevent_double_vec = &(gau_scale_sqrt);

  rep_->leaders = leaders;
  rep_->pool = pool;

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
  gen_max=gen_max_;
  if (debugging_flag) stage_diagnostics();
  bool relay_underway=true;
  learn_first_leg(gen_max_, verbose_);
}

void relay::learn_first_leg(int gen_max_, bool verbose_)
{
  bool training_underway=true;
  bool first2finish = true;

  // finding first events
#pragma omp parallel
  {
    runner *rt = runners[thread_num()];
    rt->clear_event_data();
#pragma omp for reduction(+:success_local) nowait
    for (int i = 0; i < npool; i++)
    {
      rt->detect_events(pool[i), 0, Frames);
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
  gau_scale_sqrt = noise_tol*sqrt((double)(event_observations-1));
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
        success_local += rt->run_relay(pool[i], 0, event_end, latest_event, residual_worst);
      }
    }
    gen_count++;
    printf("gen %d: %d candidates. ", gen_count, success_local);
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
  int nresamp = npool-nlead;
  printf("(gen 0: First events identified. Event frames - ");
  latest_event = event_observations = 0;
  for (int i = 0; i < n; i++) for (int j = 0; j < Frames; j++)
    if (global_event_frame_count[i][j])
    {
      event_end[i] = j-1;
      event_observations+=2*event_end[i];
      if (event_end[i]>latest_event) latest_event = event_end[i];
      printf("%d ", j-1);
      break;
    }
  printf("\n. Resampling %d weakest particles - \n", nresamp);

  for (int i = 0; i < nlead; i++) leader_board[i] = pool[i];
  int worst_best = find_worst_record(pool, nlead);
  for (int i = nlead; i < nresamp; i++)
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

double relay::compute_leader_statistics()
{
  double res_scale = sqrt(residual_best)/gau_scale_sqrt;
  max_weight_factor = log(max_weight_ceiling/((double)(leader_count)))+0.5*(res_scale)*(res_scale);

  double w_sum = 0.0;
  for (int i = 0; i < leader_count; i++)
  {
    w_sum += sample_weights[i] = leaders[i]->w(max_weight_factor,gau_scale_sqrt);
    lead_dup_count[i] = 0;
  }

  for (int i = 0; i < param_len; i++) lead_par_w_mean[i] = lead_par_w_var[i] = 0.0;
  #pragma omp parallel
  {
    runner *rt = runners[thread_num()];
    double * param_buf = rt->param_acc;
    for (int i = 0; i < param_len; i++) param_buf[i] = 0.0;
    #pragma omp for nowait
    for (int i = 0; i < leader_count; i++)
    {
      for (int j = 0; j < param_len; j++) param_buf[j]+=(sample_weights[i]*(leaders[i]->params[j]))/w_sum;
    }
    #pragma omp critical
    {
      for (int i = 0; i < param_len; i++) lead_par_w_mean[i] += param_buf[i];
    }

    for (int i = 0; i < param_len; i++) param_buf[i] = 0.0;

    #pragma omp barrier

    #pragma omp for nowait
    for (int i = 0; i < leader_count; i++)
    {
      double * params_i = leaders[i]->params;
      for (int j = 0; j < param_len; j++)
      {
        double z = (params_i[j] - lead_par_w_mean[j]);
        param_buf[j]+=(sample_weights[i])*(z*z)/(w_sum);
      }
    }
    #pragma omp critical
    {
      for (int i = 0; i < param_len; i++) lead_par_w_var[i] += param_buf[i];
    }
  }
  return w_sum;
}

void relay::resample_pool()
{
  double w_sum = compute_leader_statistics();
  dup_count=resample_count=dup_unique=0;
  #pragma omp parallel
    {
      int t = thread_num();
      AYrng *r = rng[t];
      int *dup_t = runners[t]->lead_dup_count;
      for (int i = 0; i < leader_count; i++) dup_t[i] = 0;
      #pragma omp for reduction(+:dup_count) reduction(+:resample_count) nowait
      for (int i = 0; i < npool; i++)
      {
        int j = 0;
        double uni = w_sum*(r->rand_uni_gsl(0.0, 1.0));
        while ((j<leader_count)&&(uni>0.0)) uni -= sample_weights[j++];
        double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
        if (j > 0) // if we have particles worth resampling
        {
          // resample the particle, add gaussian noise
          if (uni<0.0)
          {dup_count++; dup_t[--j]++; pool[i]->duplicate(leaders[j], gen_count, dmin, dmax,r, lead_par_w_var);}
          // we hit the resampling pool
          else
            {resample_count++; pool[i]->resample(gen_count, dmin, dmax, r);}
        }
        // we currently have no leaders (particles worth resampling)
        else
          {resample_count++; pool[i]->resample(gen_count, dmin, dmax, r);}
      }
      #pragma omp critical
      {
        for (int i = 0; i < leader_count; i++) lead_dup_count[i]+=dup_t[i];
      }
    }

    for (int i = 0; i < leader_count; i++) if (lead_dup_count[i])
    {leaders[i]->dup_count+=lead_dup_count[i]; dup_unique++;}

    printf("%d duplicates (%d unique)\n", dup_count, dup_unique);
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

    for (int i = 0; i < nlead; i++)
      if (leaders[i]->global_index != leader_board[i]->global_index)
      {
        repl_count++;
        leaders[i]->take_vals(leader_board[i]);
        leader_board[i] = leaders[i];
      }
  }
  // otherwise, we can just fill in the leaderboard
  else for (int i = 0; i < pool_candidates; i++)
  {
    repl_count++;
    leader_board[leader_count] = leaders[leader_count];
    leaders[leader_count++]->take_vals(candidates[i]);
  }

  best_leader = find_best_record(leaders, leader_count);
  residual_best = leaders[best_leader]->residual;
  printf("Best: (ID, gen, parents, residual) = (%d %d %d %e), %d replacements. ", best_leader, leaders[best_leader]->gen, leaders[best_leader]->parent_count, residual_best, repl_count);
  if (sqrt(residual_best)<gau_scale_sqrt) return true;
  return false;
}
