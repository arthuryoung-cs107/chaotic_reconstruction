#include "particle_relay.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

extern "C"
{
  #include "AYaux.h"
}

relay::relay(referee &ref_,swirl_param &sp_min_,swirl_param &sp_max_,wall_list &wl_,double t_phys_, reporter *rep_, int ic_index_): referee(ref_), sp_min(sp_min_), sp_max(sp_max_), t_phys(t_phys_), rep(rep_), ic_index(ic_index_), n(rep_->P), Frames(rep_->Frames), ts(new double[Frames]), xs(new double[2*n*Frames]), d_ang(new double[Frames]), wl(wl_),
#ifdef _OPENMP
nt(omp_get_max_threads()), // each thread is a runner
#else
nt(1), // only one runner
#endif
pg(new proximity_grid*[nt]), rng(new AYrng*[nt]), runners(new runner*[nt])
{
  alloc_grades(nt, Frames);
  rep->load_filter(ts, xs, d_ang);



  // Set up each runner's personal data
#pragma omp parallel
  {
    int t=thread_num();
    pg[t]=new proximity_grid();
    rng[t]=new AYrng();
    rng[t]->rng_init_gsl(t+1);
    runners[t]=new runner(sp_min, pg[t], wl, n, t, param_len, Frames, sp_max.cl_im, t_phys, dt_sim, t_wheels, ts, xs, d_ang, ic_index, nlead, npool);
  }
}

relay::~relay()
{
  delete [] ts;
  delete [] xs;
  delete [] d_ang;
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

void relay::stage_diagnostics()
{
  rep->nlead = nlead;
  rep->npool = npool;
  rep->nA = nA;
  rep->param_len = param_len;
  rep->Frames = Frames;

  rep->lead_dup_count = lead_dup_count;
  rep->frame_kill_count = frame_kill_count;

  rep->sample_weights = sample_weights;
  rep->gen_frame_res_data = gen_frame_res_data;
  rep->gen_param_mean = gen_param_mean;
  rep->gen_param_var = gen_param_var;

  rep->leaders = leaders;

  rep->staged_flag = true;
}

void relay::init_relay()
{
  // initialize the pool of testing particles
  leader_count=gen_count=0;
  frscore_worst=1; l2score_worst=DBL_MAX;
#pragma omp parallel
  {
    AYrng *r = rng[thread_num()];
    double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
#pragma omp for
    for (int i = 0; i < npool; i++)
      pool[i]->init_pool(param_len, Frames, dmin, dmax, r);
  }
  for (int i = 0; i < nlead; i++)
    leaders[i]->init_leader(param_len, Frames);
}

void relay::start_relay(int gen_max_, bool verbose_)
{
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
        for (int i = 0; i < n*Frames; i++) global_event_frame_count[i] = rt->event_frame_count[i];
        first2finish=false;
      }
      else for (int i = 0; i < n*Frames; i++) global_event_frame_count[i] += rt->event_frame_count[i];
    }
  }

  int latest_event = 0;
  for (int i = 0; i < n; i++) for (int j = 0; j < Frames; j++)
    if (global_event_frame_count[i][j])
    {
      event_end[i] = j-1;
      if (event_end[i]>latest_event) latest_event = event_end[i];
      break;
    }

  resample_gen0_pool();

  // train off of the first events
  do
  {
    int success_local=0;
    for (int i = 0; i < 4*Frames; i++) gen_frame_res_data[i]=0.0;
    for (int i = 0; i < Frames; i++) frame_kill_count[i] = 0;
#pragma omp parallel
    {
      runner *rt = runners[thread_num()];
      rt->reset_diagnostics();
#pragma omp for reduction(+:success_local) nowait
      for (int i = 0; i < npool; i++)
      {
        success_local += rt->run_relay(pool[i], frscore_worst, l2score_worst);
      }
      rt->consolidate_diagnostics();
#pragma omp critical
      {
        rt->update_diagnostics(frame_kill_count,gen_frame_res_data);
      }
    }
    gen_count++;
    printf("gen %d: %d candidates. ", gen_count, success_local);
    if (check_pool_results()) training_underway=false; // we win
    else if (gen_count == gen_max_) training_underway=false; // we give up
    else resample_pool(); // we try again
    if (debugging_flag) rep->write_gen_diagnostics(gen_count, leader_count, worst_leader, best_leader, dup_count, resample_count, dup_unique, repl_count, t_wheels, min_res);
  } while (training_underway);
  if (debugging_flag) rep->close_diagnostics(gen_count, worst_leader, best_leader, t_wheels, min_res);
  printf("\n");
}

void relay::resample_gen0_pool()
{
  
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
    int worst_best = find_worst_grade(candidates, nlead);
    for (int i = nlead; i < pool_success_count; i++)
      if (candidates[worst_best]->isworse(candidates[i]))
      {
        candidates[worst_best] = candidates[i];
        worst_best = find_worst_grade(candidates, nlead);
      }
    return nlead;
  }
  else return pool_success_count;
}

bool relay::check_pool_results()
{
  pool_candidates = collect_candidates();
  // if we now have a full leader roster
  repl_count=0;
  if (leader_count + pool_candidates >= nlead)
  {
    // fill up remainder of leaders
    for (int i = 0, gap = nlead-leader_count; i < gap; i++)
    {
      leader_board[leader_count] = leaders[leader_count];
      leaders[leader_count++]->take_vals(candidates[--pool_candidates]);
    }

    int worst_best = find_worst_grade(leader_board, nlead);

    // consider the pool candidates, which are positioned adjacent to the current leaders on the leaderboard
    for (int i = nlead; i < nlead + pool_candidates; i++)
      if (leader_board[worst_best]->isworse(leader_board[i]))
      {
        leader_board[worst_best] = leader_board[i];
        worst_best = find_worst_grade(leader_board, nlead);
      }

    worst_leader = worst_best; // this will be the worst leader
    frscore_worst = leader_board[worst_best]->frscore;
    l2score_worst = leader_board[worst_best]->l2score;

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

  best_leader = find_best_grade(leaders, leader_count);
  l2score_best = leaders[best_leader]->l2score;
  frscore_best = leaders[best_leader]->frscore;
  printf("Best: (ID, gen, parents, frsc, l2sc) = (%d %d %d %d %e), %d replacements. ", best_leader, leaders[best_leader]->gen, leaders[best_leader]->parent_count, frscore_best, l2score_best, repl_count);
  if (frscore_best == Frames-1) return true;
  return false;
}

double relay::compute_leader_statistics()
{
  double nlead_inv = 1.0/((double)leader_count);
  min_res = DBL_MAX;
  for (int i = 0; i < param_len; i++) gen_param_mean[i] = gen_param_var[i] = 0.0;
  #pragma omp parallel
  {
    runner *rt = runners[thread_num()];
    double * param_mean = rt->param_mean;
    for (int i = 0; i < param_len; i++) param_mean[i] = 0.0;
    #pragma omp for reduction(min:min_res) nowait
    for (int i = 0; i < leader_count; i++)
    {
      double * params_i = leaders[i]->params;
      for (int j = 0; j < param_len; j++) param_mean[j]+=(nlead_inv)*params_i[j];
      if (leaders[i]->l2score < min_res) min_res = leaders[i]->l2score;
    }
    #pragma omp critical
    {
      for (int i = 0; i < param_len; i++) gen_param_mean[i] += param_mean[i];
    }

    for (int i = 0; i < param_len; i++) param_mean[i] = 0.0;

    #pragma omp barrier

    #pragma omp for nowait
    for (int i = 0; i < leader_count; i++)
    {
      double * params_i = leaders[i]->params;
      for (int j = 0; j < param_len; j++)
      {
        double z = (params_i[j] - gen_param_mean[j]);
        param_mean[j]+=nlead_inv*(z*z);
      }
    }
    #pragma omp critical
    {
      for (int i = 0; i < param_len; i++) gen_param_var[i] += param_mean[i];
    }
  }
  return min_res;
}

void relay::resample_pool()
{
  min_res = compute_leader_statistics();

  double acc = 0.0;
  for (int i = 0; i < leader_count; i++)
    acc += sample_weights[i] = leaders[i]->w(min_res);

  acc /= (leader_count<nlead)? rs_fill_factor:rs_full_factor;

  for (int i = 0; i < leader_count; i++)
    {sample_weights[i] /= acc; lead_dup_count[i] = 0;}

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
        double uni = r->rand_uni_gsl(0.0, 1.0);
        while ((j<leader_count)&&(uni>0.0)) uni -= sample_weights[j++];
        double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
        if (j > 0) // if we have particles worth resampling
        {
          // resample the particle, add gaussian noise
          if (uni<0.0)
          {dup_count++; dup_t[--j]++; pool[i]->duplicate(leaders[j], gen_count, dmin, dmax,r, gen_param_var);}
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

    printf("%d duplicates (%d unique), %d resamples\n", dup_count, dup_unique, resample_count);
}
