#include "particle_race.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

extern "C"
{
  #include "AYaux.h"
}

race::race(referee &ref_,swirl_param &sp_min_,swirl_param &sp_max_,wall_list &wl_,double t_phys_, reporter *odr_, int ic_index_): referee(ref_), sp_min(sp_min_), sp_max(sp_max_), t_phys(t_phys_), odr(odr_), ic_index(ic_index_), n(odr_->P), Frames(odr_->Frames), ts(new double[Frames]), xs(new double[2*n*Frames]), d_ang(new double[Frames]), wl(wl_),
#ifdef _OPENMP
nt(omp_get_max_threads()), // each thread is a runner
#else
nt(1), // only one runner
#endif
pg(new proximity_grid*[nt]), rng(new AYrng*[nt]), runners(new runner*[nt])
{
  alloc_records();
  sample_weights = new double[nlead];
  dup_vec = new int[nlead];
  dup_mat = AYimatrix(nt, nlead);
  odr->load_filter(ts, xs, d_ang);
  double *x_ic = xs + 2*n*ic_index;
  double t_ic=ts[ic_index];
  double comega_ic = d_ang[ic_index];
  // Set up each runner's personal data
#pragma omp parallel
  {
      int t=thread_num();
      pg[t]=new proximity_grid();
      rng[t]= new AYrng();
      rng[t]->rng_init_gsl(t+1);
      runners[t]=new runner(sp_min, pg[t], wl, n, t, param_len, Frames, sp_max.cl_im);
      runners[t]->init_ics(t_phys, x_ic, t_ic, comega_ic);
  }
}

race::~race()
{
  delete [] dup_vec;
  free_AYimatrix(dup_mat);

  delete [] ts;
  delete [] xs;
  delete [] d_ang;
  delete [] sample_weights;
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

void race::init_race()
{
// initialize the pool of testing particles
  leader_count=gen_count=0;
  frscore_min=1; l2score_min=DBL_MAX;
#pragma omp parallel
  {
    AYrng *r = rng[thread_num()];
    double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
#pragma omp for
    for (int i = 0; i < npool; i++)
      pool[i]->init_pool(pool_params[i], param_len, Frames, gau_var_high, gau_lambda, dmin, dmax, r);
  }
  for (int i = 0; i < nlead; i++)
    leaders[i]->init_leader(lead_params[i], param_len, Frames, gau_var_high, gau_lambda);
}

void race::start_race(int gen_max_, bool verbose_)
{
  if (debugging_flag) stage_diagnostics();
  bool race_underway=true;
  do
  {
    int success_local=0;
#pragma omp parallel
    {
      runner *rt = runners[thread_num()];
#pragma omp for schedule(dynamic) reduction(+:success_local)
      for (int i = 0; i < npool; i++)
      {
        rt->reset_sim(pool_params[i]);
        success_local += rt->run_race(dt_sim, ts, xs, d_ang, pool[i], frscore_min, l2score_min);
      }
    }
    gen_count++;
    printf("gen %d: %d candidates. ", gen_count, success_local);
    if (check_pool_results()) race_underway=false; // we win
    else if (gen_count == gen_max_) race_underway=false; // we give up
    else resample_pool(); // we try again
    if (debugging_flag) odr->write_gen_diagnostics(gen_count, leader_count, worst_leader, best_leader);
  } while (race_underway);
  if (debugging_flag) odr->close_diagnostics(gen_count, leader_count, worst_leader, best_leader);
  printf("\n");
}

bool race::check_pool_results()
{
  pool_candidates = collect_pool_leaders();
  // if we now have a full leader roster
  int repl_count=0;
  if (leader_count + pool_candidates >= nlead)
  {
    // fill up remainder of leaders
    for (int i = 0, gap = nlead-leader_count; i < gap; i++)
    {
      leader_board[leader_count] = leaders[leader_count];
      leaders[leader_count++]->take_vals(pool_leaders[--pool_candidates]);
    }

    int worst_best = find_worst(leader_board, nlead);

    // consider the pool candidates, which are positioned adjacent to the current leaders on the leaderboard
    for (int i = nlead; i < nlead + pool_candidates; i++)
      if (leader_board[worst_best]->isworse(leader_board[i]))
      {
        leader_board[worst_best] = leader_board[i];
        worst_best = find_worst(leader_board, nlead);
      }

    worst_leader = worst_best; // this will be the worst leader
    frscore_min = leader_board[worst_best]->frscore;
    l2score_min = leader_board[worst_best]->l2score;

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
    leaders[leader_count++]->take_vals(pool_leaders[i]);
  }

  best_leader = find_best(leaders, leader_count);
  l2score_best = leaders[best_leader]->l2score;
  frscore_best = leaders[best_leader]->frscore;
  printf("Best: (ID, gen, parents, frsc, l2sc) = (%d %d %d %d %e), %d replacements. ", best_leader, leaders[best_leader]->gen, leaders[best_leader]->parent_count, frscore_best, l2score_best, repl_count);
  if (frscore_best == Frames-1) return true;
  return false;
}

int race::collect_pool_leaders()
{
  pool_success_count = 0;
  for (int i = 0; i < npool; i++)
    if (pool[i]->success)
      pool_leaders[pool_success_count++] = pool[i];
  if (pool_success_count > nlead) // if we have lots of good candidates
  {
    // we seek to pick the nlead best records for comparison.
    int worst_best = find_worst(pool_leaders, nlead);
    for (int i = nlead; i < pool_success_count; i++)
      if (pool_leaders[worst_best]->isworse(pool_leaders[i]))
      {
        pool_leaders[worst_best] = pool_leaders[i];
        worst_best = find_worst(pool_leaders, nlead);
      }
    return nlead;
  }
  else return pool_success_count;
}

void race::resample_pool()
{
  double acc = 0.0;
  if (z_weight_flag)
  {
    double z_acc=0.0;
    for (int i = 0; i < leader_count; i++)
      z_acc+=leaders[i]->z_eval(frscore_min);
    double  z_mean=z_acc/((double)leader_count),
            lambda_z=1.0/z_mean;
    for (int i = 0; i < leader_count; i++)
      acc += sample_weights[i] = leaders[i]->w(lambda_z);
  }
  else for (int i = 0; i < leader_count; i++)
    acc += sample_weights[i] = leaders[i]->w(Frames, lambda);


  acc /= (leader_count<nlead)? rs_fill_factor:rs_full_factor;

  for (int i = 0; i < leader_count; i++)
    {sample_weights[i] /= acc; dup_vec[i] = 0;}

  int dup_count=0, res_count=0;
  #pragma omp parallel
    {
      int t = thread_num();
      AYrng *r = rng[t];
      int *dup_t = dup_mat[t];
      for (int i = 0; i < leader_count; i++) dup_t[i] = 0;
  #pragma omp for reduction(+:dup_count) reduction(+:res_count)
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
          {dup_count++; dup_t[--j]++; pool[i]->duplicate(leaders[j], gen_count, dmin, dmax,r);}
          // we hit the resampling pool
          else
            {res_count++; pool[i]->resample(gen_count, dmin, dmax, r);}
        }
        // we currently have no leaders (particles worth resampling)
        else
          {res_count++; pool[i]->resample(gen_count, dmin, dmax, r);}
      }
    }
    for (int i = 0; i < nt; i++) for (int j = 0; j < leader_count; j++)
      dup_vec[j]+=dup_mat[i][j];

    int dup_unique=0;
    for (int i = 0; i < leader_count; i++) if (dup_vec[i])
    {leaders[i]->dup_count+=dup_vec[i]; dup_unique++;}

    printf("%d duplicates (%d unique), %d resamples\n", dup_count, dup_unique, res_count);
}
