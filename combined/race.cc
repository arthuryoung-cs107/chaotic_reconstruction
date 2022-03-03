#include "particle_race.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

race::race(referee &ref_,swirl_param &sp_min_,swirl_param &sp_max_,wall_list &wl_,double t_phys_, ODR_struct &odr_, int ic_index_): referee(ref_), sp_min(sp_min_), sp_max(sp_max_), t_phys(t_phys_), odr(odr_), ic_index(ic_index_), n(odr_.P), Frames(odr_.Frames), ts(new double[Frames]), xs(new double[2*n*Frames]), d_ang(new double[Frames]),
#ifdef _OPENMP
nt(omp_get_max_threads()), // each thread is a runner
#else
nt(1), // only one runner
#endif
wl(wl_), pg(new proximity_grid*[nt]), rng(new AYrng*[nt]), runners(new runner*[nt])
{
  alloc_records();
  sample_weights = new double[nlead];
  odr.load_filter(ts, xs, d_ang);

  double *x_ic = xs + 2*n*ic_index;
  double t_ic=ts[ic_index]/t_phys;
  double comega_ic = d_ang[ic_index];
  // Set up each runner's personal data
#pragma omp parallel
  {
      int t=thread_num();
      pg[t]=new proximity_grid();
      rng[t]= new AYrng();
      rng[t]->rng_init_gsl(t+1);
      runners[t]=new runner(sp_min, pg[t], wl, n, t, param_len, Frames, sp_max.cl_im);
      runners[t]->init_ics(x_ic, t_ic, comega_ic);
  }
}

race::~race()
{

}

void race::init_race()
{
// initialize the pool of testing particles
#pragma omp parallel
  {
    AYrng *r = rng[thread_num()];
#pragma omp for
    for (int i = 0; i < npool; i++)
    {
      double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
      for (int j = 0; j < param_len; j++) pool_params[i][j] = dmin[j] + (dmax[j]- dmin[j])*(r->rand_uni_gsl(0.0, 1.0));
      pool[i]->params = pool_params[i];
    }
  }

  for (int i = 0; i < nlead; i++) leaders[i]->params = lead_params[i];

  leader_count=gen_count=0;
  frscore_min=1; l2score_min=DBL_MAX;
}

void race::start_race(int gen_max_, bool verbose_)
{
  bool race_underway=true;
  do
  {
#pragma omp parallel
    {
      runner *rt = runners[thread_num()];
#pragma omp for
      for (int i = 0; i < npool; i++)
      {
        rt->reset_sim(pool_params[i]);
        rt->run_race(t_phys, dt_sim, ts, xs, d_ang);
        pool[i]->check_success(rt->frame, rt->pos_err_acc, frscore_min, l2score_min);
      }
    }
    gen_count++;
    if (check_pool_results()) race_underway=false; // we win
    else if (gen_count == gen_max_) race_underway=false; // we give up
    else resample_pool(); // we try again

  } while (race_underway);
}

bool race::check_pool_results()
{
  int pool_candidates = collect_pool_leaders();
  // if we now have a full leader roster
  if (leader_count + pool_candidates >= nlead)
  {
    // fill up remainder of leaders
    for (int i = 0; i < (nlead-leader_count); i++)
      leader_board[leader_count++] = pool_leaders[--pool_candidates];

    int worst_best = find_worst(leader_board, nlead);

    for (int i = nlead; i < nlead + pool_candidates; i++)
      if (leader_board[worst_best]->isworse(leader_board[i]))
      {
        leader_board[worst_best] = leader_board[i];
        worst_best = find_worst(leader_board, nlead);
      }

    frscore_min = leader_board[worst_best]->frscore;
    l2score_min = leader_board[worst_best]->l2score;

    for (int i = 0; i < nlead; i++)
      if (leaders[i]->global_index != leader_board[i]->global_index)
        leaders[i]->take_vals(leader_board[i], param_len);
  }
  // otherwise, we can just fill in the leaderboard
  else for (int i = 0; i < pool_candidates; i++)
  {
    leader_board[leader_count] = leaders[leader_count];
    leaders[leader_count++]->take_vals(pool_leaders[i], param_len);
  }

  for (int i = 0; i < nlead; i++)
    if (leaders[i]->frscore == Frames - 1) return true;

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
  for (int i = 0; i < leader_count; i++)
    acc += sample_weights[i] = exp(lambda*((double)(leaders[i]->frscore-Frames+1)));

  acc /= (leader_count<nlead)? rs_fill_factor:rs_full_factor;

  for (int i = 0; i < leader_count; i++) sample_weights[i] /= acc;

  #pragma omp parallel
    {
      AYrng *r = rng[thread_num()];
  #pragma omp for
      for (int i = 0; i < npool; i++)
      {
        int j = 0;
        double uni = r->rand_uni_gsl(0.0, 1.0);
        while ( (j<leader_count)&&(uni>0.0) ) uni -= sample_weights[j++];
        double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
        if (j > 0) // if we have particles worth resampling
        {
          // resample the particle, add gaussian noise
          if (uni<0.0)
          {
            j--;
            for (int k = 0; k < param_len; k++)
            {
              /* should we be making this relative to the width of the gap?
              moreover, should we be scaling the gaussian variance by depth into the frames? */
              pool_params[i][k] = leaders[j]->params[k] + r->rand_gau_gsl(0.0, gau_var)*(dmax[k]-dmin[k]);
              if (pool_params[i][k] > dmax[k]) pool_params[i][k] = dmax[k];
              else if (pool_params[i][k] < dmin[k]) pool_params[i][k] = dmin[k];
            }
          }
          // we hit the resampling pool
          else for (int k = 0; k < param_len; k++)
            pool_params[i][k] = dmin[k] + (dmax[k]- dmin[k])*(r->rand_uni_gsl(0.0, 1.0));
        }
        // we currently have no leaders (particles worth resampling)
        else for (int k = 0; k < param_len; k++)
          pool_params[i][k] = dmin[k] + (dmax[k]- dmin[k])*(r->rand_uni_gsl(0.0, 1.0));
      }
    }
}
