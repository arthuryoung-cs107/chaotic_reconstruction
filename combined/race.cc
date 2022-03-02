#include "particle_race.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

race::race(referee &ref_,swirl_param &sp_min_,swirl_param &sp_max_,wall_list &wl_,double t_phys_, ODR_struct &odr_, int ic_index_): referee(ref_), sp_min(sp_min_), sp_max(sp_max_), sp_rnd(sp_rnd_), t_phys(t_phys_), odr(odr_), ic_index(ic_index_), n(odr_.P), nsnap(odr_.Frames), ts(new double[nsnap]), xs(new double[2*n*nsnap]), d_ang(new double[nsnap]),
#ifdef _OPENMP
nt(omp_get_max_threads()), // each thread is a runner
#else
nt(1), // only one runner
#endif
wl(wl_), pg(new proximity_grid*[nt]), rng(new gsl_rng*[nt]), runners(new runner*[nt])
{
  alloc_records();
  odr.load_filter(ts, xs, d_ang);

  double *x_ic = xs + 2*n*ic_index;
  double t_ic=ts[ic_index]/t_phys;
  double comega_ic = d_ang[ic_index];
  // Set up each runner's personal data
#pragma omp parallel
  {
      int t=thread_num();
      pg[t]=new proximity_grid();
      rng[t]=gsl_rng_alloc(gsl_rng_taus2);
      gsl_rng_set(rng[t], t+1);
      runners[t]=new runner(sp_min, pg[t], wl, n, t, param_len, nsnap, sp_max.cl_im);
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
    gsl_rng *uni = rng[thread_num()];
#pragma omp for
    for (int i = 0; i < npool; i++)
    {
      double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
      for (int j = 0; j < param_len; j++) pool_params[i][j] = dmin[j] + (dmax[j]- dmin[j])*(gsl_rng_uniform(uni));
      pool[i].index = i; pool[i].params = pool_params[i];
    }
  }

  for (int i = 0; i < nlead; i++)
  {leaders[i].index = i; leaders[i].params = leaders[i];}

  leader_count=gen_count=0;
  frscore_min=1; l2score_min=DBL_MAX;
}

void race::run()
{
  bool race_underway=true;
  do
  {
    pool_success_count = 0;
#pragma omp parallel
    {
      runner *rt = runners[thread_num()];
#pragma omp for reduction(+=)
      for (int i = 0; i < npool; i++)
      {
        rt->reset_sim(pool_params[i]);
        rt->run_race(ts, xs, d_ang);
        pool_success_count += pool_ranking[i] = pool[i].check_success(rt->frame, rt->pos_err_acc, frscore_min, l2score_min);
      }
    }
    gen_count++;
    if (check_pool_results()) race_underway=false; // we win
    else if (gen_count == gen_max) race_underway=false; // we give up
    else resample_pool(); // we try again

  } while (race_underway);
}

bool race::check_pool_results()
{
  collect_pool_leaders(); 

  if (leader_count == nlead)
  {

  }
  else if (pool_success_count >= (nlead-leader_count))
  {
    int leader_gap = nlead-leader_count;
    // fill to capacity
    for (int i = 0; i < leader_gap; i++)
    {

    }
    leader_count = nlead;

  }
  else // just enter the results. No need to worry about overflow
  { // note that this conditional also handles our startup
    // loop over the successful entries, record their scores, enter them into the leaderboard
    for (int i = 0; i < pool_success_count; i++)
    {
      for (int j = 0; j < param_len; j++) leaders[leader_count][j] = pool[sorting_pool[i]][j];
      frame_score[leader_count] = frame_pscore[sorting_pool[i]];
      l2_score[leader_count++] = l2_pscore[sorting_pool[i]];
    }
    rank_leaders();
  }
}
