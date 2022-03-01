#include "particle_race.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

race::race(referee &ref_,swirl_param &sp_min_,swirl_param &sp_max_,wall_list &wl_,double t_phys_, ODR_struct &odr_): referee(ref_), sp_min(sp_min_), sp_max(sp_max_), sp_rnd(sp_rnd_),
t_phys(t_phys_), odr(odr_), n(odr_.P), nsnap(odr_.Frames), ts(new double[nsnap]), xs(new double[2*n*nsnap]), d_ang(new double[nsnap]),
#ifdef _OPENMP
nt(omp_get_max_threads()), // each thread is a runner
#else
nt(1), // only one runner
#endif
wl(wl_), pg(new proximity_grid*[nt]), rng(new gsl_rng*[nt]), run(new runner*[nt])
{
  alloc_records();
  odr.load_filter(ts, xs, d_ang);

  // Set up each runner's personal data
#pragma omp parallel
  {
      int t=thread_num();
      pg[t]=new proximity_grid();
      rng[t]=gsl_rng_alloc(gsl_rng_taus2);
      gsl_rng_set(rng[t], t+1);
      run[t]=new runner(sp_min, pg[t], wl, n, t);
  }
}

void race::init_race()
{

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < npool; i++)
    {
      double dmin=&sp_min.Kn, dmax=&sp_max.Kn;
      for (int j = 0; j < param_len; j++) 
    }

  }

}

void race::run()
{

  startup();

  bool race_underway = true;
  while (race_underway)
  {
    #pragma omp parallel
    {

    }
  }
}
