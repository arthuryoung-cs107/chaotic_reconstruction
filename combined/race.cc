#include "race.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

race::race(race_param &rparam,swirl_param &sp_min_,swirl_param &sp_max_,swirl_param &sp_rnd_,wall_list &wl_,double t_phys_, ODR_struct &odr_): race_param(rparam), sp_min(sp_min_), sp_max(sp_max_), sp_rnd(sp_rnd_),
t_phys(t_phys_), odr(odr_), n(odr_.P), nsnap(odr_.Frames), ts(new double[nsnap]), xs(new double[2*n*nsnap]), d_ang(new double[nsnap]),
#ifdef _OPENMP
nt(omp_get_max_threads()),
#else
nt(1),
#endif
rng(new gsl_rng*[nt]), wl(wl_), pg(new proximity_grid*[nt])
{
  odr.load_filter(ts, xs, d_ang, offset);
  // Set up the random number generators
#pragma omp parallel
  {
      int t=thread_num();
      rng[t]=gsl_rng_alloc(gsl_rng_taus2);
      pg[t]=new proximity_grid();
      gsl_rng_set(rng[t],t+1);
  }
}

void race::init()
{

}

void race::run(int frames_)
{
  bool race_underway = true;
  while (race_underway)
  {
    #pragma omp parallel
    {



    }
  }
}
