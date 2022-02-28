#include "race.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

race::race(race_param &rparam,swirl_param &sp_min_,swirl_param &sp_max_,swirl_param &sp_rnd_,wall_list &wl_,double t_phys_, ODR_struct * odr_,int offset): race_param(rparam), sp_min(sp_min_), sp_max(sp_max_), sp_rnd(sp_rnd_),
t_phys(t_phys_), min_l2(0.), min_linf(0.), nfail(0), fflags(0), odr(odr_),
#ifdef _OPENMP
nt(omp_get_max_threads()),
#else
nt(1),
#endif
ttab(new int[nt+1]), rloc(new int[nt+1]), rng(new gsl_rng*[nt]), wl(wl_),
pg(new proximity_grid*[nt]), odir(NULL), fdigest(NULL)
{
  n = odr->P;
  nsnap = odr->Frames - offset;
  ts=new double[nsnap];
  xs=new double[2*n*nsnap];
  d_ang=new double[nsnap];

  odr->load_filter(ts, xs, d_ang, offset);

  // Set up the random number generators
#pragma omp parallel
  {
      int t=thread_num();
      rng[t]=gsl_rng_alloc(gsl_rng_taus2);
      pg[t]=new proximity_grid();
      gsl_rng_set(rng[t],t+1);
  }

  // Set the zeroth entries of the thread table and resampling location table
  // to zero, since they are always held at this value
  *ttab=*rloc=0;
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
