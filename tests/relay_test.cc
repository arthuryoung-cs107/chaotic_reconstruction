#include <sys/types.h>
#include <sys/stat.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "particle_relay.hh"

const int nbeads=3;
const int param_id=0;
const int relay_id=0;
const int generations=10;
const int nlead = 500;
const int npool = 1000;
const int param_len = 12;
const double dt_sim = 0.002;
const double t_wheels = 0.012; // should we crank this up?
const double noise_tol = 1e-8;
const double weight_ceiling = 1e10;
const double alpha_tol=100.0;
char proc_loc[] = "./dat_dir/";
int main()
{
  assert(nbeads<=30);
  char rydat_dir[30]; sprintf(rydat_dir, "swirl%d.odr/", nbeads);
  char file_name[10]; sprintf(file_name, "pts.%d", param_id);
  char swbest_name[30]; sprintf(swbest_name, "swirl_best_wa%d.odr/", relay_id);

  // Physical constants
  double g_phys=9.804,                 // Gravity (m/s^2)
         d_phys=0.00635,               // Diameter (m)
         t_phys=sqrt(d_phys/g_phys);   // Time unit (s)

  double  sp_min_vals[14], sp_max_vals[14];
  int idmin=set_special_params("min", sp_min_vals),
  idmax=set_special_params("max", sp_max_vals);

  swirl_param sp_min(sp_min_vals), sp_max(sp_max_vals);

  // Create the hexagonal dish
  wall_list wl;
  const double r=5.72;
  const double fa=sqrt(0.75);
  wall_floor wf(0.);
  wall_par_planes wp0(0,1,0,r),wp1(fa,0.5,0,r),wp2(fa,-0.5,0,r);
  wl.add_wall(&wf);
  wl.add_wall(&wp0);
  wl.add_wall(&wp1);
  wl.add_wall(&wp2);

  reporter rep;
  rep.init_relay(proc_loc, rydat_dir, file_name, relay_id);

  referee ref(nlead, npool, param_len, dt_sim, noise_tol, alpha_tol, weight_ceiling);
  relay prelay(ref, sp_min,sp_max,wl,t_phys,&rep);
  prelay.init_relay();
  prelay.start_relay(generations);
  // prelay.make_best_swirl(swbest_name);

  return 0;
}
