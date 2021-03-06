#include <sys/types.h>
#include <sys/stat.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "particle_relay.hh"

const int nbeads=3;
const int param_id=0;
const int relay_id=1;
const int generations=200;
const int nlead = 1;
const int npool = 2000;
const int param_len = 12;
const double dt_sim = 0.002;
// const double t_wheels = 0.012;
const double t_wheels = 0.0;
const double noise_tol = 1e-2;
const double alpha_tol=10.0;
const double rs_full_factor=1.0;
const bool test_generation=true;
const bool run_relay=false;
const bool run_block_relay=false;
const bool walk_relay=false;
const bool noise_data=true;
char proc_loc[] = "./dat_dir/";
int main()
{
  assert(nbeads<=30);
  char rydat_dir[30]; sprintf(rydat_dir, "swirl%d.odr/", nbeads);
  char file_name[10]; sprintf(file_name, "pts.%d", param_id);
  char swbest_name[30]; sprintf(swbest_name, "swirl_best_re%d.odr/", relay_id);

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
  wl.add_wall(&wf); wl.add_wall(&wp0); wl.add_wall(&wp1); wl.add_wall(&wp2);

  double cl = sp_min.cl_im; // this will do for now, assuming that we aren't learning cl_im
  reporter rep; rep.init_relay(proc_loc, rydat_dir, file_name, relay_id, noise_data, noise_tol*cl);
  referee ref(nlead, npool, param_len, dt_sim, noise_tol, alpha_tol, rs_full_factor, cl);

  if (test_generation)
  {
    int test_id=2, Frames_test=1201;
    // int test_relay_id=(noise_data)?1:0;
    int test_relay_id=2;
    int npool_test = 1000;
    referee doc_ref(nlead, npool_test, param_len, dt_sim, noise_tol, alpha_tol, rs_full_factor, cl);
    printf("Testing %d bead generation. Test id: %d, Frames : %d\n", nbeads, test_id, Frames_test);
    doctor doc(doc_ref, sp_min,sp_max,wl,t_phys,&rep);
    doc.init_test(test_id, test_relay_id);
    doc.test_run(Frames_test);
  }
  if (run_relay)
  {
    printf("Running %d bead relay. Relay id: %d, parameter id: %d,  maximum generations: %d, leaders: %d, pool: %d \n", nbeads, relay_id, param_id, generations, nlead, npool);

    relay prelay(ref, sp_min,sp_max,wl,t_phys,&rep, t_wheels);
    prelay.init_relay();
    prelay.start_relay(generations);
  }
  else if (run_block_relay)
  {
    printf("Running %d bead block relay. Relay id: %d, parameter id: %d,  maximum generations: %d, leaders: %d, pool: %d \n", nbeads, relay_id, param_id, generations, nlead, npool);

    relay prelay(ref, sp_min,sp_max,wl,t_phys,&rep, t_wheels);
    prelay.init_relay();
    prelay.start_block_relay(generations);
  }
  else if (walk_relay)
  {
    printf("Walking %d bead relay. Relay id: %d, parameter id: %d,  maximum generations: %d, leaders: %d, pool: %d, training wheels: %e \n", nbeads, relay_id, param_id, generations, nlead, npool, t_wheels);

    relay prelay(ref, sp_min,sp_max,wl,t_phys,&rep);
    prelay.init_relay();
    prelay.t_wheels = t_wheels;
    prelay.start_relay(generations);
  }


  return 0;
}
