#include <sys/types.h>
#include <sys/stat.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "MH_solvers.hh"

const int param_id=0;
const int MH_id=0;

const int ulen = 12;
const int nbeads=3;
const int nlead = 12;
const int npool = 1000;
const double dt_sim = 0.002;
const double t_wheels = 0.012;

const double noise_sigma = 1e-1;
const double alpha_tol=10.0;
const double rs_full_factor=1.0;

const int generations=100;

const bool test_generation=true;
const bool learn_data=false;
const bool noise_data=true;

char proc_loc[] = "./dat_dir/";

int main()
{
  assert(nbeads<=30);
  char test_dir[30]; sprintf(test_dir, "swirl%d.odr/", nbeads);
  char data_name[10]; sprintf(data_name, "pts.%d", param_id);

  // Physical constants
  double g_phys=9.804,                 // Gravity (m/s^2)
         d_phys=0.00635,               // Diameter (m)
         t_phys=sqrt(d_phys/g_phys);   // Time unit (s)

  double  sp_min_vals[14], sp_max_vals[14];
  int idmin=set_special_u("min", sp_min_vals),
  idmax=set_special_u("max", sp_max_vals);

  swirl_param sp_min(sp_min_vals), sp_max(sp_max_vals);
  double noise_scale=sp_min.cl_im;

  // Create the hexagonal dish
  wall_list wl;
  const double r=5.72;
  const double fa=sqrt(0.75);
  wall_floor wf(0.);
  wall_par_planes wp0(0,1,0,r),wp1(fa,0.5,0,r),wp2(fa,-0.5,0,r);
  wl.add_wall(&wf); wl.add_wall(&wp0); wl.add_wall(&wp1); wl.add_wall(&wp2);

  // test parameters
  int test_id=2,
      Frames_test=1201,
      test_relay_id=2,
      verbose=true;

  // sampling parameters
  MH_io mh_io(proc_loc, test_dir, data_name, MH_id, noise_data, noise_sigma*noise_scale);
  MH_params mh_par(&mh_io, ulen, nlead, npool, dt_sim, t_phys, noise_sigma);
  MH_train_struct mhts(&mh_par, &sp_min, &sp_max, &wl);

  MH_doctor mh_doc(mhts, test_id, test_relay_id, Frames_test, alpha_tol);
  mh_doc.run(verbose);

  return 0;
}
