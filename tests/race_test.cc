#include <sys/types.h>
#include <sys/stat.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "particle_race.hh"

const int nbeads=3;
const int param_id=0;
const int generations=3;
char proc_loc[] = "./dat_dir/";
int main()
{
  assert(nbeads<=30);
  char rydat_dir[30]; sprintf(rydat_dir, "swirl%d.odr/", nbeads);
  char file_name[10]; sprintf(file_name, "pts.%d", param_id);

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
    rep.init_race(proc_loc, rydat_dir, file_name);
    //          nl, np, par_len, dt_sim, gau_var_h, gau_var_l, lambda_coeff, rs_fill_factor_, rs_full_factor
    referee ref(100, 1000, 12, 0.002, 1e-2, 1e-3, 1.0, 0.5, 0.8);
    race prace(ref,sp_min,sp_max,wl,t_phys,&rep);
    prace.init_race();
    prace.start_race(generations);
    prace.make_best_swirl("swirl_best.odr/");

    return 0;
}
