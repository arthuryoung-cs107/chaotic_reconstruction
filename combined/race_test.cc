#include <sys/types.h>
#include <sys/stat.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "particle_race.hh"

int main() {
    // Physical constants
    double g_phys=9.804,                 // Gravity (m/s^2)
           d_phys=0.00635,               // Diameter (m)
           t_phys=sqrt(d_phys/g_phys);   // Time unit (s)

    double  sp_min_vals[] = {0.5,1.0,500.0 ,5.0  ,5.0  ,5.0  ,0.1,0.1 ,0.1,1.8,203.0,178.0,27.6,1.0},
            sptrue_vals[] = {0.5,1.0,1000.0,40.0 ,40.0 ,40.0 ,0.5,0.25,0.5,1.8,203.0,178.0,27.6,1.0},
            sp_max_vals[] = {0.5,1.0,5000.0,120.0,120.0,120.0,1.0,1.0 ,1.0,1.8,203.0,178.0,27.6,1.0};


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

    ODR_struct odr("./dat_dir/race_3beads.odr/pts");
    referee ref;
    race prace(ref,sp_min,sp_max,wl,t_phys,odr);

    // prace.setup_output_info(255,200);

    // Solve the system
    // prace.init_race(65536);
    // prace.run(1200);

    return 0;
}
