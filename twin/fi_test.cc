#include <sys/types.h>
#include <sys/stat.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "filter.hh"

int main() {

    // Physical constants
    double g_phys=9.804,                 // Gravity (m/s^2)
           d_phys=0.00635,               // Diameter (m)
           t_phys=sqrt(d_phys/g_phys);   // Time unit (s)

    // Filter parameters
    fil_param fparam(0.002,0.1,1,0.01,0.01,0.012,3);

    // Minimum and maximum parameters
    swirl_param sp_min(0.5,1,500,5,5,5,0.1,0.1,0.1,1.8,402,380,37.6,1),
                sp_max(0.5,1,5000,120,120,120,1.,1.,1.,1.8,402,380,37.6,1),
                sp_rnd(sp_min,sp_max,0.005);

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

    filter fi(fparam,sp_min,sp_max,sp_rnd,wl,t_phys,"synth25.dat",0);
    //filter fi(fparam,sp_min,sp_max,sp_rnd,wl,t_phys,"data/50_ramp.dat",2500);
    //filter fi(fparam,sp_min,sp_max,sp_rnd,wl,t_phys,"data/100_ramp.dat",2450);
    fi.setup_output_info(255,"sy2.odr",200);

    // Solve the system
    fi.init(65536);
    fi.run(1200);
}
