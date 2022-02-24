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
             //       (rad, mass, Kn  , gnb, gnf, gnw, mub, muf , muw, amp, cx , cy , cl  , sca)
             // sp_min(0.5, 1   , 500 , 5  , 5  , 5  , 0.1, 0.1 , 0.1, 1.8, 402, 380, 37.6, 1  ) MIN
             // sp_max(0.5, 1   , 5000, 120, 120, 120, 1.0, 1.0 , 1.0, 1.8, 402, 380, 37.6, 1  ) MAX
             // sptrue(0.5, 1   , 1000, 40 , 40 , 40 , 0.5, 0.25, 0.5, 1.8, 203, 178, 27.6, 1.0) TRUE
             //       ( c , c   , *   , *  , *  , *  , *  , *   , *  , c  , ?  , ?  , ?   , c  ) KEY

    //original
    swirl_param sp_min(0.5,1,500,5,5,5,0.1,0.1,0.1,1.8,402,380,37.6,1),
                sp_max(0.5,1,5000,120,120,120,1.,1.,1.,1.8,402,380,37.6,1),
                sp_rnd(sp_min,sp_max,0.005);
    //true
    swirl_param sptrue(0.5,1,1000,40,40,40,0.5,0.25,0.5,1.8,203,178,27.6,1.);
    //corrected?
    // swirl_param sp_min(0.5,1,500,5,5,5,0.1,0.1,0.1,1.8,203.0,178.0,27.6,1.0),
    //             sp_max(0.5,1,5000,120,120,120,1.,1.,1.,1.8,203.0,178.0,27.6,1.0),
    //             sp_rnd(sp_min,sp_max,0.005);

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

    ODR_struct odr("./dat_dir/circ6_swrl.odr/pts");
    filter fi(fparam,sp_min,sp_max,sp_rnd,wl,t_phys,&odr,0);

    fi.setup_output_info(255,200);

    // Solve the system
    fi.init(65536);
    fi.run(1200);

    return 0;
}
