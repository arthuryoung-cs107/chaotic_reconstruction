#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "swirl.hh"

int main() {

    // Physical constants
    double g_phys=9.804,                 // Gravity (m/s^2)
           d_phys=0.00635,               // Diameter (m)
           t_phys=sqrt(d_phys/g_phys);   // Time unit (s)

    // double sptrue_vals[] = {0.5,1.0,1000.0,40.0 ,40.0 ,40.0 ,0.5,0.25,0.5,1.8,203.0,178.0,27.6,1.0};
    double sptrue_vals[] = {0.5,1.0,900.0,35.0 ,45.0 ,40.0 ,0.4,0.20,0.4,1.8,203.0,178.0,27.6,1.0};

    swirl_param sparam(sptrue_vals);

    // Create the hexagonal dish
    wall_list wl;
    const double r=5.72,fa=sqrt(0.75);
    wall_floor wf(0.);
    wall_par_planes wp0(0,1,0,r),wp1(fa,0.5,0,r),wp2(fa,-0.5,0,r);
    wl.add_wall(&wf);
    wl.add_wall(&wp0);
    wl.add_wall(&wp1);
    wl.add_wall(&wp2);

    // Set the initial positions of the splines
    proximity_grid pg;
    swirl sw(sparam,&pg,wl,3);
    // swirl sw(sparam,&pg,wl,10);
    // swirl sw(sparam,&pg,wl,30);

    // sw.import("input_dir/input3.dat");
    // sw.import("input_dir/input10.dat");
    // sw.import("input_dir/input30.dat");
    sw.import_true("input_dir/input3_race.dat");

    // Solve the system
    ODR_struct odr;
    odr.init_swirl("./dat_dir/", "race_3beads.odr/", "pts");
    // ODR_struct odr("./dat_dir/", "race_10beads.odr/", "pts");
    // ODR_struct odr("./dat_dir/", "race2_30beads.odr/", "pts");

    odr.set_vidspecs(t_phys, sptrue_vals[10], sptrue_vals[11], sptrue_vals[12]);
    sw.setup_output_dir(&odr);

    sw.solve(120,0.0005,1200);

    odr.end_writing();
    odr.print_time_rotation();
}
