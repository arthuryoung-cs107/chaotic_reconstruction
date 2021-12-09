#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include "swirl.hh"

int main() {

    // Physical constants
    double g_phys=9.804,                 // Gravity (m/s^2)
           d_phys=0.00635,               // Diameter (m)
           t_phys=sqrt(d_phys/g_phys);   // Time unit (s)

    swirl_param sparam(0.5,1,1000,40,40,40,0.5,0.25,0.5,1.8,203,178,27.6,1.);

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
    // swirl sw(sparam,&pg,wl,25);
    swirl sw(sparam,&pg,wl,3);
    sw.import("input3.dat");

    // Solve the system
    sw.setup_output_dir("circ6.odr");
    dat_store dstore(t_phys,402,380,37.6);
    sw.dstore=&dstore;

    sw.solve(120,0.0005,1200);

    // dstore.write("synth25.dat");
    dstore.write("synth3.dat");
    dstore.print_theta_info();
}
