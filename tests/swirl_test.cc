#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>
#include <assert.h>

#include "swirl.hh"

const int nbeads=10;
const int param_id=0;
const bool write_split = false;
char proc_loc[] = "./dat_dir/";
int main()
{
    assert(nbeads<=30);
    double sp_vals[14];
    int param_use = set_special_params(param_id, sp_vals);
    char rydat_dir[30]; sprintf(rydat_dir, "swirl%d.odr/", nbeads);
    char file_name[10]; sprintf(file_name, "pts.%d", param_use);

    // Physical constants
    double g_phys=9.804,                 // Gravity (m/s^2)
           d_phys=0.00635,               // Diameter (m)
           t_phys=sqrt(d_phys/g_phys);   // Time unit (s)


    swirl_param sparam(sp_vals);

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
    swirl sw(sparam,&pg,wl,nbeads);

    sw.import("input_dir/input30.dat");

    // Solve the system
    ODR_struct odr;
    odr.init_swirl(proc_loc, rydat_dir, file_name, write_split);

    odr.set_vidspecs(t_phys, sparam.cx_im, sparam.cy_im, sparam.cl_im);
    sw.setup_output_dir(&odr, true, true);

    sw.solve(120,0.0005,1200);

    odr.end_writing();
}
