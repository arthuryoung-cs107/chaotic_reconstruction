#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>
#include <omp.h>

#include "swirl.hh"

int main()
{
    int par_len = 14, noise_len = 100+1;
    double sptrue[] = {0.5,1.0,1000.0,40.0 ,40.0 ,40.0 ,0.5,0.25,0.5,1.8,203.0,178.0,27.6,1.0};
    double sp_min[] = {0.5,1.0,500.0 ,5.0  ,5.0  ,5.0  ,0.1,0.1 ,0.1,1.8,203.0,178.0,27.6,1.0};
    double sp_max[] = {0.5,1.0,5000.0,120.0,120.0,120.0,1.0,1.0 ,1.0,1.8,203.0,178.0,27.6,1.0};
    double * sp_del = new double[par_len];
    for (int i = 0; i < par_len; i++)
    {
      sp_del[i] = sp_max[i]-sptrue[i];
      if ((sptrue[i]-sp_min[i]) < sp_del[i]) sp_del[i] = sptrue[i]-sp_min[i];
    }

    AYmat par_mat(par_len, noise_len); memcpy(par_mat.AT[0], sptrue, (size_t)(par_mat.M*sizeof(double)));
    normal gen(-1.0, 1.0);
    for (int i = 1; i < noise_len; i++) for (int j = 0; j < par_len; j++) par_mat.AT[i][j] = sptrue[j] + (sp_del[j]*gen.rand_gen());

    delete sp_del;
    #pragma omp parallel for
      for (int i = 0; i < noise_len; i++)
      {
        char proc_loc[] = "./dat_dir/";
        char rydat_dir[] = "stat3.odr/";
        char file_name[10]; sprintf(file_name, "rand.%d", i);
        // Physical constants
        double g_phys=9.804,                 // Gravity (m/s^2)
               d_phys=0.00635,               // Diameter (m)
               t_phys=sqrt(d_phys/g_phys);   // Time unit (s)


        swirl_param sparam(par_mat.AT[i]);
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
        sw.import("input_dir/input3.dat");

        // Solve the system
        ODR_struct odr(proc_loc, rydat_dir, file_name, false);
        odr.set_vidspecs(t_phys);
        sw.setup_output_dir(&odr, false, false);

        sw.solve(120,0.0005,1200);
        odr.end_writing();
        printf("(thread %d of %d): completed %s\n", omp_get_thread_num(), omp_get_max_threads(), file_name);
      }

}
