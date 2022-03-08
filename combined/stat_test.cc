#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>
#include <omp.h>

#include "swirl.hh"

const int test = 2;

int main()
{
    int par_len = 14, noise_len = 100+1;
    // indices         0   1   2      3     4     5     6   7    8   9   10    11    12   13
    double sp_min[] = {0.5,1.0,500.0 ,5.0  ,5.0  ,5.0  ,0.1,0.1 ,0.1,1.8,203.0,178.0,27.6,1.0};
    double sptrue[] = {0.5,1.0,1000.0,40.0 ,40.0 ,40.0 ,0.5,0.25,0.5,1.8,203.0,178.0,27.6,1.0};
    double sp_max[] = {0.5,1.0,5000.0,120.0,120.0,120.0,1.0,1.0 ,1.0,1.8,203.0,178.0,27.6,1.0};
    // double sp_min[] = {0.5,1.0,500.0 ,5.0  ,5.0  ,5.0  ,0.1,0.1 ,0.1,1.7,203.0,178.0,27.6,1.0};
    // double sptrue[] = {0.5,1.0,1000.0,40.0 ,40.0 ,40.0 ,0.5,0.25,0.5,1.8,203.0,178.0,27.6,1.0};
    // double sp_max[] = {0.5,1.0,5000.0,120.0,120.0,120.0,1.0,1.0 ,1.0,1.9,203.0,178.0,27.6,1.0};
    AYvec sp_del(par_len);
    AYmat par_mat(par_len, noise_len); memcpy(par_mat.AT[0], sptrue, (size_t)(par_mat.M*sizeof(double)));

    char rydat_dir[20];
    switch(test)
    {
      case 1:
      {
        for (int i = 0; i < par_len; i++)
        {
          sp_del.A_ptr[i] = sp_max[i]-sptrue[i];
          if ((sptrue[i]-sp_min[i]) < sp_del.A_ptr[i]) sp_del.A_ptr[i] = sptrue[i]-sp_min[i];
        }
        strcpy(rydat_dir, "stat3_maxmin.odr/");
        AYuniform uni(-1.0, 1.0);
        for (int i = 1; i < noise_len; i++) for (int j = 0; j < par_len; j++) par_mat.AT[i][j] = sptrue[j] + (sp_del.A_ptr[j]*uni.rand_gen());
      }
      break;
      case 2:
      {
        sp_del.init_0();
        for (int i = 2; i < 9; i++) sp_del.A_ptr[i] = fabs(sptrue[i]);
        strcpy(rydat_dir, "stat3_gauss.odr/");
        AYnormal norm(0.0, 1e-5);
        for (int i = 1; i < noise_len; i++) for (int j = 0; j < par_len; j++) par_mat.AT[i][j] = sptrue[j] + (sp_del.A_ptr[j]*norm.rand_gen());
      }
      break;
      default:
      {
        for (int i = 0; i < par_len; i++)
        {
          sp_del.A_ptr[i] = sp_max[i]-sptrue[i];
          if ((sptrue[i]-sp_min[i]) < sp_del.A_ptr[i]) sp_del.A_ptr[i] = sptrue[i]-sp_min[i];
        }
        strcpy(rydat_dir, "stat3_maxmin.odr/");
        AYuniform uni(-1.0, 1.0);
        for (int i = 1; i < noise_len; i++) for (int j = 0; j < par_len; j++) par_mat.AT[i][j] = sptrue[j] + (sp_del.A_ptr[j]*uni.rand_gen());
      }
    }

    char proc_loc[] = "./dat_dir/";
#pragma omp parallel for schedule(dynamic)
      for (int i = 0; i < noise_len; i++)
      {
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
        sw.import_true("input_dir/input3_race.dat");

        // Solve the system
        ODR_struct odr;
        odr.init_swirl(proc_loc, rydat_dir, file_name, false);
        odr.set_vidspecs(t_phys, sparam.cx_im, sparam.cy_im, sparam.cl_im);
        sw.setup_output_dir(&odr, false, false);

        sw.solve(120,0.0005,1200);
        odr.end_writing();
        printf("(thread %d of %d): completed %s\n", omp_get_thread_num(), omp_get_max_threads(), file_name);
      }

    char specs_name[100];
    sprintf(specs_name, "%s%srand_pars", proc_loc, rydat_dir);
    par_mat.fprintf_mat(specs_name);
    sprintf(specs_name, "%s%srand_dels", proc_loc, rydat_dir);
    sp_del.fprintf_mat(specs_name);

}
