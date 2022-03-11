#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>
#include <assert.h>
#include <omp.h>

#include "swirl.hh"

char test[]="gauss";
const int nbeads=3;
const int ran_id=0;
char proc_loc[] = "./dat_dir/";
int main()
{
    assert(nbeads<=30);
    char rydat_dir[30]; sprintf(rydat_dir, "stat%d_%s%d.odr/", nbeads, test, ran_id);
    int par_len = 14, noise_len = 100+1;
    double sp_min[14], sptrue[14], sp_max[14];
    int idtrue=set_special_params("true", sptrue),
        idmin=set_special_params("min", sp_min),
        idmax=set_special_params("max", sp_max);

    AYvec sp_del(par_len);
    AYmat par_mat(par_len, noise_len); memcpy(par_mat.AT[0], sptrue, (size_t)(par_mat.M*sizeof(double)));

    for (int i = 0; i < par_len; i++)
      sp_del.A_ptr[i] = sp_max[i]-sp_min[i];

    AYrng gen; gen.rng_init_gsl();

    if (strcmp(test, "maxmin")==0)
      for (int i = 1; i < noise_len; i++) for (int j = 0; j < par_len; j++)
        par_mat.AT[i][j] = sp_min[j] + (sp_del.A_ptr[j]*gen.rand_uni_gsl(0.0, 1.0));
    else if (strcmp(test, "gauss")==0)
      for (int i = 1; i < noise_len; i++) for (int j = 0; j < par_len; j++)
      {par_mat.AT[i][j] = sptrue[j]*gen.rand_gau_gsl(1.0,1e-8);
        if (par_mat.AT[i][j] > sp_max[j]) par_mat.AT[i][j] = sp_max[j];
        else if (par_mat.AT[i][j] < sp_min[j]) par_mat.AT[i][j] = sp_min[j];}
    else
      for (int i = 1; i < noise_len; i++) for (int j = 0; j < par_len; j++)
        par_mat.AT[i][j] = sp_min[j] + (sp_del.A_ptr[j]*gen.rand_uni_gsl(0.0, 1.0));

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
        swirl sw(sparam,&pg,wl,nbeads);
        sw.import("input_dir/input30.dat");

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
