#include "particle_relay.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

extern "C"
{
  #include "AYaux.h"
}


void relay::compute_leader_statistics()
{
  lambda=1.0/((double)leader_count);
  beta=(root_res_best/tau)*(0.5*(root_res_best/tau));

  w_max = w_sum = 0.0;
#pragma omp parallel for reduction(+:w_sum) reduction(max:w_max)
  for (int i = 0; i < leader_count; i++)
  {
    w_sum += sample_weights[i] = leaders[i]->w(lambda,beta,tau);
    if (sample_weights[i]>w_max) w_max=sample_weights[i];
  }

  for (int i = 0; i < param_len; i++) lead_par_w_mean[i] = lead_par_w_var[i] = 0.0;

#pragma omp parallel
  {
    runner *rt = runners[thread_num()];
    double * param_buf = rt->param_acc;
    for (int i = 0; i < param_len; i++) param_buf[i] = 0.0;

#pragma omp for nowait
    for (int i = 0; i < leader_count; i++)
    {
      for (int j = 0; j < param_len; j++) param_buf[j]+=(sample_weights[i]*(param_acc_factor[j]/w_max))*(leaders[i]->params[j]);
    }
    for (int j = 0; j < param_len; j++) param_buf[j]*=(w_max/w_sum)/param_acc_factor[j];
#pragma omp critical
    {
      for (int i = 0; i < param_len; i++) lead_par_w_mean[i] += param_buf[i];
    }
    for (int i = 0; i < param_len; i++) param_buf[i] = 0.0;

#pragma omp barrier

#pragma omp for nowait
    for (int i = 0; i < leader_count; i++)
    {
      double * params_i = leaders[i]->params;
      for (int j = 0; j < param_len; j++)
      {
        double z = (params_i[j] - lead_par_w_mean[j]);
        param_buf[j]+=(sample_weights[i]*(param_acc_factor[j]/w_max))*(z*z);
      }
    }
    for (int j = 0; j < param_len; j++) param_buf[j]*=(w_max/(w_sum-1.0))/param_acc_factor[j];
#pragma omp critical
    {
      for (int i = 0; i < param_len; i++) lead_par_w_var[i] += param_buf[i];
    }
  }
}

void relay::resample_pool()
{
  compute_leader_statistics();

  for (int i = 0; i < leader_count; i++) lead_dup_count[i] = 0;

  dup_count=resample_count=dup_unique=0;
#pragma omp parallel
    {
      int t = thread_num();
      AYrng *r = rng[t];
      int *dup_t = runners[t]->lead_dup_count;
      for (int i = 0; i < leader_count; i++) dup_t[i] = 0;
#pragma omp for reduction(+:dup_count) reduction(+:resample_count) nowait
      for (int i = 0; i < npool; i++)
      {
        int j = 0;
        double uni = (w_sum/rs_full_factor)*(r->rand_uni_gsl(0.0, 1.0));
        while ((j<leader_count)&&(uni>0.0)) uni -= sample_weights[j++];
        double *dmin=&sp_min.Kn, *dmax=&sp_max.Kn;
        if (j > 0) // if we have particles worth resampling
        {
          // resample the particle, add gaussian noise
          if (uni<0.0)
          {dup_count++; dup_t[--j]++; pool[i]->duplicate(leaders[j], gen_count, dmin, dmax,r, lead_par_w_var);}
          // we hit the resampling pool
          else
            {resample_count++; pool[i]->resample(gen_count, dmin, dmax, r);}
        }
        // we currently have no leaders (particles worth resampling)
        else
          {resample_count++; pool[i]->resample(gen_count, dmin, dmax, r);}
      }
#pragma omp critical
      {
        for (int i = 0; i < leader_count; i++) lead_dup_count[i]+=dup_t[i];
      }
    }

    for (int i = 0; i < leader_count; i++) if (lead_dup_count[i])
    {leaders[i]->dup_count+=lead_dup_count[i]; dup_unique++;}

    printf("%d duplicates (%d unique)\n", dup_count, dup_unique);
}
