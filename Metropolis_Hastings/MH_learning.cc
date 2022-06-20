#include "MH_learning.hh"

// thread_worker

void thread_worker::reset_sim(double *utest_, double t0_, double ctheta0_, double comega0_, double *p0_)
{
  for (int i = 0; i < ulen; i++) u[i]=utest_[i];

  time=t0_;
  set_swirl(ctheta0_, comega0_);
  for (int i = 0, j = 0; i < nbeads; i++,j+=2)
  {
    q[i].set_pos(((p0_[j]-cx_im)/cl_im)+cx, ((p0_[j+1]-cy_im)/cl_im)+cy, rad);
    q[i].zero_rest();
  }
}


// MH_trainer

MH_trainer::MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_, int ichunk_width_, int dchunk_width_) : MH_params(par_),
MHT_it_ints(&leader_count), MHT_it_dubs(&rho2),
nt(get_nt()), ichunk_width(ichunk_width_), dchunk_width(dchunk_width_),
ichunk(Tmatrix<int>(nlead+npool, ichunk_width)),
ts(new double[Frames]), xs(new double[2*nbeads*Frames]), d_ang(new double[Frames]), comega_s(new double[Frames]),
dchunk(Tmatrix<double>(nlead+npool, dchunk_width)), uchunk(Tmatrix<double>(nlead+npool, ulen)),
sp_min(sp_min_), sp_max(sp_max_), wl(wl_),
pg(new proximity_grid*[nt]), rng(new MH_rng*[nt])
{
  io->load_reference(ts, xs, d_ang, comega_s, t_phys);

  // Set up each runner's personal data
#pragma omp parallel
  {
    int t=thread_num();
    pg[t]=new proximity_grid();
    rng[t]=new MH_rng(t+1);
  }
}

MH_trainer::~MH_trainer()
{
  delete [] ts;
  delete [] xs;
  delete [] d_ang;
  delete [] comega_s;
  free_Tmatrix<int>(ichunk);
  free_Tmatrix<double>(dchunk);
  free_Tmatrix<double>(uchunk);

  for (int i = 0; i < nt; i++)
  {
    delete pg[i];
    delete rng[i];
  }
  delete [] pg;
  delete [] rng;
}

double basic_MH_trainer::compute_weights(double r2_min_, basic_record ** recs_, int n_)
{
  double  wsum=0.0,
          r_min_=sqrt(r2_min_),
          rho_=sqrt(rho2);

  #pragma omp parallel for reduction(+:wsum)
  for (int i = 0; i < n_; i++)
  {
    wsum+=w_leaders[i]=recs_[i]->w=gaussian_likelihood::compute_weight(sqrt(recs_[i]->r2compare),r_min_,rho_);
  }
  return wsum;
}


void basic_MH_trainer::respawn_pool(double w_sum_, basic_thread_worker ** tws_, basic_record ** pool_, basic_record ** leaders_)
{
  memset(ndup_leaders,0,nlead*sizeof(int));
  ndup=nredraw=ndup_unique=0;

  #pragma omp_parallel
  {
    int tid = thread_num(),
        *dup_t = tws_[tid]->int_wkspc;
    MH_rng * rng_t = rng[tid];
    memset(dup_t,0,nlead*sizeof(int));
    #pragma omp for reduction(+:ndup) reduction(+:nredraw) nowait
    for (int i = 0; i < npool; i++)
    {
      int j=0;
      double uni = (w_sum_/rs_full_factor)*rng_t->rand_uni();
      while ((j<leader_count)&&(uni>0.0)) uni-=w_leaders[j++];
      if (j>0)
      {
        if (uni<0.0)
        {ndup++; dup_t[--j]++; duplicate_u(pool_[i],leaders_[j],rng_t);}
        else
        {nredraw++; redraw_u(pool_[i],rng_t);}
      }
      else
      {nredraw++; redraw_u(pool_[i],rng_t);}
    }
    #pragma omp critical
    {
      for (int i = 0; i < leader_count; i++) ndup_leaders[i]+=dup_t[i];
    }
  }

  for (int i = 0; i < leader_count; i++) if (ndup_leaders[i])
  {leaders[i]->dup_count+=ndup_leaders[i]; ndup_unique++;}
}

void basic_MH_trainer::basic_duplicate_u(basic_record *rec_child_, basic_record *rec_parent_, MH_rng *rng_t_)
{
  double  r_ = sqrt(rec_pool_->r2compare),
          r_rat = r_/sqrt(rho2),
          sigma_fac = max(0.25*(1.0-exp(0.5*(1.0-r_rat)*(1.0+r_rat))), 0.0),
          *u_child=rec_child_->u,
          *u_parent=rec_parent_->u;
}
