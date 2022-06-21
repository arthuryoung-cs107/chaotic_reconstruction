#include "MH_learning.hh"

// thread_worker

void thread_worker::reset_sim(double *utest_, double t0_, double ctheta0_, double comega0_, double *p0_)
{
  for (int i = 0; i < ulen; i++) u[i]=utest_[i];

  time=t0_;
  set_swirl(ctheta0_, comega0_);
  for (int i = 0, j = 0; i < nbeads; i++,j+=2)
  {
    double  x_it=((p0_[j]-cx_im)/cl_im)+cx,
            y_it=((p0_[j+1]-cy_im)/cl_im)+cy;
    psim[j]=x_it; psim[j+1]=y_it;
    q[i].set_pos(x_it,y_it,rad);
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

int MH_trainer::find_worst_record(record ** r_, int ncap_)
{
  int worst_index = 0;
  for (int i = 1; i < ncap_; i++)
    if (r_[worst_index]->isbetter(r_[i]))
      worst_index = i;
  return worst_index;
}

int MH_trainer::find_best_record(record ** r_, int ncap_)
{
  int best_index = 0;
  for (int i = 1; i < ncap_; i++)
    if (r_[best_index]->isworse(r_[i]))
      best_index = i;
  return best_index;
}

void MH_trainer::pick_nworst_records(record ** rin_, record ** rout_, int n_, int ncap_)
{
  for (int i = 0; i < n_; i++) rout_[i]=rin_[i];
  int i_best_worst = find_best_record(rout_,n_);
  for (int i = n_; i < ncap_; i++)
    if (rout_[i_best_worst]->isbetter(rin_[i]))
    {
      rout_[i_best_worst] = rin_[i];
      i_best_worst=find_best_record(rout_,n_);
    }
}

void MH_trainer::pick_nbest_records(record ** rin_, record ** rout_, int n_, int ncap_)
{
  for (int i = 0; i < n_; i++) rout_[i]=rin_[i];
  int i_worst_best = find_worst_record(rout_,n_);
  for (int i = n_; i < ncap_; i++)
    if (rout_[i_worst_best]->isworse(rin_[i]))
    {
      rout_[i_worst_best] = rin_[i];
      i_worst_best=find_worst_record(rout_,n_);
    }
}

void MH_trainer::pick_nworst_records(record ** rin_, int n_, int ncap_)
{
  int i_best_worst = find_best_record(rin_,n_);
  for (int i = n_; i < ncap_; i++)
    if (rin_[i_best_worst]->isbetter(rin_[i]))
    {
      rin_[i_best_worst] = rin_[i];
      i_best_worst=find_best_record(rin_,n_);
    }
}

void MH_trainer::pick_nbest_records(record ** rin_, int n_, int ncap_)
{
  int i_worst_best = find_worst_record(rin_,n_);
  for (int i = n_; i < ncap_; i++)
    if (rin_[i_worst_best]->isworse(rin_[i]))
    {
      rin_[i_worst_best] = rin_[i];
      i_worst_best=find_worst_record(rin_,n_);
    }
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

void basic_MH_trainer::duplicate_u(basic_record *rec_child_, basic_record *rec_parent_, MH_rng *rng_t_)
{
  double  r_ = sqrt(rec_pool_->r2compare),
          r_rat = r_/sqrt(rho2),
          sigma_fac = max(0.25*(1.0-exp(0.5*(1.0-r_rat)*(1.0+r_rat))), 0.0),
          *u_child=rec_child_->u,
          *u_parent=rec_parent_->u;
  for (int i = 0; i < ulen; i++)
  {
    double z = rng_t_->rand_gau();
    u_child[i]=u_parent[i] + z*sigma_fac*((z>0.0)?(umax[i]-u_parent[i]):(u_parent[i]-umin[i]));
    if (u_child[i]>umax[i]) u_child[i]=umax[i];
    if (u_child[i]<umin[i]) u_child[i]=umin[i];
  }
  rec_child_->init_basic_record(gen_count);
}
