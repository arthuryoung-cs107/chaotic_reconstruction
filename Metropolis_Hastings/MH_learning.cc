#include "MH_learning.hh"

// record

record::record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record_struct(rs_), rid(rid_), ichunk(ichunk_), dchunk(dchunk_), u(u_) {}

// thread_worker

thread_worker::thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_): swirl(sp_, pg_, wl_, tws_.nbeads), thread_worker_struct(tws_),
thread_id(thread_id_),
u(&Kn), psim(new double[2*nbeads*Frames]) {}

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
sp_min(sp_min_), sp_max(sp_max_),
umin(&sp_min.Kn), umax(&sp_max.Kn),
wl(wl_), pg(new proximity_grid*[nt]), rng(new MH_rng*[nt])
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

// basic_record

basic_record::basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record(rs_, rid_, ichunk_, dchunk_, u_),
basic_rec_ints(&gen), basic_rec_dubs(&r2) {init_basic_record();}

basic_record::basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double * umin_, double * umax_): basic_record(rs_, rid_, ichunk_, dchunk_, u_)
{draw_ranuni(ran_,umin_,umax_);}

// basic_thread_worker

basic_thread_worker::basic_thread_worker(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_): thread_worker(sp_, pg_, wl_, tws_, thread_id_), event_detector(nbeads, Frames, 2, alpha_tol_),
basic_tw_ints(&nf_obs), basic_tw_dubs(&net_r2),
int_wkspc(new int[npool]), dub_wkspc(new double[ulen]) {}


// basic_MH_trainer

basic_MH_trainer::basic_MH_trainer(MH_train_struct &mhts_, int ichunk_width_, int dchunk_width_, double t_wheels0_): MH_trainer(mhts_, ichunk_width_, dchunk_width_), gaussian_likelihood(sigma, mhts_.sp_min->cl_im),
apply_training_wheels(t_wheels0_>0.0), t_wheels0(t_wheels0_),
ndup_leaders(new int[nlead]), irepl_leaders(new int[nlead]), isuccess_pool(new int[npool]),
w_leaders(new double[nlead]),
u_mean(new double[ulen]), u_var(new double[ulen]),
u_wmean(new double[ulen]), u_wvar(new double[ulen]) {}


void basic_MH_trainer::duplicate_u(basic_record *rec_child_, basic_record *rec_parent_, MH_rng *rng_t_)
{
  double  r_ = sqrt(rec_parent_->r2),
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
