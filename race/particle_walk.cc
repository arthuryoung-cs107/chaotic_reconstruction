#include "particle_walk.hh"

extern "C"
{
  #include "AYaux.h"
}

void grade::reset_grade(int gen_)
{
  frscore=0;
  gen=gen_;
  parent_gen=-1; // negative value implies no parent
  parent_count=0;
  parent_global_index=-1; // negative value implies no parent
  dup_count=0;

  l2score=0.0;
}
void grade::resample(int gen_, double * dmin_, double *dmax_, AYrng * r_)
{
  reset_grade(gen_);
  for (int i = 0; i < len; i++)
    params[i] = dmin_[i]+(dmax_[i]-dmin_[i])*(r_->rand_uni_gsl(0.0,1.0));
}
void grade::duplicate(grade *parent_, int gen_, double *dmin_, double *dmax_, AYrng * r_, double * var_)
{
  frscore=0;
  gen=gen_;
  parent_gen=parent_->gen;
  parent_count=parent_->parent_count+1;
  parent_global_index=parent_->global_index;
  dup_count=0;

  l2score=0;

  for (int i = 0; i < len; i++)
  {
    params[i] = parent_->params[i]*(r_->rand_gau_gsl(1.0, var_[i])); // do we account for the displacement relative to mean?
    if (params[i] > dmax_[i]) params[i]=dmax_[i];
    else if (params[i] < dmin_[i]) params[i]=dmin_[i];
  }
}
void grade::take_vals(grade * gtake_)
{
  frscore=gtake_->frscore;
  gen=gtake_->gen;
  parent_gen=gtake_->parent_gen;
  parent_count=gtake_->parent_count;
  parent_global_index=gtake_->parent_global_index;
  dup_count=gtake_->dup_count;

  l2score = gtake_->l2score;

  for (int i = 0; i < len; i++) params[i] = gtake_->params[i];
}





walker::walker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, double tol_, double t_phys_, double dt_sim_, double t_wheels_, double *ts_, double *xs_, double *d_ang_, int ic_index_, int nlead_, int npool_) : swirl(sp_, pg_, wl_, n_),
thread_id(thread_id_), param_len(param_len_), Frames(Frames_-ic_index_), tol(tol_),
t_phys(t_phys_),
dt_sim(dt_sim_), t_wheels(t_wheels_),
ts(ts_+ic_index_), xs(xs_+2*n_*ic_index_), d_ang(d_ang_+ic_index_),
nlead(nlead_), npool(npool_),
lead_dup_count(new int[nlead_]), frame_kill_count(new int[Frames_]),
frame_res_data(new double[4*Frames_]), param_mean(new double[param_len_])
{
  pvals = &Kn;
  t0_raw=*ts;
  t0=t0_raw/t_phys, x0=xs, ctheta0=*d_ang;
}
walker::~walker()
{
  delete [] lead_dup_count; delete [] frame_kill_count;
  delete [] frame_res_data; delete [] param_mean;
}
void walker::reset_sim(double *ptest_)
{
  // load the set of parameters we are to test
  for (int i = 0; i < param_len; i++) pvals[i] = ptest_[i];

  time=t0;
  set_swirl(ctheta0, 0.0);
  for (int i=0,j=0; i < n; i++,j+=2)
  {
    q[i].set_pos(((x0[j]-cx_im)/cl_im)+cx, ((x0[j+1]-cy_im)/cl_im)+cy, rad);
    q[i].zero_rest();
  }
}
int walker::start_walking(grade * gra_, int frscore_worst_, double l2score_worst_)
{
  reset_sim(gra_->params);
  int frame_kill=Frames;
  dead=false;
  pos_res_acc=0.0;
  for (frame = 0; frame < Frames-1; frame++)
  {
    double dur=(ts[frame+1]-ts[frame])/t_phys, ctheta=d_ang[frame],comega=d_ang[frame+1]-d_ang[frame];
    if(comega>M_PI) comega-=2*M_PI; else if(comega<-M_PI) comega+=2*M_PI;
    comega/=dur;
    double *f = xs + (2*n*(frame+1));
    advance(dur, ctheta, comega, dt_sim);
    compute_error(f,dur);
  }
  return (int)(gra_->check_success(frame_kill, pos_res_acc, frscore_worst_, l2score_worst_));
}
void walker::compute_error(double *f_, double dur_)
{
  double pos_res=0.0, max_err=0.0, fac=t_wheels/cl_im, idur=1.0/dur_;
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double xt=(q[i].x-cx)*cl_im + cx_im - f_[j], yt=(q[i].y-cy)*cl_im + cy_im - f_[j+1], rsq = xt*xt+yt*yt, r = sqrt(rsq);

    q[i].tweak_pos(-fac*xt,-fac*yt,idur); // training wheels

    pos_res+=rsq;
    if (r>max_err) max_err=r;
  }
  frame_res_data[frame*4]+=pos_res;
  frame_res_data[frame*4+1]+=max_err/=cl_im;
  pos_res_acc+=pos_res;

  if (dead) frame_res_data[frame*4+2]+= pos_res;
  else
  {
    frame_res_data[frame*4+3]+=pos_res;
    if (max_err>cl_im*cl_im) {frame_kill_count[frame]++; frame_kill=frame; dead=true;}
  }
}
void walker::reset_diagnostics()
{
  n_test=0;
  for (int i = 0; i < 4*Frames; i++) frame_res_data[i] = 0.0;
  for (int i = 0; i < Frames; i++) frame_kill_count[i] = 0;
}
void walker::consolidate_diagnostics()
{
  int n_alive=n_test,n_dead=0;
  for (int i = 0; i < Frames; i++)
  {
    frame_res_data[4*i] /= (double)n_test;
    frame_res_data[4*i+1] /= (double)n_test;

    n_dead+=frame_kill_count[i];
    if (n_dead) frame_res_data[4*i+2] /= (double)n_dead;
    frame_res_data[4*i+3] /= (double)n_alive;
    n_alive-=frame_kill_count[i];
  }
}
void walker::update_diagnostics(int * frame_kill_count_, double * gen_frame_res_data_)
{
  double w = ((double)n_test)/((double)npool);
  for (int i = 0; i < Frames; i++)
  {
    frame_kill_count_[i]+=frame_kill_count[i];
    gen_frame_res_data_[4*i] += w*frame_res_data[4*i];
    gen_frame_res_data_[4*i+1] += w*frame_res_data[4*i+1];
    gen_frame_res_data_[4*i+2] += w*frame_res_data[4*i+2];
    gen_frame_res_data_[4*i+3] += w*frame_res_data[4*i+3];
  }
}






guide::~guide()
{
  if (alloc_flag)
  {
    free_AYdmatrix(param_chunk);

    delete [] leader_board;
    for (int i = 0; i < nlead+npool; i++) delete grades[i];
    delete [] grades;

    delete [] lead_dup_count; delete [] frame_kill_count;

    delete [] sample_weights;
    delete [] gen_frame_res_data;
    delete [] gen_param_mean;
    delete [] gen_param_var;
  }
}
void guide::alloc_grades(int nt_, int Frames_)
{
  param_chunk = AYdmatrix(nlead+npool, param_len);

  grades = leaders = classA = new grade*[nlead+npool];
  leader_board = new grade*[nlead+npool];

  for (int i = 0; i < nlead+npool; i++) grades[i] = new grade(i,param_chunk[i]);
  pool = grades + nlead;
  classB = grades + nA;
  candidates = leader_board + nlead;

  lead_dup_count = new int[nlead];
  frame_kill_count = new int [Frames_];

  sample_weights = new double[nlead];
  gen_frame_res_data = new double[4*Frames_];
  gen_param_mean = new double[param_len];
  gen_param_var = new double[param_len];

  alloc_flag = true;
}

int find_worst_grade(grade ** r, int ncap)
{
  int worst_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[worst_index]->isbetter(r[i]))
      worst_index = i;
  return worst_index;
}

int find_best_grade(grade ** r, int ncap)
{
  int best_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[best_index]->isworse(r[i]))
      best_index = i;
  return best_index;
}
