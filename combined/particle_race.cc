#include "particle_race.hh"

extern "C"
{
  #include "AYaux.h"
}

referee::~referee()
{
  if (alloc_flag)
  {
    delete leader_board;
    for (int i = 0; i < 2*nlead; i++) delete leaders[i];
    delete leaders;
    free_AYdmatrix(lead_params);
  }
}
void referee::alloc_records()
{
  leaders = new record*[nlead+npool];
  leader_board = new record*[nlead+npool];
  for (int i = 0; i < nlead+npool; i++) leaders[i] = new record(i);
  pool = leaders + nlead;
  pool_leaders = leader_board + nlead;

  lead_params = AYdmatrix(nlead+npool, param_len);
  pool_params = lead_params + nlead;

  alloc_flag = true;
}

runner::runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, double tol_): swirl(sp_, pg_, wl_, n_), thread_id(thread_id_), param_len(param_len_), Frames(Frames_), tol(tol_), x0(new double[2*n_])
{
  pvals = &Kn;
}
runner::~runner()
{delete x0;}

void runner::init_ics(double *x0_, double t0_, double ctheta0_)
// sets the initial conditions that we will return to for each simulation startup
{
  t0=t0_, ctheta0=ctheta0_; // t0 is given normalized by t_phys
  for (int i=0; i < 2*n; i++) x0[i] = x0_[i]; // read directly from training data
}

void runner::reset_sim(double *ptest_)
{
  // load the set of parameters we are to test
  for (int i = 0; i < param_len; i++) pvals[i] = ptest_[i];
  time=t0;
  set_swirl(ctheta0, 0.0);
  for (int i=0,j=0; i < n; i++,j+=2) q[i].set_pos(((x0[j]-cx_im)/cl_im)+cx, ((x0[j+1]-cy_im)/cl_im)+cy, rad);
}

void runner::run_race(double t_phys_, double dt_sim_, double *ts_, double *xs_, double *d_ang_)
{
  bool run_on=true;
  frame = 0; pos_err_acc=0.0;
  do
  {
    double dur=(ts_[frame+1]-ts_[frame])/t_phys_, ctheta=d_ang_[frame],comega=d_ang_[frame+1]-d_ang_[frame];
    if(comega>M_PI) comega-=2*M_PI;else if(comega<-M_PI) comega+=2*M_PI;
    comega/=dur;
    double *f = xs_ + (2*n*(frame+1));
    advance(dur, ctheta, comega, dt_sim_);
    if (is_lost(f)) run_on = false;
    else if (++frame == Frames-1) run_on = false;
  } while(run_on);
}

bool runner::is_lost(double *f_)
{
  double pos_err=0.0, max_err=0.0;
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double xt=(q[i].x-cx)*cl_im + cx_im - f_[j], yt=(q[i].y-cy)*cl_im + cy_im - f_[j+1], rsq = xt*xt+yt*yt, r = sqrt(rsq);
    if (r>max_err) max_err = r;
    pos_err+=r;
  }
  pos_err_acc+=pos_err;
  if (max_err > tol) return true; // lost
  else return false; // carry on then
}

int find_worst(record ** r, int ncap)
{
  int worst_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[worst_index]->isbetter(r[i]))
      worst_index = i;
  return worst_index;
}
