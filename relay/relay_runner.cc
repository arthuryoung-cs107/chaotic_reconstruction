#include "particle_relay.hh"

extern "C"
{
  #include "AYaux.h"
}
runner::runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, int nlead_, int npool_, double t_phys_, double dt_sim_, double alpha_tol_, double *ts_, double *xs_, double *d_ang_, double *comega_s_) : swirl(sp_, pg_, wl_, n_),
thread_id(thread_id_), param_len(param_len_), Frames(Frames_), nlead(nlead_), npool(npool_),
t_phys(t_phys_), dt_sim(dt_sim_), alpha_tol(alpha_tol_),
ts(ts_), xs(xs_), d_ang(d_ang_), comega_s(comega_s_),
lead_dup_count(new int[nlead_]), kill_frames(new int[n_]), event_frame_count(AYimatrix(n_, Frames_)),
param_acc(new double[param_len_]), pos_res(AYdmatrix(Frames_, n_)), alpha_INTpos_res(AYdmatrix(Frames_, n_)), INTpos_res(AYdmatrix(n_, 3))
{
  pvals = &Kn;
  t0_raw=*ts;
  t0=t0_raw/t_phys, x0=xs, ctheta0=*d_ang;
}
runner::~runner()
{
  delete [] lead_dup_count; delete [] kill_frames; delete [] param_acc;
  free_AYimatrix(event_frame_count); free_AYdmatrix(pos_res);
  free_AYdmatrix(alpha_INTpos_res); free_AYdmatrix(INTpos_res);
}
void runner::reset_sim(double *ptest_, double t0_, double ctheta0_, double comega0_, double *x0_)
{
  // load the set of parameters we are to test
  for (int i = 0; i < param_len; i++) pvals[i] = ptest_[i];

  time=t0_;
  set_swirl(ctheta0_, comega0_);
  for (int i=0,j=0; i < n; i++,j+=2)
  {
    q[i].set_pos(((x0_[j]-cx_im)/cl_im)+cx, ((x0_[j+1]-cy_im)/cl_im)+cy, rad);
    q[i].zero_rest();
  }
}

// take two steps so that we can start off the event diagnostics collection
int runner::start_detection(int start_)
{
  int frame_local = start_+1;
  double dur, *f = xs+(2*n*frame_local), res_acc_local=0.0;
  advance(dur=(ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double xt=(q[i].x-cx)*cl_im+cx_im-f[j], yt=(q[i].y-cy)*cl_im+cy_im-f[j+1], rsq=xt*xt+yt*yt;
    pos_res[start_][i]=INTpos_res[start_][i]=0.0; // zero out the first frame
    res_acc_local += pos_res[frame_local][i]=rsq;
    INTpos_res[i][2] = 0.5*rsq*dur;
    // for now, still assuming that we start on the dot
  }
  frame_local++;
  advance(dur=(ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  f = xs+(2*n*frame_local);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double xt=(q[i].x-cx)*cl_im+cx_im-f[j], yt=(q[i].y-cy)*cl_im+cy_im-f[j+1], rsq=xt*xt+yt*yt;
    res_acc_local += pos_res[frame_local][i] = rsq;
    INTpos_res[i][1] = INTpos_res[i][2] + 0.5*(rsq)*dur;
  }
  pos_res_acc+=res_acc_local;
  return frame_local+1;
}

void runner::detect_events(record * rec_, int start_, int end_)
{
  reset_sim(rec_->params, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  for (int i = 0; i < n; i++) kill_frames[i] = 0;

  frame = start_detection(start_);
  double t_i=ts[frame-1], t_m1=ts[frame-2], res_acc_local=pos_res_acc;;
  do
  {
    double t_p1 = ts[frame];
    advance(dur=(t_p1-t_i)/t_phys, d_ang[frame-1], comega_s[frame], dt_sim);
    double * f = xs+(2*n*frame);

    bool all_dead=true;
    for (int i=0, j=0; i < n; i++, j+=2)
    {
      if (!(kill_frames[i])) // if this bead does not yet have a kill frame, continue search
      {
        all_dead=false; // then we still have at least one
        double xt=(q[i].x-cx)*cl_im + cx_im - f[j], yt=(q[i].y-cy)*cl_im + cy_im - f[j+1], rsq=xt*xt+yt*yt;

        pos_res[frame][i] = rsq;
        INTpos_res[i][0] = INTpos_res[i][1] + 0.5*(rsq)*dur;
        double alpha_it=alpha_INTpos_res[frame-1][i]=alpha_comp(INTpos_res[i], t_m1, t_p1);
        if (alpha_it > alpha_tol) // if we just had a collision
        {
          kill_frames[i] = frame-1;
          alpha_kill[i] = alpha_it;
        }
        else res_acc_local+=rsq;
      }
    }
    if (all_dead) break;
    else if (++frame==end_)
    {
      for (int i=0; i < n; i++) if (!(kill_frames[i]))
      {
        kill_frames[i] = end_-1;
        alpha_kill[i] = alpha_INTpos_res[end_-2][i];
      }
      break;
    }
    t_m1 = t_i; t_i = t_p1;
  } while(true);
  for (int i = 0; i < n; i++) event_frame_count[i][kill_frames[i]]++;
  rec_->record_event_data(pos_res_acc=res_acc_local, kill_frames, INTpos_res, alpha_kill);
}

void runner::run_relay(record * rec_, int start_, int * end_, int latest_, double residual_worst_)
{
  reset_sim(rec_->params, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  for (int i = 0; i < n; i++) pos_res[start_][i]= 0.0;
  double res_acc_local = 0.0;
  for (frame = start_+1; frame < latest_; frame++)
  {
    advance((ts[frame]-ts[frame-1])/t_phys, d_ang[frame-1], comega_s[frame], dt_sim);
    double * f = xs+(2*n*frame);
    for (int i=0, j=0; i < n; i++, j+=2) if (i<end_[i])
      {
        double xt=(q[i].x-cx)*cl_im + cx_im - f[j], yt=(q[i].y-cy)*cl_im + cy_im - f[j+1];
        res_acc_local += pos_res[frame][i] = xt*xt+yt*yt;
      }
    if (res_acc_local>residual_worst_) break;
  }
  return (int)rec_->check_success(pos_res_acc=res_acc_local,residual_worst_);
}
