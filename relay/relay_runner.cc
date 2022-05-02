#include "particle_relay.hh"

extern "C"
{
  #include "AYaux.h"
}
runner::runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, int nlead_, int npool_, double t_phys_, double dt_sim_, double alpha_tol_, double *ts_, double *xs_, double *d_ang_, double *comega_s_) : swirl(sp_, pg_, wl_, n_),
thread_id(thread_id_), param_len(param_len_), Frames(Frames_), nlead(nlead_), npool(npool_),
t_phys(t_phys_), dt_sim(dt_sim_), alpha_tol(alpha_tol_),
ts(ts_), xs(xs_), d_ang(d_ang_), comega_s(comega_s_),
lead_dup_count(new int[nlead_]), kill_frames(new int[n_]),
event_frame_count(AYimatrix(n_, Frames_)),
param_acc(new double[param_len_]), alpha_kill(new double[n_]),
pos_res(AYdmatrix(Frames_, n_)), alpha_INTpos_res(AYdmatrix(Frames_, n_)), INTpos_res(AYdmatrix(n_, 3))
{pvals = &Kn;}

runner::~runner()
{
  delete [] lead_dup_count; delete [] kill_frames; delete [] param_acc; delete [] alpha_kill;
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
int runner::start_detection(int start_, double * params_, double *t_history_)
{
  for (int i = 0; i < n*Frames; i++) pos_res[0][i]=alpha_INTpos_res[0][i]=0.0;

  reset_sim(params_, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  int frame_local = start_, poffset = 2*n*frame_local, foffset=n*frame_local;
  double res_acc_local = 0.0, *f = xs+poffset;

  // initialize t=tstart data
  t_history_[0]=t_history_[1]=ts[frame_local];
  for (int i=0,j=0; i < n; i++,j+=2)
  {
    kill_frames[i] = 0;
    pos_res[frame_local][i]=alpha_INTpos_res[frame_local][i]=INTpos_res[i][0]=INTpos_res[i][1]=INTpos_res[i][2]=0.0;

    // write test outputs
    double  x_sim = q[i].x, y_sim = q[i].y,
            x_now = (x_sim-cx)*cl_im + cx_im, y_now = (y_sim-cy)*cl_im + cy_im;
  }
  // step to frame 1, begin computing divergence from data
  frame_local++; poffset+=2*n; foffset+=n; f = xs+poffset;
  advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  t_history_[2]=t_history_[1]; t_history_[1]=t_history_[0]; t_history_[0]=ts[frame_local];
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  xt=(q[i].x-cx)*cl_im + cx_im-f[j],yt=(q[i].y-cy)*cl_im + cy_im-f[j+1], rsq=xt*xt+yt*yt;
    res_acc_local+=pos_res[frame_local][i]=rsq;
    INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
    INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history_[0]-t_history_[1]);
    alpha_INTpos_res[frame_local][i]=0.0;
  }

  // step to frame 2, compute divergence from data as we will in the body of the detection
  frame_local++; poffset+=2*n; foffset+=n; f = xs+poffset;
  advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
  t_history_[2]=t_history_[1]; t_history_[1]=t_history_[0]; t_history_[0]=ts[frame_local];
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double  xt=(q[i].x-cx)*cl_im + cx_im-f[j],yt=(q[i].y-cy)*cl_im + cy_im-f[j+1], rsq=xt*xt+yt*yt;
    res_acc_local+=pos_res[frame_local][i]=rsq;
    INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
    INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history_[0]-t_history_[1]);
    alpha_INTpos_res[frame_local][i]=0.0;
  }
  pos_res_acc=res_acc_local;
  return frame_local;
}

void runner::detect_events(record * rec_, int start_, int end_)
{
  double t_history[3];
  int frame_local = start_detection(start_, rec_->params, t_history);
  int poffset = 2*n*frame_local, foffset=n*frame_local;
  double res_acc_local = pos_res_acc, res_acc_local_full = pos_res_acc, *f;
  do
  {
    frame_local++; poffset+=2*n; foffset+=n; f = xs+poffset;
    advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
    t_history[2]=t_history[1]; t_history[1]=t_history[0]; t_history[0]=ts[frame_local];

    bool all_dead=true;

    for (int i=0, j=0; i < n; i++, j+=2)
    {
      if (!(kill_frames[i])) // if this bead does not yet have a kill frame, continue search
      {
        all_dead=false;

        double  xt=(q[i].x-cx)*cl_im + cx_im-f[j],yt=(q[i].y-cy)*cl_im + cy_im-f[j+1], rsq=xt*xt+yt*yt;
        res_acc_local+=pos_res[frame_local][i]=rsq;
        INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
        INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history[0]-t_history[1]);
        alpha_INTpos_res[frame_local-1][i]=alpha_comp(INTpos_res[i], t_history[0], t_history[2]);
        if (alpha_INTpos_res[frame_local-1][i] > alpha_tol) // if we just had a collision
        {
          kill_frames[i] = frame_local-1;
          alpha_kill[i] = alpha_INTpos_res[frame_local-1][i];
        }
        else res_acc_local+=rsq;
      }
    }
    if (all_dead) break;
    else if ((frame_local+1)==end_)
    {
      for (int i=0; i < n; i++) if (!(kill_frames[i]))
        {kill_frames[i] = end_-1; alpha_kill[i] = alpha_INTpos_res[frame_local-2][i];}
      break;
    }
  } while(true);
  for (int i = 0; i < n; i++) event_frame_count[i][kill_frames[i]]++;
  rec_->record_event_data(pos_res_acc=res_acc_local, kill_frames, INTpos_res, alpha_kill);
}

int runner::run_relay(record * rec_, int start_, int * end_, int earliest_, int latest_, double residual_worst_)
{
  reset_sim(rec_->params, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  for (int i = 0; i < n; i++) pos_res[start_][i]= 0.0;
  double res_acc_local = 0.0;
  for (frame = start_+1; frame < latest_; frame++)
  {
    advance((ts[frame]-ts[frame-1])/t_phys, d_ang[frame-1], comega_s[frame], dt_sim);
    double * f = xs+(2*n*frame);
    for (int i=0, j=0; i < n; i++, j+=2)
    {
      if (i<end_[i])
      {
        double xt=(q[i].x-cx)*cl_im + cx_im - f[j], yt=(q[i].y-cy)*cl_im + cy_im - f[j+1];
        res_acc_local += pos_res[frame][i] = xt*xt+yt*yt;        
      }
    }
    if (res_acc_local>residual_worst_) break;
  }
  return (int)rec_->check_success(pos_res_acc=res_acc_local,residual_worst_);
}
