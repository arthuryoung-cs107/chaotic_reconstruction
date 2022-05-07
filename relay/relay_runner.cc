#include "particle_relay.hh"

extern "C"
{
  #include "AYaux.h"
}
runner::runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, int nlead_, int npool_, double t_phys_, double dt_sim_, double alpha_tol_, double *ts_, double *xs_, double *d_ang_, double *comega_s_) : swirl(sp_, pg_, wl_, n_),
thread_id(thread_id_), param_len(param_len_), Frames(Frames_), nlead(nlead_), npool(npool_),
t_phys(t_phys_), dt_sim(dt_sim_), alpha_tol(alpha_tol_),
ts(ts_), xs(xs_), d_ang(d_ang_), comega_s(comega_s_),
lead_dup_count(new int[nlead_]),
event_frame_count(AYimatrix(n_, Frames_)),
param_acc(new double[param_len_]),
pos_res(AYdmatrix(Frames_, n_)), alpha_INTpos_res(AYdmatrix(Frames_, n_)), INTpos_res(AYdmatrix(n_, 3)), sim_pos(new double[2*n_*Frames_]), bead_event_res_mat(new double[n_*n_])
{pvals = &Kn;}

runner::~runner()
{
  delete [] lead_dup_count; delete [] param_acc; delete [] bead_event_res_mat; delete [] sim_pos;
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
int runner::start_detection(int start_, int min_, double * params_, double *t_history_, int * event_frames, double * smooth_residual, double * net_residual_)
{
  for (int i = 0; i < n*Frames; i++) pos_res[0][i]=alpha_INTpos_res[0][i]=0.0;

  reset_sim(params_, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  int frame_local = start_, poffset = 2*n*frame_local, foffset=n*frame_local;
  double net_residual = 0.0, *f = xs+poffset;

  // initialize t=tstart data
  t_history_[0]=t_history_[1]=ts[frame_local];
  for (int i=0,j=0; i < n; i++,j+=2)
  {
    event_frames[i] = 0;
    smooth_residual[i]=pos_res[frame_local][i]=alpha_INTpos_res[frame_local][i]=INTpos_res[i][0]=INTpos_res[i][1]=INTpos_res[i][2]=0.0;

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
    net_residual+=pos_res[frame_local][i]=rsq;
    smooth_residual[i]+=rsq;
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
    net_residual+=pos_res[frame_local][i]=rsq;
    smooth_residual[i]+=rsq;
    INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
    INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history_[0]-t_history_[1]);
    alpha_INTpos_res[frame_local][i]=0.0;
  }

  if (min_>start_)
  {
    for ( frame_local = frame_local+1 ; frame_local <= min_; frame_local++)
    {
      poffset+=2*n; foffset+=n; f = xs+poffset;
      advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
      t_history_[2]=t_history_[1]; t_history_[1]=t_history_[0]; t_history_[0]=ts[frame_local];
      for (int i=0, j=0; i < n; i++, j+=2)
      {
        double  xt=(q[i].x-cx)*cl_im + cx_im-f[j],yt=(q[i].y-cy)*cl_im + cy_im-f[j+1], rsq=xt*xt+yt*yt;
        net_residual+=pos_res[frame_local][i]=rsq;
        smooth_residual[i]+=rsq;
        INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
        INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history_[0]-t_history_[1]);
        alpha_INTpos_res[frame_local-1][i]=alpha_comp(INTpos_res[i], t_history_[0], t_history_[2]);
      }
    }
    frame_local--;
  }

  *net_residual_=net_residual;
  return frame_local;
}

void runner::detect_events(record * rec_, int start_, int min_, int end_)
{
  double t_history[3], net_residual;
  int frame_local = start_detection(start_, min_, rec_->params, t_history, rec_->event_positions, rec_->smooth_residual, &net_residual);
  int poffset = 2*n*frame_local, foffset=n*frame_local, *event_frames=rec_->event_positions;
  double *f, *rec_res=rec_->smooth_residual, *alpha_data=rec_->alpha_data, net_smooth_residual=net_residual;

  do
  {
    frame_local++; poffset+=2*n; foffset+=n; f = xs+poffset;
    advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
    t_history[2]=t_history[1]; t_history[1]=t_history[0]; t_history[0]=ts[frame_local];

    bool all_dead=true;

    for (int i=0, j=0; i < n; i++, j+=2)
    {
      double  xt=(q[i].x-cx)*cl_im + cx_im-f[j],yt=(q[i].y-cy)*cl_im + cy_im-f[j+1], rsq=xt*xt+yt*yt;
      net_residual+=pos_res[frame_local][i]=rsq;
      if (!(event_frames[i])) // if this bead does not yet have a kill frame, continue search
      {
        all_dead=false;
        INTpos_res[i][2]=INTpos_res[i][1]; INTpos_res[i][1]=INTpos_res[i][0];
        INTpos_res[i][0]+=0.5*(rsq+pos_res[frame_local-1][i])*(t_history[0]-t_history[1]);
        alpha_INTpos_res[frame_local-1][i]=alpha_comp(INTpos_res[i], t_history[0], t_history[2]);
        if (alpha_INTpos_res[frame_local-1][i] > alpha_tol) // if we just had a collision
        {
          event_frames[i] = frame_local-1;
          alpha_data[i] = alpha_INTpos_res[frame_local-1][i];
          rec_res[i+n]+=rsq;
        }
        else
        {rec_res[i]+=rsq; net_smooth_residual+=rsq;}
      }
      else rec_res[i+n]+=rsq;
    }
    if (all_dead) break;
    else if ((frame_local+1)==end_)
    {
      for (int i=0; i < n; i++) if (!(event_frames[i]))
        {event_frames[i] = end_-1; alpha_data[i] = alpha_INTpos_res[frame_local-2][i];}
      break;
    }
  } while(true);
  for (int i = 0; i < n; i++) event_frame_count[i][event_frames[i]]++;
  bool done = rec_->check_success(net_residual, net_smooth_residual, DBL_MAX);
}

int runner::run_relay(record * rec_, int start_, int earliest_, int latest_, int * event_frames_ordered_, int * bead_order, double residual_worst_)
{
  int poffset = 2*n*start_;

  reset_sim(rec_->params, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + poffset);

  // set the starting frame position
  double *tsp = sim_pos+(poffset);
  for (int i=0, j=0; i < n; i++, j+=2)
  {tsp[j]=q[i].x; tsp[j+1]=q[i].y;}


  for (int i = 0; i < n*n; i++) bead_event_res_mat[i] = 0.0;

  int stretch_start = start_+1;
  for (int si = 0; si < n; si++)
  {
    double * berm_si = bead_event_res_mat+(si*n);
    for (int frame_local = stretch_start; frame_local < event_frames_ordered_[si]; frame_local++)
    {
      advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
      poffset+=2*n;
      double *f = xs+(poffset);
      for (int i=0, j=0; i < n; i++, j+=2)
      {
        double xt=(q[i].x-cx)*cl_im + cx_im - f[j], yt=(q[i].y-cy)*cl_im + cy_im - f[j+1];
        berm_si[i]+=xt*xt+yt*yt; // accumulate the event vs bead matrix.
      }
    }
    stretch_start=event_frames_ordered_[si];
  }

  double * rec_res=rec_->smooth_residual, net_residual=0.0, net_smooth_residual=0.0;

  // clear this record's residual data, which is split between stiff and smooth residual
  for (int i = 0; i < 2*n; i++) rec_res[i] = 0.0;

  // accumulate residuals across the full event block
  for (int event_i = 0; event_i < n; event_i++) // from earliest event to latest event
  {
    double  *berm_si = bead_event_res_mat + bead_order[event_i], // maps to the starting column index of the bead_res_mat
            *beadi_res = rec_res + bead_order[event_i];

    // accumulate smooth event residual for current bead
    for (int si = 0; si < event_i+1; si++) *(beadi_res)+= berm_si[si*n];

    net_smooth_residual += *(beadi_res);
    net_residual += *(beadi_res);
    beadi_res+=n;

    // accumulate stiff event residual for current bead
    for (int si = event_i+1; si < n; si++) *(beadi_res)+= berm_si[si*n];
    net_residual += *(beadi_res);
  }

  rec_res=rec_->stiff_residual;
  // if applicable, train on more of the data
  for (int frame_local = stretch_start; frame_local < latest_; frame_local++)
  {
    advance((ts[frame_local]-ts[frame_local-1])/t_phys, d_ang[frame_local-1], comega_s[frame_local], dt_sim);
    poffset+=2*n;
    double *f = xs+(poffset);
    for (int i=0, j=0; i < n; i++, j+=2)
    {
      double xt=(q[i].x-cx)*cl_im + cx_im - f[j], yt=(q[i].y-cy)*cl_im + cy_im - f[j+1], rsq=xt*xt+yt*yt;
      net_residual+=rsq;
      rec_res[i]+=rsq;
    }
  }
  return (int)rec_->check_success(net_residual, net_smooth_residual, residual_worst_);
}
