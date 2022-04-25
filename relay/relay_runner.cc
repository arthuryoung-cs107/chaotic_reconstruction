#include "particle_walk.hh"

extern "C"
{
  #include "AYaux.h"
}
runner::runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_, int param_len_, int Frames_, double tol_, double t_phys_, double dt_sim_, double t_wheels_, double *ts_, double *xs_, double *d_ang_, int ic_index_, int nlead_, int npool_) : swirl(sp_, pg_, wl_, n_),
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
runner::~runner()
{
  delete [] lead_dup_count; delete [] frame_kill_count;
  delete [] frame_res_data; delete [] param_mean;
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
  frame = start_+1;
  double dur1 = advance_runner();
  double *f = xs+(2*n*frame);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double xt=(q[i].x-cx)*cl_im + cx_im - f[j], yt=(q[i].y-cy)*cl_im + cy_im - f[j+1];
    pos_res[start_][i]=INTpos_res[start_][i]=0.0; // zero out the first frame
    pos_res[frame][i] = xt*xt+yt*yt;
    INTpos_res[i][2] = 0.5*(pos_res[frame][i])*dur1;
    // for now, still assuming that we start on the dot
  }
  frame++;
  double dur2 = advance_runner();
  f = xs+(2*n*frame);
  for (int i=0, j=0; i < n; i++, j+=2)
  {
    double xt=(q[i].x-cx)*cl_im + cx_im - f[j], yt=(q[i].y-cy)*cl_im + cy_im - f[j+1];
    pos_res[frame][i] = xt*xt+yt*yt;
    INTpos_res[i][1] = INTpos_res[i][2] + 0.5*(pos_res[frame][i])*dur2;
  }
  return frame+1;
}

void runner::detect_events(record * rec_, int start_, int end_)
{
  reset_sim(rec_->params, ts[start_]/t_phys, d_ang[start_], comega_s[start_], xs + 2*n*start_);
  for (int i = 0; i < n; i++) kill_frames[i] = 0;

  frame = start_detection(start_);

  double t_i=ts[start_+1], t_m1=ts[start_];
  do
  {
    double t_p1 = ts[frame], dur=(t_p1-t_i)/t_phys, ctheta=d_ang[frame-1], comega=d_ang[frame]-d_ang[frame-1];

    if(comega>M_PI) comega-=2*M_PI; else if(comega<-M_PI) comega+=2*M_PI; comega/=dur;
    advance(dur, ctheta, comega, dt_sim);
    double * f = xs+(2*n*frame);

    bool all_dead=true;
    for (int i=0, j=0; i < n; i++, j+=2)
    {
      if (!(kill_frames[i])) // if this bead does not yet have a kill frame, continue search
      {
        all_dead=false; // then we still have at least one
        double xt=(q[i].x-cx)*cl_im + cx_im - f[j], yt=(q[i].y-cy)*cl_im + cy_im - f[j+1], alpha_it;
        pos_res[frame][i] = xt*xt+yt*yt;
        INTpos_res[i][0] = INTpos_res[i][1] + 0.5*(pos_res[frame][i])*dur;
        alpha_INTpos_res[frame-1][i] = alpha_comp(INTpos_res[i], t_m1, t_p1);
        if (alpha_INTpos_res[frame-1][i] > alpha_tol)
        {
          kill_frames[i] = frame-1;
          alpha_kill[i] = alpha_INTpos_res[kill_frames[i]][i];
        }
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
  rec_->record_event_data(kill_frames, INTpos_res, alpha_kill);
}
