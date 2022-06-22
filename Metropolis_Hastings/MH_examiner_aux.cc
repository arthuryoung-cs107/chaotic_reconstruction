#include "MH_workers.hh"

//MH_examiner

void MH_examiner::start_detecting_events(event_record * rec_, double * t_history_ double &net_r2_local_)
{
  thread_worker::reset_sim(rec_->u, ts[0]/t_phys, d_ang[0], comega_s[0], xs);

  int f_local=0,
      poffset=0,
      foffset=0;
  double net_r2_local=0.0, *pref = xs;

  int * f_event=rec_->evframe_bead;
  double  *r2_stable=rec_->r2stable_bead,
          *r2_unstable=rec_->r2unstable_bead;

  // initialize t=tstart data
  t_history_[0]=t_history_[1]=ts[0];
  for (int i = 0, j = 0; i < nbeads; i++,j+=2)
  {
    double  x_sim=q[i].x, y_sim=q[i].y,
            x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im;
    net_r2_local+=r2_state_comp[0][i]=0.0;

    INTr2_comp_history[i][0]=INTr2_comp_history[i][1]=0.0;
    alpha_state_comp[0][i]=NAN;

    // initialize component aggregates
    f_event[i]=0;
    r2_stable[i]=r2_unstable[i]=0.0;
  }
  // step to frame 1, begin computing divergence from data
  f_local++; poffset+=ndof; foffset+=nbeads; pref+=ndof;
  advance((ts[f_local]-ts[f_local-1])/t_phys, d_ang[f_local-1], comega_s[f_local], dt_sim);
  t_history_[2]=t_history_[1]; t_history_[1]=t_history_[0]; t_history_[0]=ts[f_local];
  for (int i = 0, j = 0; i < nbeads; i++,j+=2)
  {
    double  x_sim=q[i].x, y_sim=q[i].y,
            x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im,
            x_ref=pref[j], y_ref=pref[j+1],
            xerr=x_now-x_ref, yerr=y_now-y_ref, rsq=xerr*xerr+yerr*yerr;
    net_r2_local+=r2_state_comp[f_local][i]=rsq;

    INTr2_comp_history[i][2]=INTr2_comp_history[i][1]; INTr2_comp_history[i][1]=INTr2_comp_history[i][0];
    INTr2_comp_history[i][0]+=0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history_[0]-t_history_[1]);
    alpha_state_comp[f_local][i]=NAN;

    r2_stable[i]+=rsq;
  }

  // step to frame 2, compute divergence from data as we will in the body of the detection
  f_local++; poffset+=ndof; foffset+=nbeads; pref+=ndof;
  advance((ts[f_local]-ts[f_local-1])/t_phys, d_ang[f_local-1], comega_s[f_local], dt_sim);
  t_history_[2]=t_history_[1]; t_history_[1]=t_history_[0]; t_history_[0]=ts[f_local];
  for (int i = 0, j = 0; i < nbeads; i++,j+=2)
  {
    double  x_sim=q[i].x, y_sim=q[i].y,
            x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im,
            x_ref=pref[j], y_ref=pref[j+1],
            xerr=x_now-x_ref, yerr=y_now-y_ref, rsq=xerr*xerr+yerr*yerr;
    net_r2_local+=r2_state_comp[f_local][i]=rsq;

    INTr2_comp_history[i][2]=INTr2_comp_history[i][1]; INTr2_comp_history[i][1]=INTr2_comp_history[i][0];
    INTr2_comp_history[i][0]+=0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history_[0]-t_history_[1]);

    r2_stable[i]+=rsq;
  }
  net_r2_local_=net_r2_local;
  return f_local;
}

void MH_examiner::consolidate_examiner_event_data()
{
  int s_it=0;
  do
  {
    bool found_all=true;
    for ( i = 0; i < ncomp; i++) if (!(stev_comp[i]))
    {
      found_all=false;
      if (nev_state_comp[s_it][i]) stev_comp[i]=s_it;
    }
    if (found_all) break;
    else if (s_it++==stev_latest)
    {
      printf("(MH_examiner::consolidate_event_data) this shouldn't happen\n");
      for ( i = 0; i < ncomp; i++) if (!(order_comp[i])) stev_comp[i]=stev_latest;
      break;
    }
  } while (true);
}

bool MH_examiner::report_examiner_event_data(bool first2finish_, int &stev_earliest_, int &stev_latest_, int *stev_c_, int ** nev_s_c_, int ** nobs_s_c_, double ** r2_s_c_, double **alpha_s_c_)
{
  event_block::report_event_data(nev_s_c_,nobs_s_c_,r2_s_c_,alpha_s_c_);
  if (first2finish_)
  {
    stev_earliest_=stev_earliest;
    stev_latest_=stev_latest;
    for (int i = 0; i < ncomp; i++) stev_c_[i]=stev_comp[i];
  }
  else
  {
    if (stev_earliest<stev_earliest_) stev_earliest_=stev_earliest;
    if (stev_latest<stev_latest_) stev_latest_=stev_latest;
    for (int i = 0; i < ncomp; i++) if (stev_comp[i]<stev_c_[i]) stev_c_[i]=stev_comp[i];
  }
  return false;
}

void MH_examiner::restore_event_record(event_record *rec_, double *r2_Fb_, double *alpha_Fb_)
{
  int frame_last=stev_ordered[nbeads-1];
  double  net_r2_local=0.0,
          net_r2_stable_local=0.0,
          net_r2_regime_local=0.0,
          net_r2_unstable_local=0.0,
          *r2_stable=rec_->r2stable_bead,
          *net_r2_regime=rec_->netr2_regime,
          *r2_unstable=rec_->r2unstable_bead;

  for (int i = 0; i < nbeads; i++) r2_stable[i]=net_r2_regime[i]=r2_unstable[i]=0.0;

  for (int i_frame=0,k=0; i_frame <= frame_last; i_frame++)
    for (int i_bead = 0; i_bead < nbeads; i_bead++,k++)
    {
      double r2_it = r2_Fb_[k];
      net_r2_local+=r2_it;
      if (i_frame<=stev_comp[i_bead])
      {
        net_r2_stable_local+=r2_it; r2_stable[i_bead]+=r2_it;
        if (i_frame==stev_comp[i_bead])
        {
          rec_->evframe_bead[i_bead]=i_frame;
          rec_->alpha_bead[i_bead]=alpha_Fb_[k];
        }
      }
      else {net_r2_unstable_local+=r2_it; r2_unstable[i_bead]+=r2_it;}
    }
  net_r2_regime_local=net_r2_stable;
  rec_->record_event_data(&nf_obs,&net_r2_local);
}
