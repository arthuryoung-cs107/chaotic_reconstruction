#include "MH_workers.hh"

//MH_examiner

void MH_examiner::detect_events(event_record *rec_, double *r2i_, double *alphai_)
{
  thread_worker::reset_sim(rec_->u,ts[0]/t_phys, d_ang[0], comega_s[0], xs);

  int f_local,
      *f_event=rec_->evframe_bead;

  double  netr2_local,
          netr2_stable_local,
          netr2_unstable_local,
          t_history[3],
          *r2net_bead=rec_->r2net_bead,
          *r2stable_bead=rec_->r2stable_bead,
          *r2unstable_bead=rec_->r2unstable_bead,
          *alphaev_bead=rec_->alpha_bead,
          *pref;

  start_detecting_events(f_local,f_event,netr2_local,netr2_stable_local,netr2_unstable_local,t_history,r2net_bead,r2stable_bead,r2unstable_bead,alphaev_bead);

  do
  {
    pref=thread_worker::advance_sim(++f_local,t_history);
    bool all_events_detected=true;
    for (int i=0,j=0; i < nbeads; i++,j+=dim)
    {
      double rsq=thread_worker::compute_residual(q[i].x,q[i].y,pref[j],pref[j+1]);
      r2_state_comp[f_local][i]=rsq;
      netr2_local+=rsq; r2net_bead[i]+=rsq;

      basic_thread_worker::update_integral_history(0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history[0]-t_history[1]),i);

      double alpha_it = alpha_state_comp[f_local-1][i]=alpha_comp(INTr2_comp_history[i], t_history[0], t_history[2]);
      if (!(f_event[i]))
      {
        all_events_detected=false;
        if (alpha_it>alpha_tol)
        {
          f_event[i]=f_local-1; alphaev_bead[i]=alpha_it;
          netr2_unstable_local+=r2unstable_bead[i]=rsq;
        }
        else
          {netr2_stable_local+=rsq; r2stable_bead[i]+=rsq;}
      }
      else
        {netr2_unstable_local+=rsq; r2unstable_bead[i]+=rsq;}
    }
    if (all_events_detected) break;
    else if ((f_local+1)==Frames)
    {
      for (int i = 0; i < nbeads; i++) if (!(f_event[i]))
      {
        f_event[i] = f_local-1;
        alphaev_bead[i] = alpha_state_comp[f_local-1][i];
      }
      break;
    }
  } while(true);
  netr2=netr2_local;
  netr2_stable=netr2_stable_local;
  netr2_unstable=netr2_unstable_local;

  // set thread worker statistics data
  update_event_data(f_local, f_event, r2i_, alphai_);
  // set record data
  rec_->record_event_data(&nf_netobs,&netr2);
}

bool MH_examiner::examine_u(event_record *rec_, int i_, double r2success_threshold_)
{
  thread_worker::reset_sim(rec_->u, ts[0]/t_phys, d_ang[0], comega_s[0], xs);
  // run the length of the event block, storing bead positions
  for (int i_f = 1; i_f <= stev_ordered[nbeads-1]; i_f++)
  {
    advance((ts[i_f]-ts[i_f-1])/t_phys, d_ang[i_f-1], comega_s[i_f], dt_sim);
    for (int i=0,j=dof*i_f; i < nbeads; i++,j+=dim) {psim[j]=q[i].x; psim[j+1]=q[i].y;}
  }

  // process the event block
  double  netr2_local=0.0,
          netr2_stable_local=0.0,
          netr2_unstable_local=0.0,
          *pref=xs,
          *r2net_bead=rec_->r2net_bead,
          *r2stable_bead=rec_->r2stable_bead,
          *r2unstable_bead=rec_->r2unstable_bead;
  rec_->clear_residuals();

  // purely stable residuals
  for (int i_f = 1; i_f <= stev_ordered[0]; i_f++)
  {
    for (int i=0,j=dof*i_f; i < nbeads; i++,j+=dim)
    {
      double rsq=thread_worker::compute_residual(psim[j], psim[j+1], pref[j], pref[j+1]);
      netr2_local+=rsq; r2net_bead[i]+=rsq;
      netr2_stable_local+=rsq; r2stable_bead[i]+=rsq;
    }
  }
  // entering unstable sequence components
  for (int i_f = stev_ordered[0]+1; i_f <= stev_ordered[nbeads-1]; i_f++)
  {
    for (int i=0,j=dof*i_f; i < nbeads; i++,j+=dim)
    {
      double rsq=thread_worker::compute_residual(psim[j], psim[j+1], pref[j], pref[j+1]);
      netr2_local+=rsq; r2net_bead[i]+=rsq;
      if (i_f<=stev_comp[i]) // still stable
      {netr2_stable_local+=rsq; r2stable_bead[i]+=rsq;}
      else // unstable
      {netr2_unstable_local+=rsq; r2unstable_bead[i]+=rsq;}
    }
  }
  netr2=netr2_local;
  netr2_stable=netr2_stable_local;
  netr2_unstable=netr2_unstable_local;

  bool success_local = update_training_data(i_,r2success_threshold_);
  rec_->record_training_data(&netr2,success_local);
  return success_local;
}
