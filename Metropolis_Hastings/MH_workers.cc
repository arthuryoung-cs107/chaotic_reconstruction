#include "MH_workers.hh"

//MH_examiner

void MH_examiner::detect_events(event_record *rec_, double *r2i_, double *alphai_)
{
  thread_worker::reset_sim(rec_->u,ts[0]/t_phys, d_ang[0], comega_s[0], xs);

  int f_local,
      iregime_local,
      *f_event=rec_->evframe_bead;

  double  netr2_local,
          netr2_stable_local,
          netr2_unstable_local,
          t_history[3],
          *r2stable_bead=rec_->r2stable_bead,
          *netr2_regime=rec_->netr2_regime,
          *r2unstable_bead=rec_->r2unstable_bead,
          *alphaev_bead=rec_->alpha_bead,
          *pref;

  start_detecting_events(f_local,iregime_local,f_event,t_history,netr2_local,netr2_stable_local,netr2_unstable_local,r2stable_bead,netr2_regime,r2unstable_bead,alphaev_bead);

  do
  {
    pref=thread_worker::advance_sim(++f_local,t_history);
    bool all_events_detected=true;
    for (int i=0,j=0; i < nbeads; i++,j+=2)
    {
      double rsq=thread_worker::compute_residual(q[i].x,q[i].y,pref[j],pref[j+1]);
      net_r2_local+=r2_state_comp[f_local][i]=rsq;

      basic_thread_worker::update_integral_history(0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history[0]-t_history[1]),i);

      double alpha_it = alpha_state_comp[f_local-1][i]=alpha_comp(INTr2_comp_history[i], t_history[0], t_history[2]);
      if (!(f_event[i]))
      {
        all_events_detected=false;
        if (alpha_it>alpha_tol)
        {
          f_event[i]=f_local-1;
          alphaev_bead[i]=alpha_it;
          net_r2_unstable_local+=r2_unstable[i]=netr2_regime[++iregime_local]=rsq;
        }
        else
        {
          r2_stable[i]+=rsq;
          net_r2_stable_local+=rsq;
          netr2_regime[iregime_local]+=rsq;
        }
      }
      else
      {
        r2_unstable[i]+=rsq;
        net_r2_unstable_local+=rsq;
        netr2_regime[iregime_local]+=rsq;
      }
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
  net_r2=net_r2_local;
  net_r2_stable=net_r2_stable_local;
  net_r2_unstable=net_r2_unstable_local;

  // set thread worker statistics data
  update_event_data(f_local, f_event, r2i_, alphai_);
  // set record data
  rec_->record_event_data(&nf_obs,&net_r2);
}

bool MH_examiner::examine_u(event_record *rec_, int i_, double r2success_threshold_)
{
  thread_worker::reset_sim(rec_->u, ts[0]/t_phys, d_ang[0], comega_s[0], xs);
  int poffset=0;
  // run the length of the event block, storing bead positions
  for (int i_f = 1; i_f <= stev_latest; i_f++)
  {
    poffset+=ndof;
    advance((ts[i_f]-ts[i_f-1])/t_phys, d_ang[i_f-1], comega_s[i_f], dt_sim);
    for (int i = 0, j = poffset; i < nbeads; i++,j+=2) {psim[j]=q[i].x; psim[j+1]=q[i].y;}
  }

  // process the event block
  double  net_r2_local=0.0,
          net_r2_stable_local=0.0,
          net_r2_regime_local=0.0,
          net_r2_unstable_local=0.0,
          *pref=xs,
          *r2_stable_comp=rec_->r2stable_bead,
          *netr2_regime=rec_->netr2_regime,
          *r2_unstable_comp=rec_->r2unstable_bead;
  poffset=0;
  clear_examiner_residuals(r2_stable_comp, netr2_regime, r2_unstable_comp);

  // purely stable residuals
  for (int i_f = 1; i_f <= stev_earliest; i_f++)
  {
    pref+=ndof;
    for (int i_c = 0; i_c < nbeads; i_c++)
    {
      double  x_sim=psim[j], y_sim=psim[j+1],
              x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im,
              x_ref=pref[j], y_ref=pref[j+1],
              xerr=x_now-x_ref, yerr=y_now-y_ref, net_r2_local+=rsq=xerr*xerr+yerr*yerr;
      net_r2_stable_local+=rsq; r2_stable_comp[i_c]+=rsq;
    }
  }
  netr2_regime[0]=net_r2_local;
  if (iregime_active==0) net_r2_regime_local=net_r2_local;
  // regime region
  for (int i_regime = 1; i_regime < nbeads; i_regime++)
  {
    double r2_regime_it=0.0;
    for (int i_f = stev_ordered[i_regime-1]+1; i_f <= stev_ordered[i_regime]; i_f++)
    {
      pref+=ndof;
      for (int i_c = 0, j = 0; i_c < nbeads; i_c++,j+=2)
      {
        double  x_sim=psim[j], y_sim=psim[j+1],
                x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im,
                x_ref=pref[j], y_ref=pref[j+1],
                xerr=x_now-x_ref, yerr=y_now-y_ref, net_r2_local+=rsq=xerr*xerr+yerr*yerr;
        if (i_f<=stev_comp[i_c]) // still stable
        {net_r2_stable_local+=rsq; r2_stable_comp[i_c]+=rsq;}
        else // unstable
        {r2_regime_it+=rsq; net_r2_unstable_local+=rsq; r2_unstable_comp[i_c]+=rsq;}
      }
    }
    netr2_regime[i_regime]=r2_regime_it;
    if (i_regime<=iregime_active) net_r2_regime_local+=r2_regime_it;
  }
  net_r2=net_r2_local;
  net_r2_stable=net_r2_stable_local;
  net_r2_regime=net_r2_regime_local;
  net_r2_unstable=net_r2_unstable_local;

  bool success_local = update_training_data(i_,r2success_threshold_);
  rec_->record_training_data(&net_r2,success_local);
  return success_local;
}
