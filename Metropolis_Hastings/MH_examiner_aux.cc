#include "MH_workers.hh"

//MH_examiner

void MH_examiner::start_detecting_events(int &f_local_,int &iregime_local_,int *f_event_,double &netr2_local_,double &netr2_stable_local_,double &netr2_unstable_local_,double *t_history_,double *r2stable_bead_,double *netr2_regime_,double *r2unstable_bead_,double *alphaev_bead_)
{
  // initialize everything
  f_local_=iregime_local_=0;
  netr2_local_=netr2_stable_local_=netr2_unstable_local_=0.0;
  t_history_[0]=t_history_[1]=t_history_[2]=ts[f_local_];
  for (int i = 0; i < nbeads; i++)
  {
    r2_state_comp[f_local_][i]=INTr2_comp_history[i][0]=INTr2_comp_history[i][1]=r2stable_bead_[i]=netr2_regime_[i]=r2unstable_bead_[i]=0.0;
    alpha_state_comp[f_local_][i]=alphaev_bead_[i]=NAN;
    f_event_[i]=0;
  }

  // resolve frames 1 and 2
  double *pref;
  for (int iframe = 1; iframe <= 2; iframe++)
  {
    pref=thread_worker::advance_sim(++f_local_,t_history_);
    for (int i=0,j=0; i < nbeads; i++,j+=dof)
    {
      double rsq=thread_worker::compute_residual(q[i].x,q[i].y,pref[j],pref[j+1]);
      netr2_local_+=r2_state_comp[f_local_][i]=rsq;

      basic_thread_worker::update_integral_history(0.5*(rsq+r2_state_comp[f_local_-1][i])*(t_history_[0]-t_history_[1]),i);

      alpha_state_comp[f_local_][i]=NAN;

      r2stable_bead_[i]+=rsq;
    }
  }
  // set stable residual, assuming that frames 0,1,2 were all stable
  netr2_stable_local_=netr2_regime_[0]=netr2_local_;
}

void MH_examiner::update_event_data(int final_frame_, int *f_event_, double *r2i_, double *alphai_)
{

  ntest++;
  stev_early=stev_late=f_event_[0];
  nf_stable=0;
  for (int i=0, k=0; i < nbeads; i++)
  {
    int fevent_it = f_event_[i];
    nf_stable+=fevent_it;
    nev_state_comp[fevent_it][i]++;
    if (fevent_it<stev_early) stev_early = fevent_it;
    if (fevent_it>stev_late) stev_late = fevent_it;
    for (int j = 0; j < final_frame_; j++,k++)
    {
      mur2_state_comp[0][k]+=r2i_[k]=r2_state_comp[0][k];
      mualpha_state_comp[0][k]+=alphai_[k]=alpha_state_comp[0][k];
      nobs_state_comp[0][k]++;
    }
  }
  nf_obs=stev_late*nbeads;
  nf_regime=nf_stable;
  nf_unstable=nf_obs-nf_stable;
  if (stev_early<stev_earliest) stev_earliest=stev_early;
  if (stev_late>stev_latest) stev_latest=stev_late;
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


bool MH_examiner::update_training_data(int i_, double r2success_threshold_)
{
  ntest++;
  bool success_local=(net_r2_regime)<r2success_threshold_;
  if (success_local) int_wkspc[nsuccess_test++] = i_;
  return success_local;
}

void
