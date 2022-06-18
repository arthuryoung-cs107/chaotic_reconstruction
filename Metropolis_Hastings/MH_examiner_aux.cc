#include "MH_workers.hh"

//MH_examiner

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

void MH_examiner::restore_event_record(event_record *rec_, double *r2_Fb_, double *alpha_Fb_)
{
  int frame_last=stev_ordered[nbeads-1];
  double  net_r2_local=0.0,
          net_r2_stable_local=0.0,
          net_r2_unstable_local=0.0,
          *r2_stable=rec_->r2stable_bead,
          *r2_unstable=rec_->r2unstable_bead;

  for (int i = 0; i < nbeads; i++) r2_stable[i]=r2_unstable[i]=0.0;

  for (int i_frame=0,k=0; i_frame <= frame_last; i_frame++)
    for (int i_bead = 0; i_bead < nbeads; i_bead++,k++)
    {
      double r2_it = r2_Fb_[k];
      net_r2_local+=r2_it;
      if (i_frame<=stev_comp[i_bead])
      {
        net_r2_stable_local+=r2_it; r2_stable[i_bead]+=r2_it;
        if (i_frame==stev_comp[i])
        {
          rec_->evframe_bead[i_bead]=i_frame;
          rec_->alpha_bead[i_bead]=alpha_Fb_[k];
        }
      }
      else {net_r2_unstable_local+=r2_it; r2_unstable[i_bead]+=r2_it;}
    }
  rec_->record_event_data(&nf_obs,&net_r2_local);
}
