#include "MH_workers.hh"

//MH_examiner

void MH_examiner::detect_events(event_record *rec_, double *r2i_, double *alphai_)
{
  double t_history[3], net_r2_local;
  int f_local = start_detecting_events(rec_, t_history, net_r2_local),
      poffset=ndof*f_local,
      foffset=nbeads*f_local,
      *f_event=rec_->evframe_bead;

  double  *pref=xs+poffset,
          *r2_stable=rec_->r2stable_bead,
          *r2_regime=rec_->netr2_regime,
          *r2_unstable=rec_->r2unstable_bead,
          *alpha_bead=rec_->alpha_bead,
          net_r2_stable_local=net_r2_local,
          net_r2_regime_local=net_r2_local,
          net_r2_unstable_local=0.0;
  do
  {
    f_local++; poffset+=ndof; pref+=ndof;
    advance((ts[f_local]-ts[f_local-1])/t_phys, d_ang[f_local-1], comega_s[f_local], dt_sim);
    t_history[2]=t_history[1]; t_history[1]=t_history[0]; t_history[0]=ts[f_local];
    bool all_events_detected=true;
    for (int i = 0, j = 0; i < nbeads; i++,j+=2)
    {
      double  x_sim=q[i].x, y_sim=q[i].y,
              x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im,
              x_ref=pref[j], y_ref=pref[j+1],
              xerr=x_now-x_ref, yerr=y_now-y_ref, rsq=xerr*xerr+yerr*yerr;
      r2_state_comp[f_local][i]=rsq; net_r2_local+=rsq;
      INTr2_comp_history[i][2]=INTr2_comp_history[i][1]; INTr2_comp_history[i][1]=INTr2_comp_history[i][0];
      INTr2_comp_history[i][0]+=0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history[0]-t_history[1]);
      double alpha_it = alpha_state_comp[f_local-1][i]=alpha_comp(INTr2_comp_history[i], t_history[0], t_history[2]);
      if (!(f_event[i]))
      {
        all_events_detected=false;
        if (alpha_it>alpha_tol)
        {
          f_event[i]=f_local-1;
          alpha_bead[i]=alpha_it;
          r2_unstable[i]=rsq;
          net_r2_unstable_local+=rsq;
        }
        else
        {
          r2_stable[i]+=rsq;
          net_r2_stable_local+=rsq;
        }
      }
      else
      {
        r2_unstable[i]+=rsq;
        net_r2_unstable_local+=rsq;
      }
    }
    if (all_events_detected) break;
    else if ((f_local+1)==Frames)
    {
      for (int i = 0; i < nbeads; i++) if (!(f_event[i]))
      {
        f_event[i] = f_local-1;
        alpha_bead[i] = alpha_state_comp[f_local-1][i];
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
  for (int i_f = 1, j=dof; i_f <= stev_earliest; i_f++)
  {
    pref+=ndof;
    for (int i_c = 0; i_c < nbeads; i_c++, j+=dof)
    {
      double  x_sim=psim[j], y_sim=psim[j+1],
              x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im,
              x_ref=pref[j], y_ref=pref[j+1],
              xerr=x_now-x_ref, yerr=y_now-y_ref, rsq=xerr*xerr+yerr*yerr;
      net_r2_local+=rsq;
      net_r2_stable_local+=rsq; r2_stable_comp[i_c]+=rsq;
    }
  }
  netr2_regime[0]=net_r2_local;
  if (iregime_active==0) net_r2_regime_local=net_r2_local;
  // regime region
  for (int i_regime = 1; i_regime < nbeads; i_regime++)
  {
    double r2_regime_it=0.0;
    for (int i_f = stev_ordered[i_regime-1]+1, j=dof*i_f; i_f <= stev_ordered[i_regime]; i_f++)
    {
      pref+=ndof;
      for (int i_c = 0, j = 0; i_c < nbeads; i_c++,j+=dof)
      {
        double  x_sim=psim[j], y_sim=psim[j+1],
                x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im,
                x_ref=pref[j], y_ref=pref[j+1],
                xerr=x_now-x_ref, yerr=y_now-y_ref, rsq=xerr*xerr+yerr*yerr;
        net_r2_local+=rsq;
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

// MH_medic

MH_medic::MH_medic(swirl_param  &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_, int Frames_test_, char * test_buffer_): MH_examiner(sp_,pg_,wl_,tws_,thread_id_, alpha_tol_),
Frames_test(Frames_test_),
buf_end(strlen(test_buffer_)), mtest_buffer(new char[buf_end+50]),
TEST_p(new double[ndof*Frames_test]), TEST_r2(new double[nbeads*Frames_test]), TEST_alpha(new double[nbeads*Frames_test]), TEST_INTr2(new double[nbeads*Frames_test])
{strcpy(mtest_buffer, test_buffer_);}

MH_medic::~MH_medic()
{
  delete [] mtest_buffer;
  delete [] TEST_p;
  delete [] TEST_r2;
  delete [] TEST_alpha;
  delete [] TEST_INTr2;
}

void MH_medic::test_u(event_record * rec_, int i_, bool verbose_)
{
  double t_history[3], net_r2_local;
  int f_local = start_test_u(rec_, t_history, net_r2_local),
      poffset=ndof*f_local,
      foffset=nbeads*f_local,
      *f_event=rec_->evframe_bead;

  double  *pref=xs+poffset,
          *r2_stable=rec_->r2stable_bead,
          *r2_unstable=rec_->r2unstable_bead,
          *alpha_bead=rec_->alpha_bead,
          net_r2_stable_local=net_r2_local,
          net_r2_unstable_local=0.0;
  do
  {
    f_local++; poffset+=ndof; pref+=ndof;
    advance((ts[f_local]-ts[f_local-1])/t_phys, d_ang[f_local-1], comega_s[f_local], dt_sim);
    t_history[2]=t_history[1]; t_history[1]=t_history[0]; t_history[0]=ts[f_local];
    for (int i = 0, j = 0; i < nbeads; i++,j+=2)
    {
      double  x_sim=q[i].x, y_sim=q[i].y,
              x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im,
              x_ref=pref[j], y_ref=pref[j+1],
              xerr=x_now-x_ref, yerr=y_now-y_ref, rsq=xerr*xerr+yerr*yerr;
      psim[j+poffset]=x_sim; psim[j+1+poffset]=y_sim;
      TEST_p[j+poffset]=x_now; TEST_p[j+1+poffset]=y_now;
      net_r2_local+=rsq;

      if (!(f_event[i]))
      {
        r2_state_comp[f_local][i]=rsq;
        INTr2_comp_history[i][2]=INTr2_comp_history[i][1]; INTr2_comp_history[i][1]=INTr2_comp_history[i][0];
        INTr2_comp_history[i][0]+=0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history[0]-t_history[1]);
        double alpha_it = alpha_state_comp[f_local-1][i]=alpha_comp(INTr2_comp_history[i], t_history[0], t_history[2]);

        TEST_r2[i+foffset]=r2_state_comp[f_local][i];
        TEST_INTr2[i+foffset]=INTr2_comp_history[i][0];
        TEST_alpha[i+foffset-nbeads]=alpha_it;

        if (alpha_it>alpha_tol)
        {
          f_event[i]=f_local-1;
          alpha_bead[i]=alpha_it;
          r2_unstable[i]=rsq;
          net_r2_unstable_local+=rsq;
        }
        else
        {
          r2_stable[i]+=rsq;
          net_r2_stable_local+=rsq;
        }
      }
      else
      {
        r2_unstable[i]+=rsq;
        net_r2_unstable_local+=rsq;

        TEST_r2[i+foffset]=rsq;
        TEST_INTr2[i+foffset]=0.5*(rsq+TEST_INTr2[i+foffset-nbeads])*(t_history[0]-t_history[1]);
        TEST_alpha[i+foffset-nbeads]=alpha_comp(TEST_INTr2[i+foffset], TEST_INTr2[i+foffset-(2*nbeads)], t_history[0], t_history[2]);
      }
    }
    if ((f_local+1)==Frames_test)
    {
      for (int i = 0; i < nbeads; i++) if (!(f_event[i]))
      {
        f_event[i] = Frames_test-1;
        alpha_bead[i] = alpha_state_comp[f_local-1][i];
      }
      break;
    }
  } while(true);

  net_r2=net_r2_local;
  net_r2_stable=net_r2_stable_local;
  net_r2_unstable=net_r2_unstable_local;

  stev_early=stev_late=f_event[0];
  nf_stable=0;
  for (int i = 0; i < nbeads; i++)
  {
    int fevent_it = f_event[i];
    nf_stable+=fevent_it;
    nev_state_comp[fevent_it][i]++;
    if (fevent_it<stev_early) stev_early = fevent_it;
    if (fevent_it>stev_late) stev_late = fevent_it;
  }
  nf_obs=stev_late*nbeads;
  nf_unstable=nf_obs-nf_stable;

  // set record data
  rec_->record_event_data(&nf_obs,&net_r2);
  if (verbose_) printf("(thread %d) Tested parameter %d.\n", thread_id,i_);
  write_utest_results(rec_,i_);
}

void MH_medic::write_utest_results(event_record *rec_, int i_)
{
  int header_len = 2;
  int header[] = {header_len, basic_tw_ilen, basic_tw_dlen};
  sprintf(mtest_buffer+buf_end, "par%d.redat", i_);
  FILE * data_file = fopen(mtest_buffer, "wb");
  fwrite(header,sizeof(int),header_len+1,data_file);
  fwrite(basic_tw_ints, sizeof(int), header[1],data_file);
  fwrite(basic_tw_dubs, sizeof(double), header[2],data_file);  
  fwrite(psim,sizeof(double),ndof*Frames,data_file);
  fwrite(r2_state_comp[0],sizeof(double),nbeads*Frames,data_file);
  fwrite(alpha_state_comp[0],sizeof(double),nbeads*Frames,data_file);
  fwrite(TEST_p,sizeof(double),ndof*Frames_test,data_file);
  fwrite(TEST_r2,sizeof(double),nbeads*Frames_test,data_file);
  fwrite(TEST_alpha,sizeof(double),nbeads*Frames_test,data_file);
  fwrite(TEST_INTr2,sizeof(double),nbeads*Frames_test,data_file);
  rec_->write_event_rec_full_header(data_file);
  rec_->write_record_data(data_file);
  fclose(data_file);
}

bool MH_medic::report_results(bool first2finish_, int ** nev_state_comp_)
{
  if (first2finish_)
    for (int i = 0; i < nbeads*Frames; i++) nev_state_comp_[0][i] = nev_state_comp[0][i];
  else
    for (int i = 0; i < nbeads*Frames; i++) nev_state_comp_[0][i] += nev_state_comp[0][i];
  return false;
}
