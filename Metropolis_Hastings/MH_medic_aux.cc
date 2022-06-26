#include "MH_workers.hh"

// MH_medic

MH_medic::MH_medic(MH_examiner * ex_, int Frames_test_, char * test_buffer_, size_t test_buf_end_): ex(ex_),
buf_end(test_buf_end_), Frames_test(Frames_test_),
mtest_buffer(new char[buf_end+50]),
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

int MH_medic::start_test_u(event_record * rec_, double * t_history_, double &net_r2_local_)
{
  for (int i = 0; i < nbeads*Frames; i++) r2_state_comp[0][i] = alpha_state_comp[0][i] = 0.0;

  reset_sim(rec_->u, ts[0]/t_phys, d_ang[0], comega_s[0], xs);

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
    f_event[i]=0;
    r2_stable[i]=r2_unstable[i]=r2_state_comp[0][i]=
    INTr2_comp_history[i][0]=INTr2_comp_history[i][1]=INTr2_comp_history[i][2]=0.0;
    alpha_state_comp[0][i]=NAN;

    double  x_sim=q[i].x, y_sim=q[i].y,
            x_now=(x_sim-cx)*cl_im+cx_im, y_now=(y_sim-cy)*cl_im+cy_im;
    psim[j]=x_sim; psim[j+1]=y_sim;

    // set medic data
    TEST_p[j]=x_now; TEST_p[j+1]=y_now;
    TEST_r2[i]=r2_state_comp[0][i];
    TEST_INTr2[i]=INTr2_comp_history[i][0];
    TEST_alpha[i]=alpha_state_comp[0][i];
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
    psim[j+poffset]=x_sim; psim[j+1+poffset]=y_sim;

    net_r2_local+=r2_state_comp[f_local][i]=rsq;
    r2_stable[i]+=rsq;

    INTr2_comp_history[i][2]=INTr2_comp_history[i][1]; INTr2_comp_history[i][1]=INTr2_comp_history[i][0];
    INTr2_comp_history[i][0]+=0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history_[0]-t_history_[1]);
    alpha_state_comp[f_local][i]=NAN;

    // set medic data
    TEST_p[j+poffset]=x_now;TEST_p[j+1+poffset]=y_now;
    TEST_r2[i+foffset]=r2_state_comp[f_local][i];TEST_INTr2[i+foffset]=INTr2_comp_history[i][0];
    TEST_alpha[i+foffset]=alpha_state_comp[f_local][i];
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
    psim[j+poffset]=x_sim; psim[j+1+poffset]=y_sim;

    net_r2_local+=r2_state_comp[f_local][i]=rsq;
    r2_stable[i]+=rsq;
    INTr2_comp_history[i][2]=INTr2_comp_history[i][1]; INTr2_comp_history[i][1]=INTr2_comp_history[i][0];
    INTr2_comp_history[i][0]+=0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history_[0]-t_history_[1]);

    // set medic data
    TEST_p[j+poffset]=x_now;TEST_p[j+1+poffset]=y_now;
    TEST_r2[i+foffset]=r2_state_comp[f_local][i];TEST_INTr2[i+foffset]=INTr2_comp_history[i][0];
  }
  net_r2_local_=net_r2_local;
  return f_local;
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
      double rsq = compare_p(q[i].x, q[i].y,pref[j],pref[j+1]);
      
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
  int header[] = {header_len, int_len, double_len};
  sprintf(mtest_buffer+buf_end, "par%d.redat", i_);
  FILE * data_file = fopen(mtest_buffer, "wb");
  fwrite(header,sizeof(int),header_len+1,data_file);
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
