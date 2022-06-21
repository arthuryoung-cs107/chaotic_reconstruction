#include "MH_solvers.hh"

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

void MH_doctor::stage_diagnostics()
{
  for (int frame_i = 0; frame_i < Frames_test; frame_i++)
  {
    int foffset = 2*nbeads*frame_i;
    double *f = xs+(foffset);
    for (int bead_i = 0, j=0; bead_i < nbeads; bead_i++, j+=2)
    {TEST_refp[j+foffset]=f[j]; TEST_refp[j+1+foffset]=f[j+1];}
  }
  int MH_doc_header_len = 1;
  int header_ints[] = {MH_doc_header_len, Frames_test};
  sprintf(test_buffer+test_buf_end, "results_startspecs.redat");
  FILE * test_startspecs = fopen(test_buffer, "wb");
  write_MH_params(test_startspecs);
  fwrite(header_ints, sizeof(int), MH_doc_header_len+1, test_startspecs);
  fwrite(TEST_refp, sizeof(double), 2*nbeads*Frames_test, test_startspecs);
  fclose(test_startspecs);
}
void MH_doctor::close_diagnostics()
{
  int MH_doc_header_len = 0;
  int header_ints[] = {MH_doc_header_len};
  sprintf(test_buffer+test_buf_end, "results_endspecs.redat");
  FILE * test_endspecs = fopen(test_buffer, "wb");
  fwrite(header_ints, sizeof(int), MH_doc_header_len+1, test_endspecs);
  fwrite(evcount_bead_frame[0], sizeof(int), 2*nbeads*Frames_test, test_endspecs);
  fclose(test_endspecs);
}
