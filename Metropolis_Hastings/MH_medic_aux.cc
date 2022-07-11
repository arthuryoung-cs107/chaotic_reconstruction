#include "MH_workers.hh"

// MH_medic

MH_medic::MH_medic(MH_examiner &ex_, int Frames_test_, char * test_buffer_, size_t test_buf_end_): ex(ex_),
buf_end(test_buf_end_), Frames_test(Frames_test_),
mtest_buffer(new char[buf_end+50]),
TEST_p(new double[dof*Frames_test]), TEST_INTr2(new double[nbeads*Frames_test])
{strcpy(mtest_buffer, test_buffer_);}

MH_medic::~MH_medic()
{
  delete [] mtest_buffer;
  delete [] TEST_p;
  delete [] TEST_INTr2;
}

void MH_medic::start_test_u(int &f_local_,int *f_event_,double &netr2_local_,double &netr2_stable_local_,double &netr2_unstable_local_,double *t_history_,double *r2net_bead_,double *r2stable_bead_,double *r2unstable_bead_,double *alphaev_bead_)
{
  // initialize everything
  f_local_=0;
  netr2_local_=netr2_stable_local_=netr2_unstable_local_=0.0;
  t_history_[0]=t_history_[1]=t_history_[2]=ts[f_local_];
  for (int i = 0; i < nbeads; i++)
  {
    r2_state_comp[f_local_][i]=INTr2_comp_history[i][0]=INTr2_comp_history[i][1]=
    r2net_bead_[i]=r2stable_bead_[i]=r2unstable_bead_[i]=0.0;
    alpha_state_comp[f_local_][i]=alphaev_bead_[i]=NAN;
    f_event_[i]=0;

    // medic specific data
    TEST_p[i*dim]=xs[i*dim]; TEST_p[(i*dim)+1]=xs[(i*dim)+1];
    TEST_INTr2[i]=0.0;
  }

  // resolve frames 1 and 2
  double *pref;
  for (int iframe = 1; iframe <= 2; iframe++)
  {
    pref=advance_sim(++f_local_,t_history_);
    for (int i=0,j=0,k=f_local_*dof,l=f_local_*nbeads; i < nbeads; i++,j+=dim,k+=dim,l++)
    {
      double rsq=compute_residual(psim[k]=q[i].x,psim[k+1]=q[i].y,TEST_p[k],TEST_p[k+1],pref[j],pref[j+1]);
      r2_state_comp[f_local_][i]=rsq;
      netr2_local_+=rsq; r2net_bead_[i]+=rsq;

      update_integral_history(0.5*(rsq+r2_state_comp[f_local_-1][i])*(t_history_[0]-t_history_[1]),i);

      alpha_state_comp[f_local_][i]=NAN;

      r2stable_bead_[i]+=rsq;

      // medic specific data
      TEST_INTr2[l]=INTr2_comp_history[i][0];
    }
  }
  // set stable residual, assuming that frames 0,1,2 were all stable
  netr2_stable_local_=netr2_local_;
}

void MH_medic::test_u(event_record * rec_, int i_, double *r2i_, double *alphai_, bool verbose_)
{
  reset_sim(rec_->u,ts[0]/t_phys, d_ang[0], comega_s[0], xs);

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

  start_test_u(f_local,f_event,netr2_local,netr2_stable_local,netr2_unstable_local,t_history,r2net_bead,r2stable_bead,r2unstable_bead,alphaev_bead);

  do
  {
    pref=advance_sim(++f_local,t_history);
    for (int i=0,j=0,k=f_local*dof,l=f_local*nbeads; i < nbeads; i++,j+=dim,k+=dim,l++)
    {
      double rsq = compute_residual(psim[k]=q[i].x,psim[k+1]=q[i].y,TEST_p[k],TEST_p[k+1],pref[j],pref[j+1]);
      r2_state_comp[f_local][i]=rsq;
      netr2_local+=rsq; r2net_bead[i]+=rsq;

      update_integral_history(0.5*(rsq+r2_state_comp[f_local-1][i])*(t_history[0]-t_history[1]),i);

      double alpha_it=alpha_state_comp[f_local-1][i]=alpha_comp(INTr2_comp_history[i],t_history[0],t_history[2]);

      // medic specific data
      TEST_INTr2[l]=INTr2_comp_history[i][0];

      if (!(f_event[i]))
      {
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
    if ((f_local+1)==Frames_test)
    {
      for (int i = 0; i < nbeads; i++) if (!(f_event[i]))
      {
        f_event[i] = Frames_test-1;
        alphaev_bead[i] = alpha_state_comp[f_local-1][i];
      }
      break;
    }
  } while(true);
  netr2=netr2_local,
  netr2_stable=netr2_stable_local,
  netr2_unstable=netr2_unstable_local;

  update_event_data(f_local,f_event,r2i_,alphai_);
  rec_->record_event_data(&nf_netobs,&netr2);

  if (verbose_) printf("(thread %d) Tested parameter %d.\n", thread_id,i_);
  write_utest_results(rec_,i_);
}

void MH_medic::write_utest_results(event_record *rec_, int i_)
{
  sprintf(mtest_buffer+buf_end, "par%d.mhdat", i_);
  FILE * data_file = fopen(mtest_buffer, "wb");
  fwrite(psim,sizeof(double),dof*Frames_test,data_file);
  fwrite(r2_state_comp[0],sizeof(double),nbeads*Frames_test,data_file);
  fwrite(alpha_state_comp[0],sizeof(double),nbeads*Frames_test,data_file);
  fwrite(TEST_p,sizeof(double),dof*Frames_test,data_file);
  fwrite(TEST_INTr2,sizeof(double),nbeads*Frames_test,data_file);
  rec_->write_event_rec_full_header(data_file,1);
  rec_->write_record_data(data_file);
  fclose(data_file);
}
