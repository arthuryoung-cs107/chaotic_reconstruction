#include "MH_tools.hh"

void fseek_SAFE(FILE *fp,long int offset,int origin)
{
  if(fseek(fp,offset,origin)!=0)
  {
    printf("fseek_SAFE: error shifting file position by %ld bytes\n",offset);
    exit(1);
  }
}

void fread_SAFE(void *ptr,size_t size,size_t count,FILE *fp)
{
  if(fread(ptr,size,count,fp)!=count)
  {
    printf("fread_SAFE: can't read file\n");
    exit(1);
  }
}

void event_block::synchronise_event_data(int stev_earliest_, int stev_latest_, int *stev_c_, int *stev_o_, int *comps_o_,double *rho2s_c_, double *drho2_r_)
{
  stev_earliest=stev_earliest_;
  stev_latest=stev_latest_;
  for (int i = 0; i < ncomp; i++)
  {
    stev_comp[i]=stev_c_[i];
    stev_ordered[i]=stev_o_[i];
    comps_ordered[i]=comps_o_[i];
    rho2stable_regime[i]=rho2s_c_[i];
    delrho2_regime[i]=drho2_r_[i];
  }
}

void event_block::define_event_block(double sigma_scaled_,int ndof_)
{
  double  sigma2_=sigma_scaled_*sigma_scaled_,
          dims_dub= (double)(ndof_*ncomp);
  for (int i = 0; i < ncomp; i++) rho2stable_comp[i]=sigma2_*dims_dub*((double)stev_comp[i]);
  delrho2_regime[0]=0.0; // for stable regime, no additional unstable residual added onto stable residual
  for (int i = 1; i < ncomp; i++) delrho2_regime[i]=sigma2_*dims_dub*((double)(stev_ordered[i]-stev_ordered[i-1]));
}

void event_block::clear_event_data()
{
  stev_latest=0;
  stev_earliest=nstates;
  for (int i=0,k=0; i < ncomp; i++)
  {
    stev_comp[i]=0;
    for (int j = 0; j < nstates; i++,k++)
    {
      nev_state_comp[0][k]=nobs_state_comp[0][k]=0;
      mur2_state_comp[0][k]=stdr2_state_comp[0][k]=
      mualpha_state_comp[0][k]=stdalpha_state_comp[0][k]=0.0;
    }
  }
}

void event_block::consolidate_event_data()
{
  // assumes that stev_comp, stev_ordered, and comps_ordered are already initialized.
  int early_index, early_frame, index_temp, frame_temp;
  early_frame = earliest(stev_ordered,early_index,0);
  index_temp=0; frame_temp=stev_ordered[index_temp];

  comps_ordered[index_temp]=early_index; stev_ordered[index_temp]=early_frame;
  comps_ordered[early_index]=index_temp; stev_ordered[early_index]=frame_temp;
  if (ncomp>1) earliest_recursive(stev_ordered,comps_ordered,1);
}
