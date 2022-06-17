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

bool event_block::report_event_data(bool first2finish_, int &stev_earliest_, int &stev_latest_, int *stev_c_, int ** nev_s_c_, int ** nobs_s_c_, double ** r2_s_c_, double **alpha_s_c_)
{
  for (int i = 0; i < ncomp*stev_latest; i++)
  {
    nev_s_c_[0][i] += nev_state_comp[0][i];
    nobs_s_c_[0][i] += nobs_state_comp[0][i];
    r2_s_c_[0][i] += mur2_state_comp[0][i];
    alpha_s_c_[0][i] += mualpha_state_comp[0][i];
  }
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
