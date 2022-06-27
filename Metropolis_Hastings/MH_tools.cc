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

// event_block

event_block::event_block(int ncomp_, int nstates_): ncomp(ncomp_), nstates(nstates_),
stev_comp(new int[ncomp]), stev_ordered(new int[ncomp]), comps_ordered(new int[ncomp]),
nev_state_comp(Tmatrix<int>(nstates,ncomp)), nobs_state_comp(Tmatrix<int>(nstates,ncomp)),
rho2stable_comp(new double[ncomp]), delrho2_regime(new double[ncomp]),
mur2_state_comp(Tmatrix<double>(nstates,ncomp)), stdr2_state_comp(Tmatrix<double>(nstates,ncomp)),
mualpha_state_comp(Tmatrix<double>(nstates,ncomp)), stdalpha_state_comp(Tmatrix<double>(nstates,ncomp)) {}

event_block::~event_block()
{
  delete [] stev_comp; delete [] stev_ordered; delete [] comps_ordered;
  free_Tmatrix<int>(nev_state_comp); free_Tmatrix<int>(nobs_state_comp);
  delete[] rho2stable_comp; delete [] delrho2_regime;
  free_Tmatrix<double>(mur2_state_comp); free_Tmatrix<double>(stdr2_state_comp);
  free_Tmatrix<double>(mualpha_state_comp); free_Tmatrix<double>(stdalpha_state_comp);
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

void event_block::define_event_block(double sigma_scaled_,int dof_)
{
  double  sigma2_=sigma_scaled_*sigma_scaled_,
          dof_dub= (double)dof_;
  rho2stable=0.0;
  for (int i = 0; i < ncomp; i++)
    rho2stable+=rho2stable_comp[i]=sigma2_*dof_dub*((double)stev_comp[i]);
  delrho2_regime[0]=rho2stable;
  for (int i = 1; i < ncomp; i++)
    delrho2_regime[i]=sigma2_*dof_dub*((double)((i)*(stev_ordered[i]-stev_ordered[i-1])));
}

void event_block::synchronise_event_data(int stev_earliest_, int stev_latest_, double rho2stable_, int *stev_c_, int *stev_o_, int *comps_o_,double *rho2s_c_, double *drho2_r_)
{
  stev_earliest=stev_earliest_;
  stev_latest=stev_latest_;
  rho2stable=rho2stable_;
  for (int i = 0; i < ncomp; i++)
  {
    stev_comp[i]=stev_c_[i];
    stev_ordered[i]=stev_o_[i];
    comps_ordered[i]=comps_o_[i];
    rho2stable_comp[i]=rho2s_c_[i];
    delrho2_regime[i]=drho2_r_[i];
  }
}

// event_detector

event_detector::event_detector(int ncomp_, int nstates_, int dof_, double alpha_tol_): event_block(ncomp_,nstates_),
dof(dof_), ndof(ncomp*dof), alpha_tol(alpha_tol_),
r2_state_comp(Tmatrix<double>(nstates,ncomp)), alpha_state_comp(Tmatrix<double>(nstates,ncomp)),
INTr2_comp_history(Tmatrix<double>(ncomp,3)), r2_regime_comp(Tmatrix<double>(ncomp,ncomp)) {}

event_detector::~event_detector()
{
  free_Tmatrix<double>(r2_state_comp); free_Tmatrix<double>(alpha_state_comp);
  free_Tmatrix<double>(INTr2_comp_history); free_Tmatrix<double>(r2_regime_comp);
}
