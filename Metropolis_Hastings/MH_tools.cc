#include "MH_tools.hh"

void fseek_SAFE(FILE *fp,long int offset,int origin)
  {if(fseek(fp,offset,origin)!=0) {printf("fseek_SAFE: error shifting file position by %ld bytes\n",offset); exit(1);}}

void fread_SAFE(void *ptr,size_t size,size_t count,FILE *fp)
  {if (fread(ptr,size,count,fp)!=count) {printf("fread_SAFE: can't read file\n"); exit(1);}}

void print_row_vec(int *ivec_, int n_) {for (int i = 0; i < n_; i++) printf("%d ", ivec_[i]); printf("\n");}
void print_row_vec(double *dvec_, int n_) {for (int i = 0; i < n_; i++) printf("%e ", dvec_[i]); printf("\n");}

// event_block

event_block::event_block(int ncomp_, int nstates_, int dim_): ncomp(ncomp_), nstates(nstates_), dim(dim_), dof(ncomp*dim),
stev_comp(new int[ncomp]), stev_ordered(new int[ncomp]), comps_ordered(new int[ncomp]),
nev_state_comp(Tmatrix<int>(nstates,ncomp)), nobs_state_comp(Tmatrix<int>(nstates,ncomp)),
rho2stable_comp(new double[ncomp]), rho2unstable_comp(new double[ncomp]),
mur2_state_comp(Tmatrix<double>(nstates,ncomp)), stdr2_state_comp(Tmatrix<double>(nstates,ncomp)),
mualpha_state_comp(Tmatrix<double>(nstates,ncomp)), stdalpha_state_comp(Tmatrix<double>(nstates,ncomp)) {}

event_block::~event_block()
{
  delete [] stev_comp; delete [] stev_ordered; delete [] comps_ordered;
  free_Tmatrix<int>(nev_state_comp); free_Tmatrix<int>(nobs_state_comp);
  delete [] rho2stable_comp; delete [] rho2unstable_comp;
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
    for (int j = 0; j < nstates; j++,k++)
    {
      nev_state_comp[0][k]=nobs_state_comp[0][k]=0;
      mur2_state_comp[0][k]=stdr2_state_comp[0][k]=
      mualpha_state_comp[0][k]=stdalpha_state_comp[0][k]=0.0;
    }
  }
}

void event_block::consolidate_event_data()
{
  // assumes that stev_ordered and comps_ordered are already initialized to stev_comp order
  int early_index,
      early_frame=earliest(stev_comp,early_index,0);

  stev_ordered[early_index]=stev_ordered[0]; stev_ordered[0]=early_frame;
  comps_ordered[0]=early_index; comps_ordered[early_index]=0;

  if (ncomp>1) earliest_recursive(stev_ordered,comps_ordered,1);
}

void event_block::define_event_block(double sigma_scaled_)
{
  double  obs_scale=(sigma_scaled_*((double)dim))*sigma_scaled_,
          netrho2_comp_full=obs_scale*((double)stev_ordered[ncomp-1]);
  for (int i = 0; i < ncomp; i++)
  {
    rho2stable_comp[i]=obs_scale*((double)stev_comp[i]);
    rho2unstable_comp[i]=netrho2_comp_full-rho2stable_comp[i];
  }
}

void event_block::synchronise_event_data(int stev_e_, int stev_l_, int *stev_c_, int *stev_o_, int *comps_o_,double *rho2s_c_, double *rho2us_c_)
{
  stev_earliest=stev_e_;
  stev_latest=stev_l_;
  for (int i = 0; i < ncomp; i++)
  {
    stev_comp[i]=stev_c_[i];
    stev_ordered[i]=stev_o_[i];
    comps_ordered[i]=comps_o_[i];
    rho2stable_comp[i]=rho2s_c_[i];
    rho2unstable_comp[i]=rho2us_c_[i];
  }
}
