#include "particle_race.hh"

extern "C"
{
  #include "AYaux.h"
}

void record::init(double * params_, int len_, int Frames_, double gau_h_, double gau_lambda_)
{
  params=params_;
  Frames=Frames_;
  len=len_;
  gau_h=gau_h_;
  gau_lambda=gau_lambda_;
}
void record::reset_record(int gen_)
{
  frscore=0;
  gen=gen_;
  parent_gen=-1; // negative value implies no parent
  parent_count=0;
  parent_global_index=-1; // negative value implies no parent
  dup_count=0;

  l2score=0;
}
void record::resample(int gen_, double * dmin_, double *dmax_, AYrng * r_)
{
  reset_record(gen_);
  for (int i = 0; i < len; i++)
    params[i] = dmin_[i]+(dmax_[i]-dmin_[i])*(r_->rand_uni_gsl(0.0,1.0));
}
void record::duplicate(record *parent_, int gen_, double *dmin_, double *dmax_, AYrng * r_)
{
  frscore=0;
  gen=gen_;
  parent_gen=parent_->gen;
  parent_count=parent_->parent_count+1;
  parent_global_index=parent_->global_index;
  dup_count=0;

  l2score=0;

  for (int i = 0; i < len; i++)
  {
    params[i] = parent_->params[i]*(r_->rand_gau_gsl(1.0, parent_->var()));
    if (params[i] > dmax_[i]) params[i]=dmax_[i];
    else if (params[i] < dmin_[i]) params[i]=dmin_[i];
  }
}
void record::take_vals(record * rtake_)
{
  frscore = rtake_->frscore;
  gen=rtake_->gen;
  parent_gen=rtake_->parent_gen;
  parent_count=rtake_->parent_count;
  parent_global_index=rtake_->parent_global_index;
  dup_count=rtake_->dup_count;

  l2score = rtake_->l2score;

  for (int i = 0; i < len; i++) params[i] = rtake_->params[i];
}





referee::~referee()
{
  if (alloc_flag)
  {
    delete [] leader_board;
    for (int i = 0; i < nlead+npool; i++) delete leaders[i];
    delete [] leaders;
    free_AYdmatrix(lead_params);
  }
}
void referee::alloc_records()
{
  leaders = new record*[nlead+npool];
  leader_board = new record*[nlead+npool];
  for (int i = 0; i < nlead+npool; i++) leaders[i] = new record(i);
  pool = leaders + nlead;
  pool_leaders = leader_board + nlead;

  lead_params = AYdmatrix(nlead+npool, param_len);
  pool_params = lead_params + nlead;

  alloc_flag = true;
}





int find_worst(record ** r, int ncap)
{
  int worst_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[worst_index]->isbetter(r[i]))
      worst_index = i;
  return worst_index;
}

int find_best(record ** r, int ncap)
{
  int best_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[best_index]->isworse(r[i]))
      best_index = i;
  return best_index;
}
