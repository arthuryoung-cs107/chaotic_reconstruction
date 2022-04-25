#include "particle_walk.hh"

extern "C"
{
  #include "AYaux.h"
}

void record::record_event_data(int *kill_frames_, double ** INTpos_res_, double *alpha_kill_)
{
  residual = event_residual = 0.0;
  for (int i = 0; i < beads; i++)
  {
    event_positions[i] = kill_frames_[i];
    residual += residual_data[i] = INTpos_res_[i][2];
    event_residual += INTpos_res_[i][0];
    alpha_data[i] = alpha_kill_[i];
  }
}
void record::reset_record(int gen_)
{
  frscore=0;
  gen=gen_;
  parent_gen=-1; // negative value implies no parent
  parent_count=0;
  parent_global_index=-1; // negative value implies no parent
  dup_count=0;

  l2score=0.0;
}
void record::resample(int gen_, double * dmin_, double *dmax_, AYrng * r_)
{
  reset_record(gen_);
  for (int i = 0; i < len; i++)
    params[i] = dmin_[i]+(dmax_[i]-dmin_[i])*(r_->rand_uni_gsl(0.0,1.0));
}
void record::duplicate(record *parent_, int gen_, double *dmin_, double *dmax_, AYrng * r_, double * var_)
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
    params[i] = parent_->params[i]*(r_->rand_gau_gsl(1.0, var_[i])); // do we account for the displacement relative to mean?
    if (params[i] > dmax_[i]) params[i]=dmax_[i];
    else if (params[i] < dmin_[i]) params[i]=dmin_[i];
  }
}
void record::take_vals(record * gtake_)
{
  frscore=gtake_->frscore;
  gen=gtake_->gen;
  parent_gen=gtake_->parent_gen;
  parent_count=gtake_->parent_count;
  parent_global_index=gtake_->parent_global_index;
  dup_count=gtake_->dup_count;

  l2score = gtake_->l2score;

  for (int i = 0; i < len; i++) params[i] = gtake_->params[i];
}




referee::~referee()
{
  if (alloc_flag)
  {
    free_AYdmatrix(param_chunk);

    delete [] leader_board;
    for (int i = 0; i < nlead+npool; i++) delete records[i];
    delete [] records;

    delete [] lead_dup_count; delete [] frame_kill_count;

    delete [] sample_weights;
    delete [] gen_frame_res_data;
    delete [] gen_param_mean;
    delete [] gen_param_var;
  }
}
void referee::alloc_records(int nt_, int Frames_)
{
  param_chunk = AYdmatrix(nlead+npool, param_len);

  records = leaders = classA = new record*[nlead+npool];
  leader_board = new record*[nlead+npool];

  for (int i = 0; i < nlead+npool; i++) records[i] = new record(i,param_chunk[i]);
  pool = records + nlead;
  classB = records + nA;
  candidates = leader_board + nlead;

  lead_dup_count = new int[nlead];
  frame_kill_count = new int [Frames_];

  sample_weights = new double[nlead];
  gen_frame_res_data = new double[4*Frames_];
  gen_param_mean = new double[param_len];
  gen_param_var = new double[param_len];

  alloc_flag = true;
}

int find_worst_record(record ** r, int ncap)
{
  int worst_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[worst_index]->isbetter(r[i]))
      worst_index = i;
  return worst_index;
}

int find_best_record(record ** r, int ncap)
{
  int best_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[best_index]->isworse(r[i]))
      best_index = i;
  return best_index;
}
