#include "particle_relay.hh"

extern "C"
{
  #include "AYaux.h"
}

void record::record_event_data(double res_acc_, int *kill_frames_, double ** INTpos_res_, double *alpha_kill_)
{
  residual = res_acc_; event_residual = 0.0;
  for (int i = 0; i < beads; i++)
  {
    event_positions[i] = kill_frames_[i];
    residual_data[i] = INTpos_res_[i][2];
    event_residual += INTpos_res_[i][0];
    alpha_data[i] = alpha_kill_[i];
  }
}
void record::reset_record(int gen_, int p_gen_, int p_count_, int p_gi_)
{
  gen=gen_;
  parent_gen=p_gen_; // negative value implies no parent
  parent_count=p_count_;
  parent_global_index=p_gi_; // negative value implies no parent
  dup_count=0;
  for (int i = 0; i < record_double_len; i++) double_params[i] = 0.0;
}
void record::resample(int gen_, double * dmin_, double *dmax_, AYrng * r_)
{
  reset_record(gen_);
  for (int i = 0; i < len; i++)
    params[i] = dmin_[i]+(dmax_[i]-dmin_[i])*(r_->rand_uni_gsl(0.0,1.0));
}
void record::duplicate(record *parent_, int gen_, double *dmin_, double *dmax_, AYrng * r_, double * var_)
{
  reset_record(gen_, parent_->gen, parent_->parent_count+1, parent_->global_index);
  double sfac = 0.25*(1.0-parent_->px);
  for (int i = 0; i < len; i++)
  {
    double z = r_->rand_gau_gsl(0.0, 1.0);
    params[i] = parent_->params[i] + sfac*z*((z>0.0)?(dmax_[i]-parent_->params[i]):(parent_->params[i]-dmin_[i]));
    if (params[i] > dmax_[i]) params[i]=dmax_[i];
    else if (params[i] < dmin_[i]) params[i]=dmin_[i];
  }
}
void record::take_vals(record * rtake_)
{
  for (int i = 0; i < record_int_len-1; i++) int_params[i] = rtake_->int_params[i];
  for (int i = 0; i < record_double_len; i++) double_params[i] = rtake_->double_params[i];
  for (int i = 0; i < record_int_chunk_count*beads; i++) int_chunk[i]=rtake_->int_chunk[i];
  for (int i = 0; i < record_double_chunk_count*beads; i++) double_chunk[i]=rtake_->double_chunk[i];
  for (int i = 0; i < len; i++) params[i] = rtake_->params[i];
}

void events::define_relay_leg(int leg_id_)
{
  printf("(leg %d): Events identified. Event frames -", leg_id_);

  // Update the event frames for each bead. Also identify the earliest and latest events.
  int latest_event = 0, earliest_event=Frames;
  for (int i = 0; i < beads; i++) for (int j = 0; j < Frames; j++)
    if (global_event_frame_count[i][j])
    {
      event_frames[i] = j-1;
      if (event_frames[i]>latest_event) latest_event = event_frames[i];
      if (event_frames[i]<earliest_event) earliest_event = event_frames[i];
      printf(" %d", j-1);
      break;
    }
  printf(". Earliest: %d, latest: %d ", earliest_event, latest_event);
  *ref_earliest_event=earliest_event; *ref_latest_event=latest_event;

  // sort the events
  int i_it, early_it, i_temp, early_temp;
  for (int i = 0; i < beads; i++)
  {
    event_sorted[i] = event_frames[i];
    bead_order[i] = i;
  }
  early_it = earliest(&i_it, 0);
  early_temp = event_sorted[0]; i_temp = 0;
  event_sorted[0]=early_it; bead_order[0]=i_it;
  event_sorted[i_it]=early_temp; bead_order[i_it]=i_temp;
  earliest_recursive(1);

  printf("Events, ordered:\n");
  for (int i = 0; i < beads; i++)
    printf("(event %d) bead %d - frame %d\n",i,bead_order[i],event_sorted[i]);

  tau = (noise_tol*data_scale)*sqrt((double)(event_observations));
  tau_sqr = tau*tau;
}

referee::~referee()
{
  if (alloc_flag)
  {
    free_AYimatrix(record_int_chunk);
    free_AYimatrix(global_event_frame_count);
    free_AYdmatrix(param_chunk);
    free_AYdmatrix(record_double_chunk);

    delete [] leader_board;
    for (int i = 0; i < nlead+npool; i++) delete records[i];
    delete [] records;

    delete [] lead_dup_count;
    delete [] event_frames;
    delete [] ev;

    delete [] sample_weights;
    delete [] lead_par_w_mean;
    delete [] lead_par_w_var;
  }
}

void referee::alloc_records(int nt_, int Frames_, int beads_)
{
  record_int_chunk = AYimatrix(nlead+npool, record_int_chunk_count*beads_);
  global_event_frame_count = AYimatrix(beads_,Frames_);
  record_double_chunk = AYdmatrix(nlead+npool, record_double_chunk_count*beads_);
  param_chunk = AYdmatrix(nlead+npool, param_len);

  records = leaders = new record*[nlead+npool];
  leader_board = new record*[nlead+npool];

  for (int i = 0; i < nlead+npool; i++)
    records[i] = new record(beads_,Frames_,param_len,i,record_int_chunk[i],param_chunk[i],record_double_chunk[i]);

  pool = records + nlead;
  candidates = leader_board + nlead;

  lead_dup_count = new int[nlead];
  event_frames = new int[beads_];
  ev = new events(beads_, event_frames);


  sample_weights = new double[nlead];
  lead_par_w_mean = new double[param_len];
  lead_par_w_var = new double[param_len];

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
