#include "MH_solvers.hh"

// MH_genetic
void MH_genetic::stage_diagnostics()

{
  int header_len = 2;
  int header[] = {header_len, genetic_train_const_ilen, genetic_train_const_dlen};
  sprintf(obuf+obuf_end, "startspecs.mhdat");
  FILE * startspecs_file = fopen(obuf, "wb");
  write_MH_params(startspecs_file);
  fwrite(header_ints, sizeof(int), header_len+1, startspecs_file);
  fwrite(genetic_train_const_ints, sizeof(int), genetic_train_const_ilen, startspecs_file);
  fwrite(genetic_train_const_dubs, sizeof(double), genetic_train_const_dubs, startspecs_file);
  fclose(startspecs_file);
}

void MH_genetic::clear_genetic_event_data()
{
  event_block::clear_event_data();
  #pragma omp parallel
  {
    double * clear_buf;
    #pragma omp for nowait
    for (int i = 0; i < npool; i++)
    {clear_buf=r2_pool_Framebead[i]; for (int j = 0; j < nbeads*Frames; j++) clear_buf[j]=0.0;}

    #pragma omp for nowait
    for (int i = 0; i < npool; i++)
    {clear_buf=alpha_pool_Framebead[i]; for (int j = 0; j < nbeads*Frames; j++) clear_buf[j]=0.0;}
  }
}

void MH_genetic::consolidate_genetic_event_data()
{
  for (int i = 0; i < nbeads; i++)
  {
    comps_ordered[i]=i;
    stev_ordered[i]=stev_comp[i];
  }
  event_block::consolidate_event_data();
}

void MH_genetic::consolidate_genetic_event_data(int bleader_index_)
{
  pool[bleader_index_]->determine_event_block(stev_earliest,stev_latest,stev_comp,stev_ordered,comps_ordered);
  event_block::consolidate_event_data();
}

void MH_genetic::synchronise_genetic_event_data()
{
  int nf_obs, nf_stable, nf_unstable;
  event_block::set_state_counts(nf_obs, nf_stable, nf_unstable);
  #pragma omp parallel
  {
    MH_examiner *ex_t=examiners[thread_num()];
    ex_t->synchronise_examiner_event_data(&nf_obs,stev_earliest,stev_latest,stev_comp,stev_ordered,comps_ordered,rho2stable_comp,delrho2_regime);
  }
}

void MH_genetic::report_genetic_event_data(event_record **recs_, int n_)
{
  // begin by finishing computation of event detection statistics
  #pragma omp parallel
  {
    // finish computing mean residuals for each component of the state space
    #pragma omp for
    for (int i = 0; i < nbeads*stev_latest; i++)
    {
      double  nobs_inv = 1.0/((double)(nobs_state_comp[0][i])),
              nobsm1_inv = 1.0/((double)(nobs_state_comp[0][i]-1)),
              mur2_i = mur2_state_comp[0][i]*=nobs_inv,
              mualpha_i = mualpha_state_comp[0][i]*=nobs_inv,
              var_r2=0.0,
              var_alpha=0.0;

      for (int j = 0; j < n_; j++) if ((recs_[j]->nf_obs/nbeads)>=(i/nbeads))
      {
        double  diffr2 = r2_pool_Framebead[j][i]-mur2_i,
                diffalpha = alpha_pool_Framebead[j][i]-mualpha_i;
        var_r2+=diffr2*diffr2;
        var_alpha+=diffalpha*diffalpha;
      }
      stdr2_state_comp[0][i]=sqrt(var_r2*nobsm1_inv);
      stdalpha_state_comp[0][i]=sqrt(var_alpha*nobsm1_inv);
    }
  }

  // write event data
  int hlen=2;
  int header[] = {hlen, n_, stev_latest};
  sprintf(obuf+obuf_end, "event_block%d.mhdat",event_block_count);
  FILE * data_file = fopen(obuf, "wb");
  fwrite(stev_comp, sizeof(int), nbeads, data_file);
  fwrite(stev_ordered, sizeof(int), nbeads, data_file);
  fwrite(comps_ordered, sizeof(int), nbeads, data_file);
  fwrite(nev_state_comp, sizeof(int), nbeads*stev_latest, data_file);
  fwrite(nobs_state_comp, sizeof(int), nbeads*stev_latest, data_file);
  fwrite(rho2stable_comp, sizeof(double), nbeads, data_file);
  fwrite(delrho2_regime, sizeof(double), nbeads, data_file);
  fwrite(mur2_state_comp, sizeof(double), nbeads*stev_latest, data_file);
  fwrite(stdr2_state_comp, sizeof(double), nbeads*stev_latest, data_file);
  fwrite(mualpha_state_comp, sizeof(double), nbeads*stev_latest, data_file);
  fwrite(stdalpha_state_comp, sizeof(double), nbeads*stev_latest, data_file);
  recs_[0]->write_event_record_full_header(data_file,n_);
  for (int i = 0; i < n_; i++) recs_[i]->write_record_data(data_file);
  fclose(data_file);
  event_block_count++;
}

void MH_genetic::post_event_resampling(event_record ** recs_, int n_)
{
  // start by using the r2_pool_Framebead data to adjust the records back to the correct event data
  double r2_min=DBL_MAX;
  #pragma omp parallel
  {
    MH_examiner *ex_t = examiners[thread_num()];
    #pragma omp for reduction(min:r2_min)
    for (int i = 0; i < n_; i++)
    {
      ex_t->restore_event_record(recs_[i],r2_pool_Framebead[i]);
      double r2_it = recs_[i]->set_stable_objective();
      if (r2_it<r2_min) r2_min=r2_it;
      leader_board[i]=recs_[i];
    }
  }
  if (n_!=npool) printf("(MH_genetic::post_event_resampling) Warning: n_!=npool\n");

  // collect leaders
  MH_trainer::pick_nbest_records(leader_board,nlead,n_);
  set_leader_records();
  write_Class_diagnostics();
  Class_count++;

  // resample pool
  double w_sum = basic_MH_trainer::compute_weights(r2_min,rho2,leaders,nlead);
  basic_MH_trainer::respawn_pool(w_sum,examiners,pool,leaders);
  gen_count++;
}

void MH_genetic::set_leader_records()
{
  nreplace=MH_trainer::take_records(leader_board,leaders,irepl_leaders,nlead);
  for (int i = 0; i < nreplace; i++)
  {
    int repl_index = irepl_leaders[i];
    // make sure that first bit of leaderboard is always pointing to leaders
    leader_board[repl_index]=leaders[repl_index];
    leaders[repl_index]->init_basic_record(gen_count,Class_count);
  }

  leader_count=nlead;
  bleader_rid=wleader_rid=0;
  for (int i = 1; i < nlead; i++)
  {
    if (leaders[bleader_rid]->isworse(leaders[i])) bleader_rid=i; // leader i is better than the current best record
    if (leaders[wleader_rid]->isbetter(leaders[i])) wleader_rid=i; // leader i is worse than the current worst record
  }
  br2=leaders[bleader_rid]->r2compare;
  wr2=leaders[wleader_rid]->r2compare;
}

void MH_genetic::consolidate_genetic_training_data()
{
  ncandidates=nsuccess;
  for (int i = 0; i < nsuccess; i++) candidates[i]=pool[isuccess_pool[i]];
  // if we have more successful particles than we can store, we have to narrow down the candidates
  if (nsuccess>nlead) MH_trainer::pick_nbest_records(candidates,ncandidates=nlead,nsuccess);
  MH_trainer::pick_nbest_records(leader_board,nlead,leader_count+ncandidates);
  set_leader_records();
  prob_best=gaussian_likelihood::compute_prob(sqrt(br2),sqrt(rho2));
  prob_worst=gaussian_likelihood::compute_prob(sqrt(wr2),sqrt(rho2));
}

void MH_genetic::report_genetic_training_data()
{
  if (nreplace)
  {
    int hlen_Class=2;
    int header_Class[] = {hlen_Class, genetic_it_ilen_full(), genetic_it_dlen_full()};
    sprintf(obuf+obuf_end, "Class%d.mhdat",Class_count);
    FILE * Class_file = fopen(obuf, "wb");
    fwrite(header_Class,sizeof(int),hlen_Class+1,Class_file);
    write_genetic_it_ints(Class_file);
    write_genetic_it_dubs(Class_file);
    leaders[0]->write_event_rec_training_header(Class_file,nlead);
    for (int i = 0; i < nlead; i++) leaders[i]->write_event_rec_training_data(Class_file);
    fclose(Class_file);
    Class_count++;
  }
  int hlen_gen=2;
  int header_gen[] = {hlen_gen, genetic_it_ilen_full(), genetic_it_dlen_full()};
  sprintf(obuf+obuf_end, "gen%d.mhdat",gen_count);
  FILE * gen_file = fopen(obuf, "wb");
  write_genetic_it_ints(gen_file);
  write_genetic_it_dubs(gen_file);
  basic_MH_trainer::write_ustats(gen_file);
  fclose(gen_file);
  gen_count++;
}

void MH_genetic::close_diagnostics()
{
  int header_len = 2;
  int header[] = {header_len, genetic_it_ilen_full(), genetic_it_dlen_full()};
  sprintf(obuf+obuf_end, "endspecs.mhdat");
  FILE * endspecs_file = fopen(obuf, "wb");
  fwrite(header_ints, sizeof(int), header_len+1, endspecs_file);
  write_genetic_it_ints(endspecs_file);
  write_genetic_it_dubs(endspecs_file);
  records[0]->write_event_rec_full_header(endspecs_file, nlead+npool);
  for (int i = 0; i < nlead+npool; i++) records[i]->write_record_data(endspecs_file);
  fclose(endspecs_file);
}
