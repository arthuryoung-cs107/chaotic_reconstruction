#include "MH_solvers.hh"

// MH_genetic

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
  int nf_obs, nf_stable, nf_regime, nf_unstable;
  event_block::set_state_counts(nf_obs, nf_stable, nf_regime, nf_unstable);
  #pragma omp parallel
  {
    examiners[thread_num()]->synchronise_examiner_event_data(&nf_obs,stev_earliest,stev_latest,rho2stable,stev_comp,stev_ordered,comps_ordered,rho2stable_comp,delrho2_regime);
  }
}

void MH_genetic::report_genetic_event_data()
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

      for (int j = 0; j < npool; j++) if ((pool[j]->nfobs/nbeads)>=(i/nbeads))
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
  write_event_diagnostics(event_block_count);
}

void MH_genetic::write_event_diagnostics(int &event_block_count_)
{
  // write event data
  int hlen=2;
  int header[] = {hlen, npool, stev_latest};
  sprintf(obuf+obuf_end, "event_block%d.mhdat",event_block_count_);
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
  pool[0]->write_event_rec_full_header(data_file,npool);
  for (int i = 0; i < npool; i++) pool[i]->write_record_data(data_file);
  fclose(data_file);
  event_block_count_++;
}

void MH_genetic::set_stable_objective()
{
  rho2=rho2stable;
  // start by using the r2_pool_Framebead data to adjust the records back to the correct event data
  double r2_min=DBL_MAX;
  #pragma omp parallel
  {
    MH_examiner *ex_t = examiners[thread_num()];
    #pragma omp for reduction(min:r2_min) nowait
    for (int i = 0; i < npool; i++)
    {
      ex_t->restore_event_record(pool[i],r2_pool_Framebead[i]);
      double r2_it = pool[i]->set_record_stable();
      if (r2_it<r2_min) r2_min=r2_it;
      leader_board[i]=pool[i];
    }
    #pragma omp for
    for (int i = 0; i < nlead; i++)
    {
      double r2_work = ex_t->set_record_stable();
    }
  }
  double w_sum_full = compute_weights(r2_min,rho2stable,pool,npool);
  compute_weighted_ustats(w_sum_full,pool,npool);

  // collect leaders
  pick_nbest_records(leader_board,nlead,npool);
  set_leader_records();
  report_genetic_training_data();

  // resample pool
  double w_sum=0.0; // given that weights have already been computed, just sum up leader terms

  #pragma omp for reduction(+:w_sum)
  for (int i = 0; i < nlead; i++)
    w_sum+=leaders[i]->w;

  respawn_pool(w_sum);
}

void MH_genetic::set_unstable_objective()
{
  double  rho2full = ((double)2*stev_latest*ncomp)*(sigma_scaled*sigma_scaled),
          rho2unstable=rho2=rho2full-rho2stable; // rho2unstable

  double r2_min=DBL_MAX;
  #pragma omp parallel
  {
    MH_examiner *ex_t = examiners[thread_num()];
    #pragma omp for reduction(min:r2_min) nowait
    for (int i = 0; i < npool+nlead; i++)
    {
      double r2_it = records[i]->set_record_unstable();
      if (r2_it<r2_min) r2_min=r2_it;
      leader_board[i]=records[i];
    }
  }
  double w_sum_full = compute_weights(r2_min,rho2unstable,records,npool+nlead);
  compute_weighted_ustats(w_sum_full,records,npool+nlead);

  // collect leaders
  pick_nbest_records(leader_board,nlead,npool+nlead);
  set_leader_records();
  report_genetic_training_data();

  // resample pool
  double w_sum=0.0; // given that weights have already been computed, just sum up leader terms

  #pragma omp for reduction(+:w_sum)
  for (int i = 0; i < nlead; i++)
    w_sum+=leaders[i]->w;

  respawn_pool(w_sum);
}

double MH_genetic::compute_weights(double r2_min_, double rho2in_, event_record ** recs_, int n_)
{
  double  wsum=0.0,
          r_min_=sqrt(r2_min_),
          rho_=sqrt(rho2in_);

  #pragma omp parallel for reduction(+:wsum)
  for (int i = 0; i < n_; i++)
  {
    wsum+=recs_[i]->w=gaussian_likelihood::compute_weight(sqrt(recs_[i]->get_r2()),r_min_,rho_);
  }
  return wsum;
}

double MH_genetic::compute_weights(double r2_min_, double rho2in_)
{
  double  wsum=0.0,
          r_min_=sqrt(r2_min_),
          rho_=sqrt(rho2in_);

  for (int i = 0; i < nlead; i++)
    wsum+=leaders[i]->w=gaussian_likelihood::compute_weight(sqrt(leaders[i]->get_r2()),r_min_,rho_);

  return wsum;
}

void MH_genetic::compute_weighted_ustats(double wsum_, event_record ** recs_, int n_)
{
  bool first2finish=true;
  #pragma omp parallel
  {
    int tid = thread_num();
    double *uwkspc_t = examiners[tid]->dub_wkspc;
    for (int i = 0; i < ulen; i++) uwkspc_t[i]=0.0;

    // compute the weighted mean
    #pragma omp for nowait
    for (int i = 0; i < n_; i++)
    {
      double wrati=(recs_[i]->w)/wsum_;
      for (int j = 0; j < ulen; j++) uwkspc_t[j]+=wrati_*(recs_[i]->u[j]);
    }
    #pragma omp critical
    {
      if (first2finish)
      {
        for (int j = 0; j < ulen; j++) u_wmean[j]=uwkspc_t[j];
        first2finish=false;
      }
      else for (int j = 0; j < ulen; j++) u_wmean[j]+=uwkspc_t[j];
    }

    for (int i = 0; i < ulen; i++) uwkspc_t[i]=0.0;

    #pragma omp barrier

    #pragma omp single
    {
      first2finish=true;
    }

    // compute the weighted variance (covariance diagonal)
    #pragma omp for nowait
    for (int i = 0; i < n_; i++)
    {
      double wrati=(recs_[i]->w)/wsum_;
      for (int j = 0; j < ulen; j++)
      {
        double udiff = recs_[i]->u[j]-u_wmean[j];
        uwkspc_t[j]+=wrati*udiff*udiff;
      }
    }
    #pragma omp critical
    {
      if (first2finish)
      {
        for (int j = 0; j < ulen; j++) u_wvar[j]=uwkspc_t[j];
        first2finish=false;
      }
      else for (int j = 0; j < ulen; j++) u_wvar[j]+=uwkspc_t[j];
    }
  }
}

void MH_genetic::set_leader_records()
{
  nreplace=take_records(leader_board,leaders,irepl_leaders,nlead);
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

void MH_genetic::clear_genetic_training_data()
{
  clear_basic_trainer_training_data();
}

void MH_genetic::consolidate_genetic_training_data()
{  

  double wsum = compute_weights();
  ncandidates=nsuccess;
  for (int i = 0; i < nsuccess; i++) candidates[i]=pool[isuccess_pool[i]];
  // if we have more successful particles than we can store, we have to narrow down the candidates
  if (nsuccess>nlead) pick_nbest_records(candidates,ncandidates=nlead,nsuccess);
  pick_nbest_records(leader_board,nlead,leader_count+ncandidates);
  set_leader_records();
  prob_best=gaussian_likelihood::compute_prob(sqrt(br2),sqrt(rho2));
  prob_worst=gaussian_likelihood::compute_prob(sqrt(wr2),sqrt(rho2));
}

void MH_genetic::write_Class_diagnostics(int &Class_count_)
{
  int hlen_Class=2;
  int header_Class[] = {hlen_Class, genetic_it_ilen_full(), genetic_it_dlen_full()};
  sprintf(obuf+obuf_end, "Class%d.mhdat",Class_count_);
  FILE * Class_file = fopen(obuf, "wb");
  fwrite(header_Class,sizeof(int),hlen_Class+1,Class_file);
  write_genetic_it_ints(Class_file);
  write_genetic_it_dubs(Class_file);
  leaders[0]->write_event_rec_training_header(Class_file,nlead);
  for (int i = 0; i < nlead; i++) leaders[i]->write_event_rec_training_data(Class_file);
  fclose(Class_file);
  Class_count_++;
}

void MH_genetic::write_generation_diagnostics(int &gen_count_)
{
  int hlen_gen=2;
  int header_gen[] = {hlen_gen, genetic_it_ilen_full(), genetic_it_dlen_full()};
  sprintf(obuf+obuf_end, "gen%d.mhdat",gen_count_);
  FILE * gen_file = fopen(obuf, "wb");
  write_genetic_it_ints(gen_file);
  write_genetic_it_dubs(gen_file);
  basic_MH_trainer::write_ustats(gen_file);
  fclose(gen_file);
  gen_count++;
}

void MH_genetic::write_ustats(FILE * file_)
{

}

void MH_genetic::respawn_pool(double w_sum_, int offset_)
{
  memset(ndup_leaders,0,nlead*sizeof(int));
  for (int i = 0; i < ulen; i++) u_var[i]=u_mean[i]=0.0;
  ndup=nredraw=ndup_unique=0;

  #pragma omp parallel
  {
    int tid = thread_num(),
        *dup_t = examiners[tid]->int_wkspc;
    double  *u_stat_t = examiners[tid]->dub_wkspc,
            inv_npool = 1.0/((double)npool),
            inv_npoolm1 = 1.0/((double)(npool-1));
    MH_rng * rng_t = rng[tid];

    memset(dup_t,0,nlead*sizeof(int));
    for (int i = 0; i < ulen; i++) u_stat_t[i]=0.0;

    #pragma omp for reduction(+:ndup) reduction(+:nredraw) nowait
    for (int i = 0; i < npool; i++)
    {
      if (i<offset_) pool[i]->take_record(leaders[i]); // reloading the leaders
      else
      {
        int j=0;
        double uni = (w_sum_/rs_full_factor)*rng_t->rand_uni();
        while ((j<leader_count)&&(uni>0.0)) uni-=w_leaders[j++];
        if (j>0)
        {
          if (uni<0.0)
          {ndup++; dup_t[--j]++; duplicate_u(pool[i],leaders[j],rng_t);}
          else
          {nredraw++; redraw_u(pool[i],rng_t);}
        }
        else
        {nredraw++; redraw_u(pool[i],rng_t);}
      }
      for (int i_u = 0; i_u < ulen; i_u++) u_stat_t[i_u]+=inv_npool*pool[i]->u[i_u];
    }
    #pragma omp critical
    {
      for (int i = 0; i < leader_count; i++) ndup_leaders[i]+=dup_t[i];
      for (int i = 0; i < ulen; i++) u_mean[i]+=u_stat_t[i];
    }
    for (int i = 0; i < ulen; i++) u_stat_t[i]=0.0;

    #pragma omp barrier

    #pragma omp for nowait
    for (int i = 0; i < npool; i++)
    {
      double *ui = pool[i]->u;
      for (int j = 0; j < ulen; j++)
      {
        double diff_j=ui[j]-u_mean[j];
        u_stat_t[i]+=(inv_npoolm1)*diff_j*diff_j;
      }
    }
    #pragma omp critical
    {
      for (int i = 0; i < ulen; i++) u_var[i]+=u_stat_t[i];
    }
  }

  for (int i = 0; i < leader_count; i++) if (ndup_leaders[i])
  {leaders_[i]->dup_count+=ndup_leaders[i]; ndup_unique++;}
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
