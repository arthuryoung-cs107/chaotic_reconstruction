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

void MH_genetic::synchronise_genetic_event_data()
{
  #pragma omp parallel
  {
    examiners[thread_num()]->synchronise_examiner_event_data( stev_earliest,
                                                              stev_latest,
                                                              stev_comp,
                                                              stev_ordered,
                                                              comps_ordered,
                                                              rho2stable_comp,
                                                              rho2unstable_comp);
  }
}

void MH_genetic::report_genetic_event_data()
{
  // begin by finishing computation of event detection statistics
  #pragma omp parallel
  {
    // finish computing mean residuals for each component of the state space
    #pragma omp for
    for (int i = 0; i < nbeads*max(stev_latest,stev_min); i++)
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
  write_event_diagnostics(event_block_count++);
}

void MH_genetic::write_event_diagnostics(int event_block_count_)
{
  // write event data
  int hlen=3;
  int header[] = {hlen, npool, stev_latest, write_pool_event_data};
  sprintf(obuf+obuf_end, "event_block%d.mhdat",event_block_count_);
  FILE * data_file = fopen(obuf, "wb");
  fwrite(header, sizeof(int), hlen+1, data_file);
  fwrite(stev_comp, sizeof(int), nbeads, data_file);
  fwrite(stev_ordered, sizeof(int), nbeads, data_file);
  fwrite(comps_ordered, sizeof(int), nbeads, data_file);
  fwrite(nev_state_comp, sizeof(int), nbeads*stev_latest, data_file);
  fwrite(nobs_state_comp, sizeof(int), nbeads*stev_latest, data_file);
  fwrite(rho2stable_comp, sizeof(double), nbeads, data_file);
  fwrite(rho2unstable_comp, sizeof(double), nbeads, data_file);
  fwrite(mur2_state_comp, sizeof(double), nbeads*stev_latest, data_file);
  fwrite(stdr2_state_comp, sizeof(double), nbeads*stev_latest, data_file);
  fwrite(mualpha_state_comp, sizeof(double), nbeads*stev_latest, data_file);
  fwrite(stdalpha_state_comp, sizeof(double), nbeads*stev_latest, data_file);
  pool[0]->write_event_rec_full_header(data_file,npool);
  for (int i = 0; i < npool; i++) pool[i]->write_record_data(data_file);
  if (write_pool_event_data)
  {
    for (int i = 0; i < npool; i++)
      fwrite(r2_pool_Framebead[i], sizeof(double), nbeads*stev_latest, data_file);
    for (int i = 0; i < npool; i++)
      fwrite(alpha_pool_Framebead[i], sizeof(double), nbeads*stev_latest, data_file);
  }
  fclose(data_file);
}

double MH_genetic::set_objective(bool verbose_, double &r2_scale_, bool stable_flag_)
{
  double r2_min=DBL_MAX;
  #pragma omp parallel reduction(min:r2_min)
  {
    MH_examiner *ex_t = examiners[thread_num()];
    if (stable_flag_)
    {
      ex_t->set_stable_objective();
      #pragma omp for nowait
      for (int i = 0; i < npool; i++)
      {
        // start by using the r2_pool_Framebead data to adjust the records back to the correct event data
        ex_t->restore_event_record(pool[i],r2_pool_Framebead[i],alpha_pool_Framebead[i]);
        double r2_it = pool[i]->set_record_stable();
        if (r2_it<r2_min) r2_min=r2_it;
        leader_board[i]=pool[i];
      }
      #pragma omp for
      for (int i = 0; i < nlead; i++)
      {
        double r2_work = leaders[i]->set_record_stable();
      }
    }
    else
    {
      ex_t->set_unstable_objective();
      #pragma omp for
      for (int i = 0; i < npool+nlead; i++)
      {
        double r2_it = records[i]->set_record_unstable();
        if (r2_it<r2_min) r2_min=r2_it;
        leader_board[i]=records[i];
      }
    }
  }

  event_record ** records_use;
  int n_use;
  if (stable_flag_)
  {
    rho2_objective=rho2=event_block::compute_netrho2stable();
    records_use=pool;
    n_use=npool;
  }
  else
  {
    rho2_objective=rho2=event_block::compute_netrho2unstable();
    records_use=records;
    n_use=nlead+npool;
  }

  if (verbose_) verbose_set_objective_1();

  double  rho2obj=rho2_objective,
          wsum_full = compute_weights(r2_min,rho2obj,records_use,n_use);

  compute_weighted_ustats(wsum_full,records_use,n_use);
  pick_nbest_records(leader_board,nlead,n_use);
  r2_scale_=set_leader_records();

  if (verbose_) verbose_set_objective_2();

  report_genetic_training_data(nreplace,Class_count,gen_count);
  double wsum_leaders=0.0; // given that weights have already been computed, just sum up leader terms

  #pragma omp parallel for reduction(+:wsum_leaders)
  for (int i = 0; i < nlead; i++)
    wsum_leaders+=w_leaders[i]=leaders[i]->w;

  respawn_pool(verbose_, wsum_leaders, w_leaders);
  return rho2obj;
}

double MH_genetic::compute_weights(double r2_min_, double rho2in_, event_record ** recs_, int n_)
{
  double  wsum=0.0,
          r_min_=sqrt(r2_min_),
          rho_=sqrt(rho2in_);

  #pragma omp parallel for reduction(+:wsum)
  for (int i = 0; i < n_; i++)
    wsum+=recs_[i]->w=gaussian_likelihood::compute_weight(sqrt(recs_[i]->get_r2()),r_min_,rho_);

  return wsum;
}

double MH_genetic::compute_weights(double r2_min_, double rho2in_, double *w_leaders_)
{
  double  wsum=0.0,
          r_min_=sqrt(r2_min_),
          rho_=sqrt(rho2in_);

  for (int i = 0; i < nlead; i++)
    wsum+=w_leaders_[i]=leaders[i]->w=gaussian_likelihood::compute_weight(sqrt(leaders[i]->get_r2()),r_min_,rho_);

  return wsum;
}

void MH_genetic::compute_weighted_ustats(double wsum_, event_record ** recs_, int n_)
{
  bool  first2finish1=true,
        first2finish2=true;
  #pragma omp parallel
  {
    int tid = thread_num();
    double *uwkspc_t = examiners[tid]->ustat_buf;
    for (int i = 0; i < ulen; i++) uwkspc_t[i]=0.0;

    // compute the weighted mean
    #pragma omp for nowait
    for (int i = 0; i < n_; i++)
    {
      double wrati=(recs_[i]->w)/wsum_;
      for (int j = 0; j < ulen; j++) uwkspc_t[j]+=wrati*(recs_[i]->u[j]);
    }
    #pragma omp critical
    {
      if (first2finish1)
      {
        for (int j = 0; j < ulen; j++) u_wmean[j]=uwkspc_t[j];
        first2finish1=false;
      }
      else for (int j = 0; j < ulen; j++) u_wmean[j]+=uwkspc_t[j];
    }

    for (int i = 0; i < ulen; i++) uwkspc_t[i]=0.0;

    #pragma omp barrier

    // compute the weighted variance (covariance diagonal)
    #pragma omp for nowait
    for (int i = 0; i < n_; i++)
    {
      double wrati=(recs_[i]->w)/(wsum_);
      for (int j = 0; j < ulen; j++)
      {
        double udiff = recs_[i]->u[j]-u_wmean[j];
        uwkspc_t[j]+=wrati*udiff*udiff;
      }
    }
    #pragma omp critical
    {
      if (first2finish2)
      {
        for (int j = 0; j < ulen; j++) u_wvar[j]=uwkspc_t[j];
        first2finish2=false;
      }
      else for (int j = 0; j < ulen; j++) u_wvar[j]+=uwkspc_t[j];
    }
  }
}

double MH_genetic::set_leader_records()
{
  // assumes that leader board has been already sorted to the best performing particles
  nreplace=take_records(leader_board,leaders,irepl_leaders,nlead);

  for (int i = 0; i < nreplace; i++)
  {
    int repl_index = irepl_leaders[i];
    // make sure that first bit of leaderboard is always pointing to leaders
    leader_board[repl_index]=leaders[repl_index];
    leaders[repl_index]->Class=Class_count;
  }
  leader_count=nlead;
  int bleader_rid_local=0,
      wleader_rid_local=0;
  for (int i = 1; i < nlead; i++)
  {
    if (leaders[bleader_rid_local]->isworse(leaders[i])) bleader_rid_local=i; // leader i is better than the current best record
    if (leaders[wleader_rid_local]->isbetter(leaders[i])) wleader_rid_local=i; // leader i is worse than the current worst record
  }
  bleader=leaders[bleader_rid_local]; wleader=leaders[wleader_rid_local];

  bleader_rid=bleader->rid; wleader_rid=wleader->rid;
  br2=bleader->get_r2(); wr2=wleader->get_r2();
  return wr2;
}

double MH_genetic::consolidate_genetic_training_data(double wsum_pool_,double * w_leaders_, double rho2_,int &nreplace_,double &r2_scale_)
{
  // begin by completing weighted statistics
  double inv_w= 1.0/wsum_pool_;
  for (int i = 0; i < ulen; i++) u_wmean[i]*=inv_w;

  bool first2finish=true;

  #pragma omp parallel
  {
    double *ustat_t = examiners[thread_num()]->ustat_buf;
    for (int i = 0; i < ulen; i++) ustat_t[i]=0.0;

    #pragma omp for
    for (int i_pool = 0; i_pool < npool; i_pool++)
    {
      double  wrat_i=pool[i_pool]->w*inv_w,
              *ui=pool[i_pool]->u;
      for (int i_u = 0; i_u < ulen; i_u++)
      {
        double u_diff = ui[i_u]-u_wmean[i_u];
        ustat_t[i_u]+=wrat_i*u_diff*u_diff;
      }
    }
    #pragma omp critical
    {
      if (first2finish)
      {
        for (int i = 0; i < ulen; i++) u_wvar[i]=ustat_t[i];
        first2finish=false;
      }
      else for (int i = 0; i < ulen; i++) u_wvar[i]+=ustat_t[i];
    }
  }

  // now update the leaders
  ncandidates=nsuccess;
  for (int i = 0; i < nsuccess; i++) candidates[i]=pool[isuccess_pool[i]];

  // if we have more successful particles than we can store, we have to narrow down the candidates
  if (nsuccess>nlead) pick_nbest_records(candidates,ncandidates=nlead,nsuccess);
  pick_nbest_records(leader_board,nlead,nlead+ncandidates);
  r2_scale_=set_leader_records();
  prob_best=gaussian_likelihood::compute_prob(sqrt(br2),sqrt(rho2_));
  prob_worst=gaussian_likelihood::compute_prob(sqrt(wr2),sqrt(rho2_));
  return compute_weights(r2_scale_,rho2_,w_leaders_);
}

void MH_genetic::write_Class_diagnostics(int Class_count_)
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
}

void MH_genetic::write_generation_diagnostics(int gen_count_)
{
  int hlen_gen=2;
  int header_gen[] = {hlen_gen, genetic_it_ilen_full(), genetic_it_dlen_full()};
  sprintf(obuf+obuf_end, "gen%d.mhdat",gen_count_);
  FILE * gen_file = fopen(obuf, "wb");
  fwrite(header_gen,sizeof(int),hlen_gen+1,gen_file);
  write_genetic_it_ints(gen_file);
  write_genetic_it_dubs(gen_file);
  basic_MH_trainer::write_ustats(gen_file);
  if (write_training_pool)
  {
    pool[0]->write_event_rec_training_header(gen_file,npool);
    for (int i = 0; i < npool; i++) pool[i]->write_event_rec_training_data(gen_file);
  }
  fclose(gen_file);
}

void MH_genetic::close_diagnostics()
{
  int header_len = 2;
  int header[] = {header_len, genetic_it_ilen_full(), genetic_it_dlen_full()};
  sprintf(obuf+obuf_end, "endspecs.mhdat");
  FILE * endspecs_file = fopen(obuf, "wb");
  fwrite(header, sizeof(int), header_len+1, endspecs_file);
  write_genetic_it_ints(endspecs_file);
  write_genetic_it_dubs(endspecs_file);
  records[0]->write_event_rec_full_header(endspecs_file, nlead+npool);
  for (int i = 0; i < nlead+npool; i++) records[i]->write_record_data(endspecs_file);
  fclose(endspecs_file);
}
