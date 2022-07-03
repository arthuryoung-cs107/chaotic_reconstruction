#include "MH_solvers.hh"
#include <sys/stat.h>

// MH_genetic

MH_genetic::MH_genetic(MH_train_struct &mhts_, int Class_max_, int gen_max_, int itrain_max_, double t_wheels0_, double alpha_tol_, double rs_full_factor_, double train_tol_): basic_MH_trainer(mhts_,comp_event_rec_ichunk_len(mhts_.get_par_nbeads()),comp_event_rec_dchunk_len(mhts_.get_par_nbeads()), t_wheels0_), event_block(nbeads, Frames),
obuf(new char[io->obuf_len+100]),
Class_max(Class_max_), gen_max(gen_max_), itrain_max(itrain_max_),
alpha_tol(alpha_tol_), rs_full_factor(rs_full_factor_), train_tol(train_tol_),
r2_pool_Framebead(Tmatrix<double>(npool,Frames*nbeads)), alpha_pool_Framebead(Tmatrix<double>(npool,Frames*nbeads)),
examiners(new MH_examiner*[nt]),
records(new event_record*[nlead+npool]), leaders(records), pool(records+nlead),
leader_board(new event_record*[nlead+npool]), candidates(leader_board+nlead)
{
  // constant structures for initializing thread workers and records
  thread_worker_struct tws(ulen,nbeads,Frames,nlead,npool,dt_sim,t_phys,ts,xs,d_ang,comega_s);
  record_struct rs(ulen,nbeads,Frames,ichunk_width,dchunk_width);

  #pragma omp parallel
  {
    int tid=thread_num();
    MH_rng * ran_t = rng[tid];
    MH_examiner * ex_t = examiners[tid] = new MH_examiner(sp_min,pg[tid],wl,tws,tid,alpha_tol);
    #pragma omp for nowait
    for (int i = 0; i < nlead; i++)
    {
      records[i] = new event_record(rs, i, ichunk[i], dchunk[i], uchunk[i]);
    }
    for (int i = nlead; i < nlead+npool; i++)
    {
      records[i] = new event_record(rs, i, ichunk[i], dchunk[i], uchunk[i], ran_t, umin, umax);
    }
  }
}

MH_genetic::~MH_genetic()
{
  delete [] leader_board;
  for (int i = 0; i < npool+nlead; i++) delete records[i];
  delete [] records;

  for (int i = 0; i < nt; i++) delete examiners[i];
  delete [] examiners;

  free_Tmatrix<double>(r2_pool_Framebead); free_Tmatrix<double>(alpha_pool_Framebead);

  delete [] obuf;
}

void MH_genetic::run(bool verbose_)
{
  initialize_genetic_run(); // set various counts to zero
  stage_diagnostics(); // write startup data for post processing
  do
  {
    bool stable_conv, unstable_conv;
    find_events(verbose_); // find critical frames for residual blowup
    train_event_block(verbose_,stable_conv,unstable_conv); // train around the critical frames
    if (check_run_convergence(stable_conv,unstable_conv)) break; // check convergence
    else stage_event_search(); // prepare for a new round of critical frame detection
  } while(true);
  close_diagnostics(); // write data from end of run
}

void MH_genetic::stage_diagnostics()
{
  // prepare output information
  sprintf(obuf, "%s%s.MH%d_results/", io->obuf,io->data_name,io->id);
  obuf_end = strlen(obuf);
  mkdir(obuf, S_IRWXU);
  printf("Made test directory: %s\n", obuf);

  int header_len = 2;
  int header[] = {header_len, genetic_train_const_ilen, genetic_train_const_dlen};
  sprintf(obuf+obuf_end, "startspecs.mhdat");
  FILE * startspecs_file = fopen(obuf, "wb");
  write_MH_params(startspecs_file); // write ulen, # of beads, # of frames, # of leader particles and pool particles
  fwrite(header, sizeof(int), header_len+1, startspecs_file);
  fwrite(genetic_train_const_ints, sizeof(int), genetic_train_const_ilen, startspecs_file); // write max gen and Class counts
  fwrite(genetic_train_const_dubs, sizeof(double), genetic_train_const_dlen, startspecs_file); // write alpha tolerance and the resampling full factor
  fclose(startspecs_file);
}

void MH_genetic::find_events(bool verbose_)
{
  bool first2finish=true;

  clear_genetic_event_data(); // clear event block data and clear alpha and residual statistics buffers
  #pragma omp parallel
  {
    MH_examiner *ex_t = examiners[thread_num()];
    ex_t->clear_examiner_event_data();
    #pragma omp for nowait
    for (int i = 0; i < npool; i++)
    {
      ex_t->detect_events(pool[i], r2_pool_Framebead[i], alpha_pool_Framebead[i]);
    }
    ex_t->consolidate_examiner_event_data();
    #pragma omp critical
    {
      first2finish=ex_t->report_examiner_event_data(first2finish,
        stev_earliest,stev_latest,
        stev_comp,nev_state_comp,nobs_state_comp,
        mur2_state_comp,mualpha_state_comp);
    }
  }

  // sort event states chronologically, given stev_comp is already set
  if (gen_count>0) consolidate_genetic_event_data(bleader_rid); // using our best guess
  else consolidate_genetic_event_data(); // using conservative estimates

  // compute expected residuals using presumed noise level
  define_genetic_event_block(sigma_scaled);
  synchronise_genetic_event_data(); // set event data of thread workers to the consolidated values
  report_genetic_event_data(); // finish event stats and write out results
}

void MH_genetic::train_event_block(bool verbose_, bool &stable_convergence_, bool &unstable_convergence_)
{
  int nit_train=0;

  // perform stable training
  int nit_stable_train=0;
  double rho2_stable_local=rho2=set_stable_objective(r2_scale);
  stable_convergence_=train_objective(verbose_,nit_train,nit_stable_train,rho2_stable_local);
  // perform unstable training
  int nit_unstable_train=0;
  double rho2_unstable_local=rho2=set_unstable_objective(r2_scale);
  unstable_convergence_=train_objective(verbose_,nit_train,nit_unstable_train,rho2_unstable_local);
}

bool MH_genetic::train_objective(bool verbose_, int &nit_, int &nit_objective_, double rho2_)
{
  bool training_success;
  do
  {
    int nsuccess_local=nsuccess=0;
    bool first2finish=true;
    double wsum_pool=0.0;
    clear_genetic_training_data();
    #pragma omp parallel
    {
      MH_examiner *ex_t=examiners[thread_num()];
      ex_t->clear_examiner_training_data();
      double  r_scale=sqrt(r2_scale),
              rho=sqrt(rho2_);

      #pragma omp for reduction(+:nsuccess_local) reduction(+:wsum_pool) nowait
      for (int i = 0; i < npool; i++)
      {
        if (ex_t->examine_u(pool[i],i,wr2)) nsuccess_local++;
        wsum_pool+=pool[i]->w=gaussian_likelihood::compute_weight(sqrt(pool[i]->get_r2()),r_scale,rho);
      }
      ex_t->consolidate_examiner_training_data(pool);
      #pragma omp critical
      {
        first2finish=ex_t->report_examiner_training_data(first2finish,&bpool,isuccess_pool,nsuccess,u_wmean);
      }
    }
    double wsum_leaders = consolidate_genetic_training_data(wsum_pool,rho2_,nreplace,r2_scale);
    report_genetic_training_data(nreplace,Class_count,gen_count);
    if (check_objective_convergence(++nit_, ++nit_objective_, training_success)) break;
    else respawn_pool(wsum_leaders);
  } while (true);
  return training_success;
}

bool MH_genetic::check_objective_convergence(int nit_, int nit_objective_, bool &training_success_)
{
  if (br2<rho2*(1.0+train_tol))
  {
    training_success_=true;
    return true;
  }
  else if ((nit_objective_>=itrain_max)||(gen_count>=gen_max))
  {
    training_success_=false;
    return true;
  }
  else return false;
}

bool MH_genetic::check_run_convergence(bool stable_conv_, bool unstable_conv_)
{
  if (event_block::check_stev_convergence()&&stable_conv_&&unstable_conv_)
  {
    printf("We win.\n");
    return true;
  }
  else if (gen_count>=gen_max)
  {
    printf("max generations reached.\n");
    return true;
  }
  else return false;
}
