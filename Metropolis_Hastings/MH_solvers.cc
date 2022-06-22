#include "MH_solvers.hh"

#include <assert.h>
// #include <sys/stat.h>

// MH_genetic

MH_genetic::MH_genetic(MH_train_struct &mhts_, int Class_max_, int gen_max_, double t_wheels0_, double alpha_tol_, double rs_full_factor_): basic_MH_trainer(mhts_,comp_event_rec_ichunk_len(mhts_.get_par_nbeads()),comp_event_rec_dchunk_len(mhts_.get_par_nbeads()), t_wheels0_), event_block(nbeads, Frames),
Class_max(Class_max_), gen_max(gen_max_),
alpha_tol(alpha_tol_), rs_full_factor(rs_full_factor_),
genetic_train_const_ints(&Class_max), genetic_train_const_dubs(&alpha_tol),
genetic_train_ints(&Class_count), genetic_train_dubs(&prob_best),
obuf(new char[io->obuf_len+100]),
r2_pool_Framebead(Tmatrix<double>(npool,Frames*nbeads)), alpha_pool_Framebead(Tmatrix<double>(npool,Frames*nbeads)),
examiners(new MH_examiner*[nt]), records(new event_record*[nlead+npool]),
leaders(records), pool(records+nlead),
leader_board(new event_record*[nlead+npool]), candidates(leader_board+nlead)
{
  // prepare output information
  sprintf(obuf, "%s%s.MH%d_results/", io->obuf,io->data_name,io->id);
  obuf_end = strlen(obuf);
  mkdir(obuf, S_IRWXU);
  printf("Made test directory: %s\n", obuf);

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
  initialize_genetic_run();
  stage_diagnostics();
  find_events();
  do
  {
    train_event_block(verbose_);
    if (check_convergence()) break;
    else find_events(verbose_);
  } while(true);
  close_diagnostics();
}

void MH_genetic::find_events(bool verbose_)
{
  bool first2finish=true;
  if (gen_count>0) reload_leaders(bleader_rid);

  clear_genetic_event_data();
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
      first2finish=ex_t->report_examiner_event_data(first2finish,stev_earliest,stev_latest,stev_comp_,nev_state_comp, nobs_state_comp, mur2_state_comp, mualpha_state_comp);
    }
  }

  // sort event states chronologically, given stev_comp is already set
  if (gen_count>0) consolidate_genetic_event_data(bleader_rid); // using our best guess
  else consolidate_genetic_event_data(); // using conservative estimates

  // compute expected residuals using presumed noise level
  event_block::define_event_block(sigma_scaled);
  synchronise_genetic_event_data(); // set event data of thread workers to the consolidated values
  report_genetic_event_data(); // finish event stats and write out results
  set_genetic_stable_objective(); // set the expected residual to be that expected from the stable data
  post_event_resampling(); // restore records, sort by performance, and redraw generation
}

void MH_genetic::train_event_block(bool verbose_)
{
  int nit_train=0;

  for (int i_regime = 0; i_regime < nbeads; i_regime++) // walking through event block
  {
    set_regime_objective(i_regime);
    int nit_regime=0;
    do
    {
      nit_train++;
      bool first2finish=true;
      int nsuccess_local=0;
      clear_genetic_training_data();
      #pragma omp parallel
      {
        MH_examiner * ex_t=examiners[thread_num()];
        ex_t->clear_examiner_training_data();
        #pragma omp for reduction(+:nsuccess_local) nowait
        for (int i = 0; i < npool; i++)
        {
          if (ex_t->examine_u(pool[i],i,wr2)) nsuccess_local++;
        }
        ex_t->consolidate_examiner_training_data();
        #pragma omp critical
        {
          first2finish=ex_t->report_examiner_training_data(first2finish,isuccess_pool,nsuccess);
        }
      }
      consolidate_genetic_training_data();
      report_genetic_training_data();
      if (check_regime_convergence(++nit_train_stable)) break;
      else respawn_pool();
    } while (true)
  }
}

void MH_genetic::check_event_block_convergence()
{

}

void MH_genetic::check_convergence()
{

}

// MH_doctor

MH_doctor::MH_doctor(MH_train_struct &mhts_, int test_id_, int test_relay_id_, int Frames_test_, double alpha_tol_): basic_MH_trainer(mhts_,comp_event_rec_ichunk_len(mhts_.get_par_nbeads()),comp_event_rec_dchunk_len(mhts_.get_par_nbeads())), event_block(nbeads, Frames),
test_id(test_id_), test_relay_id(test_relay_id_), Frames_test(Frames_test_),
alpha_tol(alpha_tol_),
test_buffer(new char[io->obuf_len+100]),
TEST_refp(new double[2*Frames_test*nbeads]),
medics(new MH_medic*[nt]), records(new event_record*[nlead+npool]),
leaders(records), pool(records+nlead)
{
  // write the name of the input data file
  sprintf(test_buffer, "%s%s.re%d_test%d.redat", io->obuf,io->data_name,test_relay_id,test_id);
  FILE * test_file = fopen(test_buffer, "r");
  printf("reading matlab test file: %s\n", test_buffer);

  // read the input parameters
  int header[2];
  fread_SAFE(header, sizeof(int), 2, test_file);
  assert((header[0]==ulen)&&(header[1]==npool));
  fread_SAFE(uchunk[nlead], sizeof(double), header[0]*header[1], test_file);
  fclose(test_file);

  // make the output directory
  sprintf(test_buffer, "%s%s.re%d_test%d_results/", io->obuf, io->data_name, test_relay_id_, test_id_);
  test_buf_end = strlen(test_buffer);
  mkdir(test_buffer, S_IRWXU);
  printf("Made test directory: %s\n", test_buffer);

  // constant structures for initializing thread workers and records
  thread_worker_struct tws(ulen,nbeads,Frames,nlead,npool,dt_sim,t_phys,ts,xs,d_ang,comega_s);
  record_struct rs(ulen,nbeads,Frames,ichunk_width,dchunk_width);

  #pragma omp parallel
  {
    int tid=thread_num();
    MH_rng * ran_t = rng[tid];
    MH_medic * med_t = medics[tid] = new MH_medic(sp_min,pg[tid],wl,tws,tid,alpha_tol,Frames_test,test_buffer);
    #pragma omp for
    {
      for (int i = 0; i < nlead+npool; i++)
        records[i] = new event_record(rs, i, ichunk[i], dchunk[i], uchunk[i]);
    }
  }
}

MH_doctor::~MH_doctor()
{
  for (int i = 0; i < npool+nlead; i++) delete records[i];
  delete [] records;

  for (int i = 0; i < nt; i++) delete medics[i];
  delete [] medics;

  delete [] test_buffer;
  delete [] TEST_refp;
}

void MH_doctor::run(bool verbose_)
{
  bool first2finish=true;
  initialize_run();
  stage_diagnostics();
  #pragma omp parallel
  {
    int tid=thread_num();
    MH_medic *med_t = medics[tid];
    med_t->clear_event_data();
    #pragma omp for nowait
    for (int i = 0; i < npool; i++)
    {
      med_t->test_u(pool[i],i, verbose_);
    }
    #pragma omp critical
    {
      first2finish=med_t->report_event_data(first2finish,stev_earliest,stev_latest,nev_state_comp, nobs_state_comp, mur2_state_comp, mualpha_state_comp);
    }
  }
  close_diagnostics();
}
