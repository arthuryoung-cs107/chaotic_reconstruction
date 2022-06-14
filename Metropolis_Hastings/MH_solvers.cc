#include "MH_solvers.hh"

#include <assert.h>
#include <sys/stat.h>

// MH_genetic

MH_genetic::MH_genetic(MH_train_struct &mhts_, double t_wheels0_, int gen_max_, double alpha_tol_): basic_MH_trainer(mhts_,comp_event_rec_ichunk_len(mhts_.get_par_nbeads()),comp_event_rec_dchunk_len(mhts_.get_par_nbeads()), t_wheels0_),
leader_gen_max(gen_max_), alpha_tol(alpha_tol_),
genetic_train_const_ints(&leader_gen_max), genetic_train_const_dubs(&alpha_tol),
genetic_train_ints(&leader_gen_count), genetic_train_dubs(&prob_best),
obuf(new char[io->obuf_len+100]),
r2_pool_Framebead(Tmatrix<double>(npool,Frames*nbeads)), alpha_pool_Framebead(Tmatrix<double>(npool,Frames*nbeads)),
mur2_Frame_bead(Tmatrix<double>(Frames, nbeads)), stdr2_Frame_bead(Tmatrix<double>(Frames, nbeads)),
mualpha_Frame_bead(Tmatrix<double>(Frames, nbeads)), stdalpha_Frame_bead(Tmatrix<double>(Frames, nbeads)),
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
  thread_worker_struct tws(&ulen,&dt_sim,ts,xs,d_ang,comega_s);
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
  free_Tmatrix<double>(mur2_Frame_bead); free_Tmatrix<double>(stdr2_Frame_bead);
  free_Tmatrix<double>(mualpha_Frame_bead); free_Tmatrix<double>(stdalpha_Frame_bead);

  delete [] obuf;
}

void MH_genetic::inspect_event_data()
{}

void MH_genetic::clear_event_data()
{
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

void MH_genetic::find_events()
{
  bool first2finish=true;
  int ncheck;
  event_record ** rec_check;
  if (gen_count==0)
  {
    rec_check=pool;
    ncheck=npool;
  }
  else
  {
    rec_check=leaders;
    ncheck=leader_count;
  }
  clear_event_data();
  #pragma omp parallel
  {
    MH_examiner *ex_t = examiner[thread_num()];
    ex_t->clear_event_data();
    #pragma omp for nowait
    for (int i = 0; i < ncheck; i++)
    {
      ex_t->detect_events(rec_check[i], r2_pool_Framebead[i], alpha_pool_Framebead[i]);
    }
    ex_t->consolidate_results();
    #pragma omp critical
    {
      first2finish=ex_t->report_results(mur2);
    }
  }
  inspect_event_data();
}

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

void MH_genetic::close_diagnostics()
{
  int header_len = 2;
  int header[] = {header_len, it_ilen_full(), it_dlen_full()};
  sprintf(obuf+obuf_end, "endspecs.mhdat");
  FILE * endspecs_file = fopen(obuf, "wb");
  fwrite(header_ints, sizeof(int), header_len+1, endspecs_file);
  write_it_ints(endspecs_file);
  write_it_dubs(endspecs_file);
  fwrite(genetic_train_it_ints, sizeof(int), genetic_train_it_ilen, endspecs_file);
  fwrite(genetic_train_it_dubs, sizeof(double), genetic_train_it_dlen, endspecs_file);
  records[0]->write_event_record_header(endspecs_file, nlead+npool);
  for (int i = 0; i < nlead+npool; i++) records[i]->write_record_data(endspecs_file);
  fclose(endspecs_file);
}

void MH_genetic::run(bool verbose_)
{
  initialize_run();
  stage_diagnostics();
  find_events();
  do
  {
    train_stable();
    train_unstable();
    if (check_convergence()) break;
    else find_events();
  } while(true);
  close_diagnostics();
}


// MH_doctor

MH_doctor::MH_doctor(MH_train_struct &mhts_, int test_id_, int test_relay_id_, int Frames_test_, double alpha_tol_): basic_MH_trainer(mhts_,comp_event_rec_ichunk_len(mhts_.get_par_nbeads()),comp_event_rec_dchunk_len(mhts_.get_par_nbeads())),
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
  thread_worker_struct tws(&ulen,&dt_sim,ts,xs,d_ang,comega_s);
  record_struct rs(ulen,nbeads,Frames,ichunk_width,dchunk_width);

  #pragma omp parallel
  {
    int tid=thread_num();
    MH_rng * ran_t = rng[tid];
    MH_medic * med_t = medics[tid] = new MH_medic(sp_min,pg[tid],wl,tws,tid,alpha_tol,Frames_test,test_buffer);
    #pragma omp for
      for (int i = 0; i < nlead+npool; i++)
        records[i] = new event_record(rs, i, ichunk[i], dchunk[i], uchunk[i]);
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
    med_t->initialize_utest();
    #pragma omp for nowait
      for (int i = 0; i < npool; i++)
      {
        med_t->test_u(pool[i],i, verbose_);
      }
    #pragma omp critical
    {
      first2finish = med_t->report_results(first2finish, evcount_bead_frame);
    }
  }
  close_diagnostics();
}
void MH_doctor::stage_diagnostics()
{
  for (int frame_i = 0; frame_i < Frames_test; frame_i++)
  {
    int foffset = 2*nbeads*frame_i;
    double *f = xs+(foffset);
    for (int bead_i = 0, j=0; bead_i < nbeads; bead_i++, j+=2)
    {TEST_refp[j+foffset]=f[j]; TEST_refp[j+1+foffset]=f[j+1];}
  }
  int MH_doc_header_len = 1;
  int header_ints[] = {MH_doc_header_len, Frames_test};
  sprintf(test_buffer+test_buf_end, "results_startspecs.redat");
  FILE * test_startspecs = fopen(test_buffer, "wb");
  write_MH_params(test_startspecs);
  fwrite(header_ints, sizeof(int), MH_doc_header_len+1, test_startspecs);
  fwrite(TEST_refp, sizeof(double), 2*nbeads*Frames_test, test_startspecs);
  fclose(test_startspecs);
}
void MH_doctor::close_diagnostics()
{
  int MH_doc_header_len = 0;
  int header_ints[] = {MH_doc_header_len};
  sprintf(test_buffer+test_buf_end, "results_endspecs.redat");
  FILE * test_endspecs = fopen(test_buffer, "wb");
  fwrite(header_ints, sizeof(int), MH_doc_header_len+1, test_endspecs);
  fwrite(evcount_bead_frame[0], sizeof(int), 2*nbeads*Frames_test, test_endspecs);
  fclose(test_endspecs);
}
