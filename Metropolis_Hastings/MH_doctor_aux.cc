#include "MH_solvers.hh"

// MH_doctor

MH_doctor::MH_doctor(MH_train_struct &mhts_, int test_id_, int test_relay_id_, int Frames_test_, double alpha_tol_): MH_genetic(mhts_,0,0,0.0,alpha_tol_,1.0),
test_id(test_id_), test_relay_id(test_relay_id_), Frames_test(Frames_test_),
test_buffer(new char[io->obuf_len+100]),
TEST_refp(new double[2*Frames_test*nbeads]),
medics(new MH_medic*[nt])
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
  sprintf(test_buffer, "%s%s.MH%d_test%d_results/", io->obuf, io->data_name, test_relay_id, test_id);
  test_buf_end = strlen(test_buffer);

  #pragma omp parallel
  {
    medics[tid] = new MH_medic(*(examiners[thread_num()]), Frames_test, test_buffer, test_buf_end);
  }
}

MH_doctor::~MH_doctor()
{
  for (int i = 0; i < nt; i++) delete medics[i];
  delete [] medics;

  delete [] test_buffer;
  delete [] TEST_refp;
}

void MH_doctor::run(bool verbose_)
{
  bool first2finish=true;
  initialize_doctor_run();
  stage_doctor_diagnostics();

  clear_doctor_event_data();

  #pragma omp parallel
  {
    MH_medic *med_t = medics[thread_num()];
    med_t->clear_medic_event_data();
    #pragma omp for nowait
    for (int i = 0; i < npool; i++)
    {
      med_t->test_u(pool[i],i,verbose_);
    }
    med_t->consolidate_medic_event_data();
    #pragma omp critical
    {
      first2finish=med_t->report_medic_event_data(first2finish,stev_earliest,stev_latest,stev_comp,nev_state_comp, nobs_state_comp, mur2_state_comp, mualpha_state_comp);
    }
  }
  consolidate_doctor_event_data();
  event_block::define_event_block(sigma_scaled);
  synchronise_doctor_event_data();

  close_doctor_diagnostics();
}

void MH_doctor::initialize_doctor_run()
{
    initialize_genetic_run();

    // consider the reference data
    for (int frame_i = 0; frame_i < Frames_test; frame_i++)
    {
      int foffset = ndof*frame_i;
      double *f = xs+(foffset);
      for (int bead_i = 0, j=0; bead_i < nbeads; bead_i++, j+=dof)
      {TEST_refp[j+foffset]=f[j]; TEST_refp[j+1+foffset]=f[j+1];}
    }
}

void MH_doctor::stage_doctor_diagnostics()
{
  mkdir(test_buffer, S_IRWXU);
  printf("Made test directory: %s\n", test_buffer);

  int header_genetic_len=2;
  int header_genetic[] = {header_genetic_len, genetic_train_const_ilen, genetic_train_const_dlen};
  int header_doctor_len = 1;
  int header_doctor[] = {MHdoc_header_len, Frames_test};
  sprintf(test_buffer+test_buf_end, "startspecs.mhdat");
  FILE * startspecs_file = fopen(test_buffer, "wb");
  write_MH_params(startspecs_file);
  // write MH_genetic starting data
  fwrite(header_genetic,sizeof(int),header_genetic_len+1,startspecs_file);
  fwrite(genetic_train_const_ints, sizeof(int), header_genetic[1], startspecs_file);
  fwrite(genetic_train_const_dubs, sizeof(double), header_genetic[2], startspecs_file);
  // write MH_doctor starting data
  fwrite(header_doctor, sizeof(int), header_doctor_len+1, startspecs_file);
  fwrite(TEST_refp, sizeof(double), ndof*Frames_test, startspecs_file);
  fclose(startspecs_file);
}

void MH_doctor::synchronise_doctor_event_data()
{
  int nf_obs, nf_stable, nf_regime, nf_unstable;
  event_block::set_state_counts(nf_obs, nf_stable, nf_regime, nf_unstable);
  #pragma omp parallel
  {
    medics[thread_num()]->synchronise_medic_event_data(&nf_obs,stev_earliest,stev_latest,rho2stable,stev_comp,stev_ordered,comps_ordered,rho2stable_comp,delrho2_regime);
  }
}

void MH_doctor::close_doctor_diagnostics()
{
  int hlen=2;
  int header[] = {hlen, npool, stev_latest};
  sprintf(obuf+obuf_end, "event_block%d.mhdat",event_block_count_);
  FILE * data_file = fopen(obuf, "wb");
  fwrite(stev_comp, sizeof(int), nbeads, data_file);
  fwrite(stev_ordered, sizeof(int), nbeads, data_file);
  fwrite(comps_ordered, sizeof(int), nbeads, data_file);
  fwrite(nev_state_comp, sizeof(int), nbeads*Frames_test, data_file);
  fwrite(nobs_state_comp, sizeof(int), nbeads*Frames_test, data_file);
  fwrite(rho2stable_comp, sizeof(double), nbeads, data_file);
  fwrite(delrho2_regime, sizeof(double), nbeads, data_file);
  pool[0]->write_event_rec_full_header(data_file,npool);
  for (int i = 0; i < npool; i++) pool[i]->write_record_data(data_file);
  fclose(data_file);
}
