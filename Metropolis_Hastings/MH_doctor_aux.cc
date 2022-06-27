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
  sprintf(test_buffer, "%s%s.re%d_test%d_results/", io->obuf, io->data_name, test_relay_id, test_id);
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
      first2finish=med_t->report_results(first2finish,nev_state_comp);
    }
  }
  close_diagnostics();
}

void MH_doctor::stage_diagnostics()
{
  mkdir(test_buffer, S_IRWXU);
  printf("Made test directory: %s\n", test_buffer);

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
