#include "MH_solvers.hh"
#include <assert.h>

MH_doctor::MH_doctor(MH_train_struct &mhts_, int test_id_, int test_relay_id_, int Frames_test_, double alpha_tol_): basic_MH_trainer(mhts_,0.0),
test_id(test_id_), test_relay_id(test_relay_id_), Frames_test(Frames_test_),
alpha_tol(alpha_tol_),
test_buffer(new char[mhts_.get_io_obuf_len()+100]),
TEST_ref_pos(new double[2*Frames_test*mhts_.get_nbeads()]),
medics(new MH_medic*[mhts_.get_nt()]),
records(new event_record*[mhts_.get_record_len()])
{
  // write the name of the input data file
  sprintf(test_buffer, "%s%s.re%d_test%d.redat", io->obuf,io->data_name,test_relay_id);
  FILE * test_file = fopen(test_buffer, "r");
  printf("reading matlab test file: %s\n", test_buffer);

  // read the input parameters
  int header[2];
  fread_safe(header, sizeof(int), 2, test_file);
  assert((header[0]==param_len)&&(header[1]<=npool));
  fread_safe(param_chunk[nlead], sizeof(double), header[0]*header[1], test_file);
  fclose(test_file);

  // make the output directory
  sprintf(test_buffer, "%s%s.re%d_test%d_results/", rep->out_buf, rep->file_name, test_relay_id_, test_id_);
  mkdir(test_buffer, S_IRWXU);
  printf("Made test directory: %s\n", test_buffer);

  thread_work_struct tws(&param_len, &dt_sim,ts,xs,d_ang,comega_s);

  #pragma omp parallel
  {
    int tid=thread_num();
    MH_medic * med_t = medics[tid] = new MH_medic(tws,tid,alpha_tol,test_buffer);

  }
}
