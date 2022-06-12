#include "MH_workers.hh"

MH_medic::MH_medic(swirl_param  &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_, int Frames_test_, char * test_buffer_): basic_thread_worker(sp_, pg_, wl_, tws_, thread_id_, alpha_tol_),
Frames_test(Frames_test_),
buf_end(strlen(test_buffer_)), mtest_buffer(new char[buf_end+50]),
TEST_p(new double[2*nbeads*Frames_test]), TEST_r2(new double[nbeads*Frames_test]), TEST_alpha(new double[nbeads*Frames_test]), TEST_INTr2(new double[nbeads*Frames_test])
{strcpy(mtest_buffer, test_buffer_);}

MH_medic::~MH_medic()
{
  delete [] mtest_buffer;
  delete [] TEST_p;
  delete [] TEST_r2;
  delete [] TEST_alpha;
  delete [] TEST_INTr2;
}
