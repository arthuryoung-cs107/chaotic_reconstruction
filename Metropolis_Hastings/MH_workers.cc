#include "MH_workers.hh"

MH_medic::MH_medic(thread_worker_struct &tws_, int thread_id_, double alpha_tol_, char * test_buffer_): basic_thread_worker(tws_, thread_id_, alpha_tol_),
buf_end(strlen(test_buffer_)), mtest_buffer(new char[buf_end+50]) 
{strcpy(mtest_buffer, test_buffer_);}

MH_medic
