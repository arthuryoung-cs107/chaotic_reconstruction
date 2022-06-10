#ifndef MH_WORKERS_HH
#define MH_WORKERS_HH

#include "MH_records.hh"

class MH_medic: public basic_thread_worker
{
  public:

    MH_medic(thread_work_struct &tws_, int thread_id_, double alpha_tol_, int Frames_test_, char * test_directory_);
    ~MH_medic();

  private:

    const size_t buf_end;

    const int Frames_test;

    char * const mtest_buffer;

    double  * const TEST_p,
            * const TEST_r2,
            * const TEST_alpha,
            * const TEST_INTr2;
};


#endif
