#ifndef MH_WORKERS_HH
#define MH_WORKERS_HH

#include "MH_learning.hh"

class MH_medic: public basic_thread_worker
{
  public:

    MH_medic(thread_work_struct &tws_, int thread_id_, double alpha_tol_, char * test_directory_);    
    ~MH_medic();

  private:

    const size_t buf_end;

    char * const mtest_buffer;

    double  * const TEST_pos,
            * const TEST_res,
            * const TEST_alpha,
            * const TEST_int;
};


#endif
