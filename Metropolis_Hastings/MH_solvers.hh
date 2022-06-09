#ifndef MH_SOLVERS_HH
#define MH_SOLVERS_HH

#include "MH_workers.hh"

const int event_struct_ilen=0;
const int event_struct_dlen=0;

struct event_record : public basic_record
{
  event_struct();
};

class MH_doctor : public basic_MH_trainer
{
  public:

    MH_doctor(MH_train_struct &mhts_, int test_id_, int Frames_test_, int test_relay_id_, double alpha_tol_);    
    ~MH_doctor();

  private:

    char * const test_buffer;

    const int test_id,
              test_relay_id,
              Frames_test;

    const double alpha_tol;

    double  * const TEST_ref_pos;

    MH_medic ** const medics;
    event_record ** const records;
};

#endif
