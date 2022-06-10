#ifndef MH_RECORDS_HH
#define MH_RECORDS_HH

#include "MH_learning.hh"

int comp_event_rec_ichunk_len(int nbeads_);
int comp_event_rec_dchunk_len(int nbeads_);
struct event_record : public basic_record
{
  event_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): basic_record(rs_, rid_, ichunk_, dchunk_, u_), evframe_bead(ichunk), r2stable_bead(dchunk), r2unstable_bead(dchunk+nbeads), alpha_bead(dchunk+2*nbeads) {}
  ~event_record() {}

  int * evframe_bead;

  double  r2stable,
          r2unstable,
          * r2stable_bead,
          * r2unstable_bead,
          * alpha_bead;
};

#endif
