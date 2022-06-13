#include "MH_records.hh"

int comp_event_rec_ichunk_len(int nbeads_) {return nbeads_;}
int comp_event_rec_dchunk_len(int nbeads_) {return 3*nbeads_;}

void event_record::write_ints(FILE * file_)
{
  basic_record::write_ints(file_);
  fwrite(event_rec_ints, sizeof(int), event_rec_ilen, file_);
}

void event_record::write_dubs(FILE * file_)
{
  basic_record::write_dubs(file_);
  fwrite(event_rec_dubs, sizeof(double), event_rec_dlen, file_);
}
