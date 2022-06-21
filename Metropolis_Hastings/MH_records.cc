#include "MH_records.hh"

int comp_event_rec_ichunk_len(int nbeads_) {return nbeads_;}
int comp_event_rec_dchunk_len(int nbeads_) {return 3*nbeads_;}
int comp_event_rec_it_ichunk_len(int nbeads_) {return 0;}
int comp_event_rec_it_dchunk_len(int nbeads_) {return 2*nbeads_;}

int event_record::take_record(event_record *r_)
{
  if (basic_record::take_record(r_))
  {
    memcpy(event_rec_ints, r_->event_rec_ints, event_rec_ilen*sizeof(int));
    memcpy(event_rec_dubs, r_->event_rec_dubs, event_rec_dlen*sizeof(double));
    return 1;
  }
  else return 0;
}

void event_record::write_event_rec_full_header(FILE * file_, int len_)
{
  int hlen = 5;
  int header[] = {hlen, len_, event_rec_ilen_full(), event_rec_dlen_full(), ichunk_len, dchunk_len};
  fwrite(header, sizeof(int), hlen+1, file_);
}

void event_record::write_event_rec_training_header(FILE * file_, int len_)
{
  int hlen = 5;
  int header[] = {hlen, len_, event_rec_ilen_train(), event_rec_dlen_train(), event_rec_it_ichunk_len, event_rec_it_dchunk_len};
  fwrite(header, sizeof(int), hlen+1, file_);
}

void event_record::write_event_rec_training_data(FILE *file_)
{
  write_basic_rec_ints(file_); fwrite(event_rec_ints,sizeof(int),event_rec_it_ilen,file_);
  write_basic_rec_dubs(file_); fwrite(event_rec_dubs,sizeof(double),event_rec_it_dlen,file_);
  fwrite(ichunk,sizeof(int),event_rec_it_ichunk_len,file_);
  fwrite(dchunk,sizeof(double),event_rec_it_dchunk_len,file_);
}
