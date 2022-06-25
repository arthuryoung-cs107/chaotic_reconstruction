#include "MH_records.hh"

int comp_event_rec_ichunk_len(int nbeads_) {return nbeads_;}
int comp_event_rec_dchunk_len(int nbeads_) {return 4*nbeads_;}
int comp_event_rec_it_ichunk_len(int nbeads_) {return 0;}
int comp_event_rec_it_dchunk_len(int nbeads_) {return 2*nbeads_;}

event_record::event_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): basic_record(rs_, rid_, ichunk_, dchunk_, u_),
event_rec_ints(&nfobs), event_rec_dubs(&r2stable),
event_rec_it_ichunk_len(comp_event_rec_it_ichunk_len(nbeads)), event_rec_it_dchunk_len(comp_event_rec_it_dchunk_len(nbeads)),
evframe_bead(ichunk),
r2stable_bead(dchunk), netr2_regime(dchunk+nbeads), r2unstable_bead(dchunk+2*nbeads), alpha_bead(dchunk+3*nbeads)
{}

event_record::event_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double *umin_, double *umax_): basic_record(rs_, rid_, ichunk_, dchunk_, u_, ran_, umin_, umax_),
event_rec_ints(&nfobs), event_rec_dubs(&r2stable),
event_rec_it_ichunk_len(comp_event_rec_it_ichunk_len(nbeads)), event_rec_it_dchunk_len(comp_event_rec_it_dchunk_len(nbeads)),
evframe_bead(ichunk),
r2stable_bead(dchunk), netr2_regime(dchunk+nbeads), r2unstable_bead(dchunk+2*nbeads), alpha_bead(dchunk+3*nbeads)
{}

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

void event_record::determine_event_block(int &stev_earliest_, int &stev_latest_, int *stev_comp_, int *stev_ordered_, int *comps_ordered_)
{
  /* we're only actually setting stev_earliest_, stev_latest_, and stev_comp_ here. Everything else is just staging the event block ordering
  */
  int evframe_latest=0,
      evframe_earliest=Frames;
  for (int i = 0; i < nbeads; i++)
  {
    int evframe_it = evframe_bead[i];
    stev_comp_[i]=stev_ordered_[i]=evframe_it;
    comps_ordered_[i]=i;
    if (evframe_it<evframe_earliest) evframe_earliest=evframe_it;
    if (evframe_it>evframe_latest) evframe_latest=evframe_it;
  }
  stev_earliest_=evframe_earliest;
  stev_latest_=evframe_latest;
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
