#include "MH_records.hh"

int comp_event_rec_ichunk_len(int nbeads_) {return nbeads_;}
int comp_event_rec_dchunk_len(int nbeads_) {return 4*nbeads_;}
int comp_event_rec_it_ichunk_len(int nbeads_) {return 0;}
int comp_event_rec_it_dchunk_len(int nbeads_) {return 2*nbeads_;}

event_record::event_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): basic_record(rs_, rid_, ichunk_, dchunk_, u_),
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
  if (take_basic_record(r_))
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

int find_worst_record(event_record ** r_, int ncap_)
{
  int worst_index = 0;
  for (int i = 1; i < ncap_; i++)
    if (r_[worst_index]->isbetter(r_[i]))
      worst_index = i;
  return worst_index;
}

int find_best_record(event_record ** r_, int ncap_)
{
  int best_index = 0;
  for (int i = 1; i < ncap_; i++)
    if (r_[best_index]->isworse(r_[i]))
      best_index = i;
  return best_index;
}

void pick_nworst_records(event_record ** rin_, event_record ** rout_, int n_, int ncap_)
{
  for (int i = 0; i < n_; i++) rout_[i]=rin_[i];
  int i_best_worst = find_best_record(rout_,n_);
  for (int i = n_; i < ncap_; i++)
    if (rout_[i_best_worst]->isbetter(rin_[i]))
    {
      rout_[i_best_worst] = rin_[i];
      i_best_worst=find_best_record(rout_,n_);
    }
}

void pick_nbest_records(event_record ** rin_, event_record ** rout_, int n_, int ncap_)
{
  for (int i = 0; i < n_; i++) rout_[i]=rin_[i];
  int i_worst_best = find_worst_record(rout_,n_);
  for (int i = n_; i < ncap_; i++)
    if (rout_[i_worst_best]->isworse(rin_[i]))
    {
      rout_[i_worst_best] = rin_[i];
      i_worst_best=find_worst_record(rout_,n_);
    }
}

void pick_nworst_records(event_record ** rin_, int n_, int ncap_)
{
  int i_best_worst = find_best_record(rin_,n_);
  for (int i = n_; i < ncap_; i++)
    if (rin_[i_best_worst]->isbetter(rin_[i]))
    {
      rin_[i_best_worst] = rin_[i];
      i_best_worst=find_best_record(rin_,n_);
    }
}

void pick_nbest_records(event_record ** rin_, int n_, int ncap_)
{
  int i_worst_best = find_worst_record(rin_,n_);
  for (int i = n_; i < ncap_; i++)
    if (rin_[i_worst_best]->isworse(rin_[i]))
    {
      rin_[i_worst_best] = rin_[i];
      i_worst_best=find_worst_record(rin_,n_);
    }
}
