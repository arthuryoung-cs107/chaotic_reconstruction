#include "MH_records.hh"

int comp_event_rec_ichunk_len(int nbeads_) {return nbeads_;}
int comp_event_rec_dchunk_len(int nbeads_) {return 4*nbeads_;}
int comp_event_rec_it_ichunk_len(int nbeads_) {return 0;}
int comp_event_rec_it_dchunk_len(int nbeads_) {return 3*nbeads_;}

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
  fwrite(u,sizeof(double),ulen,file_);
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
