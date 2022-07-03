#ifndef MH_RECORDS_HH
#define MH_RECORDS_HH

#include "MH_learning.hh"

int comp_event_rec_ichunk_len(int nbeads_);
int comp_event_rec_dchunk_len(int nbeads_);
int comp_event_rec_it_ichunk_len(int nbeads_);
int comp_event_rec_it_dchunk_len(int nbeads_);
const int event_rec_ilen = 4;
const int event_rec_dlen = 3;
const int event_rec_it_ilen = 0;
const int event_rec_it_dlen = 3;
struct event_record : public basic_record
{
  event_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_);
  event_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double *umin_, double *umax_);
  ~event_record() {}

  int nfobs,
      nfstable,
      nfregime,
      nfunstable,
      * const event_rec_ints = &nfobs;

  double  r2stable,
          r2regime,
          r2unstable,
          * const event_rec_dubs = &r2stable;

  int * evframe_bead;

  double  * r2stable_bead,
          * netr2_regime,
          * r2unstable_bead,
          * alpha_bead;

  // event

  inline void record_event_data(int *int_data_, double * double_data_)
  {
    nfobs=int_data_[0]; nfstable=int_data_[1]; nfregime=int_data_[2]; nfunstable=int_data_[3];
    r2=double_data_[0]; r2stable=double_data_[1]; r2regime=double_data_[2]; r2unstable=double_data_[3];
  }

  void determine_event_block(int &stev_earliest_,int &stev_latest_,int *stev_comp_,int *stev_ordered_,int *comps_ordered_);

  // training

  inline double set_record_stable() {dub_compare_bad = &r2stable; return r2stable;}
  inline double set_record_unstable() {dub_compare_bad = &r2unstable; return r2unstable;}

  inline void record_training_data(double *double_data_, bool success_)
  {
    r2=double_data_[0]; r2stable=double_data_[1]; r2regime=double_data_[2]; r2unstable=double_data_[3];
    success=success_;
  }

  // sampling
  int take_record(event_record *r_);

  // io
  void write_event_rec_full_header(FILE * file_, int len_=0);
  void write_event_rec_training_header(FILE * file_, int len_=0);
  void write_event_rec_training_data(FILE * file_);

  inline void write_event_rec_ints(FILE * file_)
  {
    write_basic_rec_ints(file_);
    fwrite(event_rec_ints,sizeof(int),event_rec_ilen,file_);
  }
  inline void write_event_rec_dubs(FILE * file_)
  {
    write_basic_rec_dubs(file_);
    fwrite(event_rec_dubs,sizeof(double),event_rec_dlen,file_);
  }

  void write_ints(FILE * file_) {write_event_rec_ints(file_);}
  void write_dubs(FILE * file_) {write_event_rec_dubs(file_);}
  inline int event_rec_ilen_full() {return basic_rec_ilen_full() + event_rec_ilen;}
  inline int event_rec_dlen_full() {return basic_rec_dlen_full() + event_rec_dlen;}
  inline int event_rec_ilen_train() {return basic_rec_ilen_full() + event_rec_it_ilen;}
  inline int event_rec_dlen_train() {return basic_rec_dlen_full() + event_rec_it_dlen;}

  private:
    const int event_rec_it_ichunk_len,
              event_rec_it_dchunk_len;
};

int find_worst_record(event_record ** r_, int ncap_);
int find_best_record(event_record **r_, int ncap_);
void pick_nworst_records(event_record ** rin_, event_record ** rout_, int n_, int ncap_);
void pick_nbest_records(event_record ** rin_, event_record ** rout_, int n_, int ncap_);
void pick_nworst_records(event_record ** rin_, int n_, int ncap_);
void pick_nbest_records(event_record ** rin_, int n_, int ncap_);

#endif
