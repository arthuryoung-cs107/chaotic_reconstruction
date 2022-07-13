#ifndef MH_RECORDS_HH
#define MH_RECORDS_HH

#include "MH_learning.hh"

int comp_event_rec_ichunk_len(int nbeads_);
int comp_event_rec_dchunk_len(int nbeads_);
int comp_event_rec_it_ichunk_len(int nbeads_);
int comp_event_rec_it_dchunk_len(int nbeads_);
const int event_rec_ilen = 3;
const int event_rec_dlen = 3;
const int event_rec_it_ilen = 0;
const int event_rec_it_dlen = 3;
struct event_record : public basic_record
{
  event_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): basic_record(rs_, rid_, ichunk_, dchunk_, u_),
  event_rec_it_ichunk_len(comp_event_rec_it_ichunk_len(nbeads)), event_rec_it_dchunk_len(comp_event_rec_it_dchunk_len(nbeads))
  {}

  event_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double *umin_, double *umax_): basic_record(rs_, rid_, ichunk_, dchunk_, u_, ran_, umin_, umax_),
  event_rec_it_ichunk_len(comp_event_rec_it_ichunk_len(nbeads)), event_rec_it_dchunk_len(comp_event_rec_it_dchunk_len(nbeads))
  {}

  ~event_record() {}

  int nfobs,
      nfstable,
      nfunstable,
      * const event_rec_ints = &nfobs;

  double  r2net,
          r2stable,
          r2unstable,
          * const event_rec_dubs = &r2net;

  int * evframe_bead=ichunk;

  double  * const r2net_bead=dchunk,
          * const r2stable_bead=dchunk+nbeads,
          * const r2unstable_bead=dchunk+2*nbeads,
          * const alpha_bead=dchunk+3*nbeads;

  // event

  inline void record_event_data(int *int_data_, double * double_data_)
  {
    nfobs=int_data_[0]; nfstable=int_data_[1]; nfunstable=int_data_[2];
    r2net=double_data_[0]; r2stable=double_data_[1]; r2unstable=double_data_[2];
  }

  // training

  inline double set_record_net() {dub_compare_bad = &r2net; return r2net;}
  inline double set_record_stable() {dub_compare_bad = &r2stable; return r2stable;}
  inline double set_record_unstable() {dub_compare_bad = &r2unstable; return r2unstable;}
  inline void clear_residuals() {for (int i = 0; i < 3*nbeads; i++) dchunk[i]=0.0;}

  inline void record_training_data(double *double_data_, bool success_)
  {
    r2net=double_data_[0]; r2stable=double_data_[1]; r2unstable=double_data_[2];
    success=success_;
  }

  int take_record(event_record *r_);

  // sampling

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

  // debugging

  inline void print_event_record(int nlead_, int npool_)
  {
    printf("(event_record) record %d of %d leaders, %d pool.\n", rid,nlead_,npool_);
    printf("(event_record) int data:\n");
    printf(" nfobs: %d\n",nfobs);
    printf(" nfstable: %d\n",nfstable);
    printf(" nfunstable: %d\n",nfunstable);

    printf("(event_record) double data:\n");
    printf(" r2net: %e\n",r2net);
    printf(" r2stable: %e\n",r2stable);
    printf(" r2unstable: %e\n",r2unstable);

    printf("(event_record) evframe_bead: "); print_row_vec(evframe_bead,nbeads);
    printf("(event_record) r2stable_bead: "); print_row_vec(r2stable_bead,nbeads);
    printf("(event_record) r2unstable_bead: "); print_row_vec(r2unstable_bead,nbeads);
    printf("(event_record) alpha_bead: "); print_row_vec(alpha_bead,nbeads);

    print_basic_record();
  }

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
