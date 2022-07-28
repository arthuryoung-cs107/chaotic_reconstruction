#ifndef MH_RECORDS_HH
#define MH_RECORDS_HH

#include "MH_learning.hh"

const int basic_rec_ilen=7;
const int basic_rec_dlen=1;
struct basic_record: public record
{
  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_): record(rs_, rid_, ichunk_, dchunk_, u_) {init_basic_record();}

  basic_record(record_struct &rs_, int rid_, int * ichunk_, double * dchunk_, double * u_, MH_rng * ran_, double * umin_, double * umax_): basic_record(rs_, rid_, ichunk_, dchunk_, u_) {draw_ranuni(ran_,umin_,umax_);}

  ~basic_record() {}

  bool success;

  int gen, // 0: generation in which this particle was generated
      Class, // 1: leader Class this particle belongs to, if any
      dup_count, // 2: number of times this particle has been duplicated
      parent_count, // 3: number of particle ancestors
      parent_rid, // 4: index position of parent particle
      parent_gen, // 5: generation this particle comes from
      parent_Class, // 6: leader Class this particle comes from
      * const basic_rec_ints = &gen;

  double  w,
          * const basic_rec_dubs = &w;

  // sampling

  inline void init_basic_record(int gen_=0,int Class_=-1,int dup_count_=0,int parent_count_=0,int parent_rid_=-1,int parent_gen_=-1,int parent_Class_=-1,double w_=0.0)
  {
    gen=gen_; Class=Class_; dup_count=dup_count_;
    parent_count=parent_count_; parent_rid=parent_rid_;
    parent_gen=parent_gen_; parent_Class=parent_Class_;
    w=w_;
  }

  // training

  inline int take_basic_record(basic_record * rec_)
  {
    if (take_record_chunks(rec_)) // successful replacement
    {
      memcpy(basic_rec_ints, rec_->basic_rec_ints, basic_rec_ilen*sizeof(int));
      memcpy(basic_rec_dubs, rec_->basic_rec_dubs, basic_rec_dlen*sizeof(double));
      return 1;
    }
    else return 0;
  }

  inline double get_r2() {return *dub_compare_bad;}

  // io
  inline int basic_rec_ilen_full() {return basic_rec_ilen;}
  inline int basic_rec_dlen_full() {return basic_rec_dlen;}

  inline void write_basic_rec_ints(FILE * file_)
    {fwrite(basic_rec_ints, sizeof(int), basic_rec_ilen, file_);}

  inline void write_basic_rec_dubs(FILE * file_)
    {fwrite(basic_rec_dubs, sizeof(double), basic_rec_dlen, file_);}

  // debugging
  void print_basic_record(const char indent_[]=" ", bool print_all_=true);

};

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
  void print_event_record(int nlead_, int npool_);

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
