#ifndef MH_SOLVERS_HH
#define MH_SOLVERS_HH

#include "MH_workers.hh"

const int genetic_train_const_ilen=2;
const int genetic_train_const_dlen=2;
const int genetic_train_it_ilen=2;
const int genetic_train_it_dlen=2;
class MH_genetic : public basic_MH_trainer, public event_block
{
  public:

    MH_genetic(MH_train_struct &mhts_, int Class_max_, int gen_max_, double t_wheels0_, double alpha_tol_, double rs_full_factor_);
    ~MH_genetic();

    void run(bool verbose_=true);

  protected:

    const int Class_max,
              gen_max;

          int Class_count,
              event_block_count;

    const double  alpha_tol,
                  rs_full_factor;

          double  prob_best,
                  prob_worst;

    // run
    bool check_convergence();

    inline void initialize_genetic_run()
    {
      initialize_basic_trainer_run();
      Class_count=event_block_count=0;
      prob_best=prob_worst=0.0;
    }

    // event
    void find_events(bool verbose_=true);
    void clear_genetic_event_data();
    void consolidate_genetic_event_data();
    void consolidate_genetic_event_data(int bleader_index_);
    void synchronise_genetic_event_data();
    void report_genetic_event_data();
    void stage_event_search();

    // training
    void train_event_block(bool verbose_=true);
    void set_regime_objective(int iregime_);
    inline void clear_genetic_training_data() {clear_basic_train_training_data();}
    bool check_regime_convergence(int iregime_count_);
    void consolidate_genetic_training_data();
    inline void report_genetic_training_data()
    {
      if (nreplace) write_Class_diagnostics(Class_count);
      write_generation_diagnostics(gen_count);
    }
    void set_leader_records();

    // sampling
    void reload_leaders(int bleader_index_);
    void post_event_resampling();
    void respawn_pool(double wsum_, int offset_=0);

    // inelegant, but I kind of screwed myself over by not planning around this polymorphism issue
    virtual int find_worst_record(event_record ** r_, int ncap_);
    virtual int find_best_record(event_record **r_, int ncap_);
    virtual void pick_nworst_records(event_record ** rin_, event_record ** rout_, int n_, int ncap_);
    virtual void pick_nbest_records(event_record ** rin_, event_record ** rout_, int n_, int ncap_);
    virtual void pick_nworst_records(event_record ** rin_, int n_, int ncap_);
    virtual void pick_nbest_records(event_record ** rin_, int n_, int ncap_);
    virtual double compute_weights(double r2_min_, event_record ** recs_, int n_);
    inline void take_records(event_record ** rin_, event_record ** rout_, int ncap_)
    {for (int i = 0; i < ncap_; i++) rout_[i]->take_record(rin_[i]);}
    inline int take_records(event_record ** rin_, event_record ** rout_, int *repl_list_, int ncap_)
    {
      int nrepl=0;
      for (int i = 0; i < ncap_; i++) if (rout_[i]->take_record(rin_[i])) repl_list_[nrepl++]=i;
      return nrepl;
    }

    // io
    void stage_diagnostics();
    void close_diagnostics();
    void write_Class_diagnostics(int &Class_count_);
    void write_generation_diagnostics(int &gen_count_);
    inline void write_genetic_it_ints(FILE * file_)
    {write_basic_train_it_ints(file_);fwrite(genetic_train_it_ints,sizeof(int),genetic_train_it_ilen,file_);}
    inline void write_genetic_it_dubs(FILE * file_)
    {write_basic_train_it_dubs(file_);fwrite(genetic_train_it_dubs,sizeof(double),genetic_train_it_dlen,file_);}
    inline void write_ustats(FILE * file_)
    {

    }
    inline int genetic_it_ilen_full() {return basic_train_it_ilen_full() + genetic_train_it_ilen;}
    inline int genetic_it_dlen_full() {return basic_train_it_dlen_full() + genetic_train_it_dlen;}

  private:

    size_t  obuf_end;

    char * const obuf;

    const int * const genetic_train_const_ints;
          int * const genetic_train_it_ints;

    const double  * const genetic_train_const_dubs;
          double  * const genetic_train_it_dubs,
                  ** const r2_pool_Framebead,
                  ** const alpha_pool_Framebead;

    MH_examiner ** const examiners;

    event_record  ** const records,
                  ** const leaders,
                  ** const pool,
                  ** const leader_board,
                  ** const candidates;
};

class MH_doctor : public MH_genetic
{
  public:

    MH_doctor(MH_train_struct &mhts_, int test_id_, int test_relay_id_, int Frames_test_, double alpha_tol_);
    ~MH_doctor();

    void run(bool verbose_=true);
    void stage_doctor_diagnostics();
    void close_doctor_diagnostics();

    inline void initialize_doctor_run() {initialize_genetic_run();}

  private:

    size_t test_buf_end;

    const int test_id,
              test_relay_id,
              Frames_test;

    const double alpha_tol;

    char * const test_buffer;
    double  * const TEST_refp;

    MH_medic ** const medics;
    event_record  ** const records,
                  ** const leaders,
                  ** const pool;
};

#endif
