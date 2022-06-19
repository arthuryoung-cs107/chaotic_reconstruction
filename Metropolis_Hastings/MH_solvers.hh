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

    void stage_diagnostics();
    void close_diagnostics();
    void clear_genetic_event_data();
    void consolidate_genetic_event_data();
    void report_genetic_event_data(event_record **recs_, int n_);
    void synchronise_genetic_event_data();
    void post_event_resampling(event_record ** recs_, int n_);

    inline void set_leader_records()
    {
      nreplace=MH_trainer::take_records(leader_board,leaders,irepl_leaders,nlead);
      for (int i = 0; i < nreplace; i++) leaders[i]->init_basic_record(gen_count,Class_count);
    }
    inline void initialize_genetic_run()
    {
      initialize_basic_trainer_run();
      Class_count=event_block_count=0;
      prob_best=prob_worst=0.0;
    }
    inline void set_stable_training()
    {rho2=gaussian_likelihood::expected_r2(stev_comp,nbeads);}
    inline int genetic_it_ilen_full() {return basic_MH_trainer::basic_train_it_ilen_full() + genetic_train_it_ilen;}
    inline int genetic_it_dlen_full() {return basic_MH_trainer::basic_train_it_dlen_full() + genetic_train_it_dlen;}
    inline void write_genetic_it_ints(FILE * file_)
    {
      write_MHT_it_ints(file_);
      fwrite(genetic_train_it_ints, sizeof(int), genetic_train_it_ilen, file_);
    }
    inline void write_genetic_it_dubs(FILE * file_)
    {
      write_MHT_it_dubs(file_);
      fwrite(genetic_train_it_dubs, sizeof(double), genetic_train_it_dlen, file_);
    }

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

class MH_doctor : public basic_MH_trainer, public event_block
{
  public:

    MH_doctor(MH_train_struct &mhts_, int test_id_, int test_relay_id_, int Frames_test_, double alpha_tol_);
    ~MH_doctor();

    void run(bool verbose_=true);
    void stage_diagnostics();
    void close_diagnostics();

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
