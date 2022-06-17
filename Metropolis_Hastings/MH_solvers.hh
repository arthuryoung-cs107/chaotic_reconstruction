#ifndef MH_SOLVERS_HH
#define MH_SOLVERS_HH

#include "MH_workers.hh"

const int genetic_train_const_ilen=1;
const int genetic_train_const_dlen=1;
const int genetic_train_it_ilen=2;
const int genetic_train_it_dlen=2;
class MH_genetic : public basic_MH_trainer, public event_block
{
  public:

    MH_genetic(MH_train_struct &mhts_, double t_wheels0_, int gen_max_, double alpha_tol_);
    ~MH_genetic();

    void run(bool verbose_=true);
    void stage_diagnostics();
    void close_diagnostics();
    void report_event_data(event_record **recs_, int n_);

    void consolidate_event_data()
    {
      for (int i = 0; i < nbeads; i++)
      {
        comps_ordered[i]=i;
        stev_ordered[i]=stev_comp[i];
      }
      event_block::consolidate_event_data();
    }
    void define_event_block()
    {event_block::define_event_block(sigma_scaled);}
    void synchronise_event_data()
    {
      int nf_obs, nf_stable, nf_unstable;
      set_state_counts(nf_obs, nf_stable, nf_unstable);
      #pragma omp parallel
      {
        MH_examiner *ex_t=examiners[thread_num()];
        ex_t->synchronise_event_data(&nf_obs, stev_earliest,stev_latest,stev_comp,stev_ordered,comps_ordered,rho2_regime);
      }
    }
    void clear_event_data()
    {
      event_block::clear_event_data();
      #pragma omp parallel
      {
        double * clear_buf;
        #pragma omp for nowait
        for (int i = 0; i < npool; i++)
        {clear_buf=r2_pool_Framebead[i]; for (int j = 0; j < nbeads*Frames; j++) clear_buf[j]=0.0;}

        #pragma omp for nowait
        for (int i = 0; i < npool; i++)
        {clear_buf=alpha_pool_Framebead[i]; for (int j = 0; j < nbeads*Frames; j++) clear_buf[j]=0.0;}
      }
    }
    void initialize_run()
    {
      basic_MH_trainer::initialize_run();
      Class_count=event_block_count=0;
      prob_best=prob_worst=0.0;
    }
    void write_it_ints(FILE * file_)
    {
      MH_trainer::write_it_ints(file_);
      fwrite(genetic_train_it_ints, sizeof(int), genetic_train_it_ilen, file_);
    }
    void write_it_dubs(FILE * file_)
    {
      MH_trainer::write_it_dubs(file_);
      fwrite(genetic_train_it_dubs, sizeof(double), genetic_train_it_dlen, file_);
    }

  protected:

    const int Class_max;

          int Class_count,
              event_block_count;

    const double  alpha_tol;

          double  prob_best,
                  prob_worst;

    int it_ilen_full() {return basic_MH_trainer::it_ilen_full() + genetic_train_it_ilen;}
    int it_dlen_full() {return basic_MH_trainer::it_dlen_full() + genetic_train_it_dlen;}

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

    void post_event_resampling(event_record ** recs_, int n_);
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
