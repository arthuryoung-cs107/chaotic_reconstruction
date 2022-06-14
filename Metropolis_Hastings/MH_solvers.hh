#ifndef MH_SOLVERS_HH
#define MH_SOLVERS_HH

#include "MH_workers.hh"

const int genetic_train_const_ilen=1;
const int genetic_train_const_dlen=1;
const int genetic_train_it_ilen=2;
const int genetic_train_it_dlen=2;
class MH_genetic : public basic_MH_trainer
{
  public:

    MH_genetic(MH_train_struct &mhts_, double t_wheels0_, int gen_max_, double alpha_tol_);
    ~MH_genetic();

    void run(bool verbose_=true);
    void stage_diagnostics();
    void close_diagnostics();

    virtual void initialize_run()
    {
      basic_MH_trainer::initialize_run();
      class_count=class_count=0;
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

    const int class_max;

          int class_count,
              event_block_count;

    const double  alpha_tol;

          double  prob_best,
                  prob_worst;

    void inspect_event_data(); 

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
                  ** const alpha_pool_Framebead,
                  ** const mur2_Frame_bead,
                  ** const stdr2_Frame_bead,
                  ** const mualpha_Frame_bead,
                  ** const stdalpha_Frame_bead;

    MH_examiner ** const examiners;

    event_record  ** const records,
                  ** const leaders,
                  ** const pool,
                  ** const leader_board,
                  ** const candidates;

    void clear_event_data();
};

class MH_doctor : public basic_MH_trainer
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
