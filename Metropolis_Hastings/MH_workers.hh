#ifndef MH_WORKERS_HH
#define MH_WORKERS_HH

#include "MH_records.hh"


class MH_examiner: public basic_thread_worker, public event_block
{
  public:

    MH_examiner(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_);
    ~MH_examiner();

    void update_event_data(double *r2i_, double alphai_);
    void detect_events(event_record *rec_, double *r2i_, double *alphai_);
    void consolidate_results();
    bool report_results(bool first2finish_);

    inline void clear_event_data()
    {
      test_count=0;
      for (int i = 0; i < nbeads*Frames; i++)
      {
        evcount_comp_state[0][i]=nobs_bead_Frame[0][i]=0;
        mur2_Frame_bead[0][i]=sigmar2_Frame_bead[0][i]=
        mualpha_Frame_bead[0][i]=sigmaalpha_Frame_bead[0][i]=0.0;
      }
    }

  private:
    int test_count;

    int ** const nobs_bead_Frame;

    double  ** const mur2_Frame_bead,
            ** const stdr2_Frame_bead,
            ** const mualpha_Frame_bead,
            ** const stdalpha_Frame_bead;

    void start_detecting_events(event_record * rec_, double * t_history_ double &net_r2_local_);
};

class MH_medic: public basic_thread_worker
{
  public:

    MH_medic(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_, int Frames_test_, char * test_directory_);
    ~MH_medic();

    void test_u(event_record *rec_, int i_, bool verbose_);
    bool report_results(bool first2finish_, int ** evcount_comp_state_agg_)

    inline void initialize_utest() {for (int i = 0; i < nbeads*Frames; i++) evcount_comp_state[0][i] = 0;}

  private:

    const size_t buf_end;

    const int Frames_test;

    char * const mtest_buffer;

    double  * const TEST_p,
            * const TEST_r2,
            * const TEST_alpha,
            * const TEST_INTr2;

    int start_test_u(event_record *rec_, double *t_history_, double &net_r2_local_);
    void write_utest_results(event_record *rec_, int i_);
};
#endif
