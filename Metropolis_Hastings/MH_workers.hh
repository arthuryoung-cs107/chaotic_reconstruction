#ifndef MH_WORKERS_HH
#define MH_WORKERS_HH

#include "MH_records.hh"

class MH_examiner: public basic_thread_worker, public event_block
{
  public:

    MH_examiner(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_): basic_thread_worker(sp_, pg_, wl_, tws_, thread_id_, alpha_tol_), event_block(nbeads, Frames) {}
    ~MH_examiner() {}

    void detect_events(event_record *rec_, double *r2i_, double *alphai_);
    void consolidate_examiner_event_data();

    inline void synchronise_examiner_event_data(int *nf_, int stev_earliest_, int stev_latest_, int *stev_c_, int *stev_o_, int *comps_o_, double *rho2_r_)
    {nf_obs=nf_[0]; nf_stable=nf_[1]; nf_unstable=nf_[2];
    event_block::synchronise_event_data(stev_earliest_,stev_latest_,stev_c_,stev_o_,comps_o_,rho2_r_);}
    inline void clear_examiner_event_data() {event_block::clear_event_data(); test_count=0;}

  protected:

    int test_count;

    void start_detecting_events(event_record * rec_, double * t_history_ double &net_r2_local_);
    void update_event_data(int final_frame_, double *r2i_, double *alphai_);
};

class MH_medic: public basic_thread_worker, public event_block
{
  public:

    MH_medic(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_, int Frames_test_, char * test_directory_);
    ~MH_medic();

    void test_u(event_record *rec_, int i_, bool verbose_);
    bool report_results(bool first2finish_, int ** nev_state_comp_)

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
