#ifndef MH_WORKERS_HH
#define MH_WORKERS_HH

#include "MH_records.hh"

class MH_examiner: public basic_thread_worker
{
  public:

    MH_examiner(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_): basic_thread_worker(sp_, pg_, wl_, tws_, thread_id_, alpha_tol_) {}
    ~MH_examiner() {}

    // event
    inline void clear_examiner_event_data() {event_block::clear_event_data(); ntest=0;}
    void detect_events(event_record *rec_, double *r2i_, double *alphai_);
    void consolidate_examiner_event_data();
    inline void synchronise_examiner_event_data(int *nf_, int stev_earliest_, int stev_latest_, int *stev_c_, int *stev_o_, int *comps_o_,double *rho2s_c_, double *drho2_r_)
    {event_block::synchronise_event_data(stev_earliest_,stev_latest_,stev_c_,stev_o_,comps_o_,rho2s_c_,drho2_r_); nf_obs=nf_[0]; nf_stable=nf_[1]; nf_unstable=nf_[2];}
    bool report_examiner_event_data(bool first2finish_, int &stev_earliest_, int &stev_latest_, int *stev_c_, int ** nev_s_c_, int ** nobs_s_c_, double ** r2_s_c_, double **alpha_s_c_);
    void restore_event_record(event_record *rec_, double *r2_Fb_);

    // training
    inline void clear_examiner_training_data() {clear_basic_tw_training_data(); ntest=0; nsuccess_test=0;}
    bool examine_u(event_record *pooli_, int i_, double r2success_threshold_);
    void consolidate_examiner_training_data();
    inline bool report_examiner_training_data(bool first2finish_,int *isuccess_pool_,int &nsuccess_)
    {
      for (int i = 0; i < nsuccess_test; i++) isuccess_pool_[i+nsuccess_]=int_wkspc[i];
      nsuccess_+=nsuccess_test;
      return false;
    }
    void set_regime_objective(int iregime_);
    void set_record_regime(event_record *rec_);

  protected:

    int ntest,
        nsuccess_test;

    // event
    void start_detecting_events(int &f_local_,int &iregime_local_,int *f_event_,double &netr2_local_,double &netr2_stable_local_,double &netr2_unstable_local_,double *t_history_,double *r2stable_bead_,double *netr2_regime_,double *r2unstable_bead_,double *alphaev_bead_);
    void update_event_data(int final_frame_, int *f_event_, double *r2i_, double *alphai_);

    // training
    bool update_training_data(int i_, double r2success_threshold_);
    inline void clear_examiner_residuals(double *r2s_c_, double *r2n_r_,double *r2us_c_)
    {for (int i = 0; i < nbeads; i++) r2_state_comp[0][i]=r2s_c_[i]=r2n_r_[i]=r2us_c_[i]=0.0;}

};

class MH_medic
{
  public:
    MH_medic(MH_examiner *ex_, int Frames_test_, char * test_buffer_);
    ~MH_medic();

    MH_examiner * const ex;

    void test_u(event_record *rec_, int i_, bool verbose_);
    bool report_results(bool first2finish_, int ** nev_state_comp_);

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
