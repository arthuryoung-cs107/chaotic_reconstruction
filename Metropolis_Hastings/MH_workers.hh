#ifndef MH_WORKERS_HH
#define MH_WORKERS_HH

#include "MH_records.hh"

class MH_examiner: public basic_thread_worker
{
  public:

    MH_examiner(swirl_param &sp_, proximity_grid *pg_, wall_list &wl_, thread_worker_struct &tws_, int thread_id_, double alpha_tol_): basic_thread_worker(sp_, pg_, wl_, tws_, thread_id_, alpha_tol_, comp_iwkspc_len(tws_),comp_dwkspc_len(tws_)) {}
    ~MH_examiner() {}

    int * const dupcount_leaders=int_wkspc,
        * const itest_list=int_wkspc,
        * const isuccess_list=int_wkspc+npool;

    double * const ustat_buf=dub_wkspc;

    // event
    inline void clear_examiner_event_data() {event_block::clear_event_data(); ntest=0;}
    void detect_events(event_record *rec_, double *r2i_, double *alphai_, int stev_min_=0);
    void consolidate_examiner_event_data();
    bool report_examiner_event_data(bool first2finish_, int &stev_earliest_, int &stev_latest_, int *stev_c_, int ** nev_s_c_, int ** nobs_s_c_, double ** r2_s_c_, double **alpha_s_c_);

    inline void synchronise_examiner_event_data(int stev_e_, int stev_l_, int *stev_c_, int *stev_o_, int *comps_o_,double *rho2s_c_,double *rho2us_c_)
    {
      event_block::synchronise_event_data(stev_e_,stev_l_,stev_c_,stev_o_,comps_o_,rho2s_c_,rho2us_c_);
      nf_netobs=event_block::nst_full();
      nf_stobs=event_block::nst_stable();
      nf_usobs=event_block::nst_unstable();
    }

    void restore_event_record(event_record *rec_, double *r2_Fb_, double *alpha_Fb_);

    // training

    inline void clear_examiner_training_data() {clear_basic_tw_training_data(); ntest=0; nsuccess_test=0;}
    bool examine_u(event_record *pooli_, int i_, double r2success_threshold_);
    void consolidate_examiner_training_data(event_record ** pool_);

    bool report_examiner_training_data(bool first2finish_, event_record ** bpool_address_, int *isuccess_pool_,int &nsuccess_,double *u_wmean_);

    // debugging

    inline void print_MH_examiner(int thread_count_, double sigma_scaled_)
    {
      printf("(MH_examiner) worker %d of %d:\n", thread_id, thread_count_);
      print_basic_tw("     ",sigma_scaled_);
    }


  protected:
    friend class MH_medic;

    int ntest,
        nsuccess_test,
        nf_netobs,  // net number of observation frames
        nf_stobs,   // net number of stable observation frames
        nf_usobs;   // net number of unstable observation frames

    event_record *btest;

    // event
    void start_detecting_events(int &f_local_,int *f_event_,double &netr2_local_,double &netr2_stable_local_,double &netr2_unstable_local_,double *t_history_,double *r2net_bead,double *r2stable_bead_,double *r2unstable_bead_,double *alphaev_bead_);
    void update_event_data(int final_frame_, int *f_event_, double *r2i_, double *alphai_);

    // training
    inline bool update_training_data(int i_, double r2success_threshold_)
    {
      itest_list[ntest++]=i_;
      if (basic_thread_worker::check_success(r2success_threshold_))
        {isuccess_list[nsuccess_test++] = i_; return true;}
      else return false;
    }

  private:

    inline int comp_iwkspc_len(thread_worker_struct &tws_) {return 2*tws_.npool;}
    inline int comp_dwkspc_len(thread_worker_struct &tws_) {return tws_.ulen;}
};

class MH_medic
{
  public:
    MH_medic(MH_examiner &ex_, int Frames_test_, char * test_buffer_, size_t test_buf_end_);
    ~MH_medic();

    MH_examiner &ex;

    inline void clear_medic_event_data()
    {ex.clear_examiner_event_data(); stev_earliest=ex.stev_earliest; stev_latest=ex.stev_latest; ntest=ex.ntest;}

    void test_u(event_record *rec_, int i_, double *r2i_, double *alphai_, bool verbose_);

    inline void consolidate_medic_event_data()
      {ex.consolidate_examiner_event_data();}

    inline bool report_medic_event_data(bool first2finish_, int &stev_earliest_, int &stev_latest_, int *stev_c_, int ** nev_s_c_, int ** nobs_s_c_, double ** r2_s_c_, double **alpha_s_c_)
      {return ex.report_examiner_event_data(first2finish_,stev_earliest_,stev_latest_,stev_c_,nev_s_c_,nobs_s_c_,r2_s_c_,alpha_s_c_);}

  protected:

    // matching parameters found in MH_examiner

    const int nbeads=ex.nbeads,
              dim=ex.dim,
              dof=ex.dof,
              thread_id=ex.thread_id;

    const double  t_phys=ex.t_phys,
                  alpha_tol=ex.alpha_tol;

    int * const stev_comp=ex.stev_comp,
        * const stev_ordered=ex.stev_ordered,
        * const comps_ordered=ex.comps_ordered,
        ** const nev_state_comp=ex.nev_state_comp,
        ** const nobs_state_comp=ex.nobs_state_comp;

    double  * const ts=ex.ts,
            * const d_ang=ex.d_ang,
            * const comega_s=ex.comega_s,
            * const xs=ex.xs,
            * const psim=ex.psim,
            * const rho2stable_comp=ex.rho2stable_comp,
            * const rho2unstable_comp=ex.rho2unstable_comp,
            ** const r2_state_comp=ex.r2_state_comp,
            ** const alpha_state_comp=ex.alpha_state_comp,
            ** const INTr2_comp_history=ex.INTr2_comp_history;

    particle * const q=ex.q;

  private:

    const size_t buf_end;

    const int Frames_test;

    char * const mtest_buffer;

    int ntest,
        stev_early,
        stev_late,
        stev_earliest,
        stev_latest,
        nf_netobs,
        nf_stobs,
        nf_usobs;

    double  netr2,
            netr2_stable,
            netr2_unstable,
            rho2stable,
            * const TEST_p,
            * const TEST_INTr2;

    void start_test_u(int &f_local_,int *f_event_,double &netr2_local_,double &netr2_stable_local_,double &netr2_unstable_local_,double *t_history_,double *r2net_bead_,double *r2stable_bead_,double *r2unstable_bead_,double *alphaev_bead_);

    inline void update_event_data(int f_local_, int *f_event_, double *r2i_, double *alphai_)
    {
      ex.update_event_data(f_local_,f_event_,r2i_,alphai_);

      // update by referencing this medic's examiner
      ntest=ex.ntest;
      stev_early=ex.stev_early;
      stev_late=ex.stev_late;
      stev_earliest=ex.stev_earliest;
      stev_latest=ex.stev_latest;
      nf_netobs=ex.nf_netobs;
      nf_stobs=ex.nf_stobs;
      nf_usobs=ex.nf_usobs;
    }

    void write_utest_results(event_record *rec_, int i_);

    // wrapper functions

    inline void reset_sim(double *utest_, double t0_, double ctheta0_, double comega0_, double *p0_)
      {ex.reset_sim(utest_,t0_,ctheta0_,comega0_,p0_);}

    inline double * advance_sim(int f_local_,double *t_history_)
      {return ex.advance_sim(f_local_,t_history_);}

    inline double compute_residual(double xs_, double ys_, double &x_now_, double &y_now_, double xr_, double yr_)
      {return ex.compute_residual(xs_,ys_,x_now_,y_now_,xr_,yr_);}

    inline void update_integral_history(double INT_now_, int ibead_)
      {ex.update_integral_history(INT_now_,ibead_);}

    inline double alpha_comp(double *a_, double t_p1, double t_m1)
      {return ex.alpha_comp(a_,t_p1,t_m1);}
};
#endif
