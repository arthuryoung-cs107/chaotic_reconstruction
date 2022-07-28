#ifndef MH_SOLVERS_HH
#define MH_SOLVERS_HH

#include "MH_workers.hh"

class basic_MH_trainer: public MH_trainer, public gaussian_likelihood
{
    public:

      basic_MH_trainer(MH_train_struct &mhts_, int ichunk_width_, int dchunk_width_, double t_wheels0_=-1.0);

      ~basic_MH_trainer()
      {
        delete [] ndup_leaders; delete [] irepl_leaders; delete [] isuccess_pool;
        delete [] w_leaders;
        delete [] u_mean; delete [] u_var;
        delete [] u_wmean; delete [] u_wvar;
      }

    protected:

      bool apply_training_wheels;

      const double t_wheels0; // initial drift fraction

      double  t_wheels, // current drift fraction
              r2_scale;

      int * const ndup_leaders,
          * const irepl_leaders,
          * const isuccess_pool;

      double  * const w_leaders,
              * const u_mean,
              * const u_var,
              * const u_wmean,
              * const u_wvar;

      // debugging
      void print_basic_MHT(const char indent_[]="     ");

      // run
      inline void initialize_basic_MHT_run() {initialize_MHT_run(); t_wheels=t_wheels0;}

      // sampling
      virtual void redraw_u(basic_record * rec_pool_, MH_rng * rng_t_)
        {redraw_u_uni(rec_pool_, rng_t_); rec_pool_->init_basic_record(gen_count);}

      void duplicate_u_basic(basic_record *child_, basic_record *parent_, MH_rng *rng_t_, double fac_low_=0.0);

      // training
      inline void clear_basic_MHT_training_data() {memset(isuccess_pool,0,npool*sizeof(int));nsuccess=0;}

      // io
      inline int basic_MHT_it_ilen_full() {return MHT_it_ilen;}
      inline int basic_MHT_it_dlen_full() {return MHT_it_dlen;}
      inline void write_basic_MHT_it_ints(FILE * file_) {write_MHT_it_ints(file_);}
      inline void write_basic_MHT_it_dubs(FILE * file_) {write_MHT_it_dubs(file_);}
      inline void write_ustats(FILE *file_)
      {
        fwrite(u_mean,sizeof(double),ulen,file_);
        fwrite(u_wmean,sizeof(double),ulen,file_);
        fwrite(u_var,sizeof(double),ulen,file_);
        fwrite(u_wvar,sizeof(double),ulen,file_);
      }
};

const int genetic_train_const_ilen=3;
const int genetic_train_const_dlen=3;
const int genetic_train_it_ilen=4;
const int genetic_train_it_dlen=2;
class MH_genetic : public basic_MH_trainer, public event_block
{
  public:

    MH_genetic(MH_train_struct &mhts_, int Class_max_, int gen_max_, int itrain_max_, double t_wheels0_, double alpha_tol_, double rs_full_factor_, double train_tol_);
    ~MH_genetic();

    void run(bool verbose_=true);

    bool write_full_training_data=false;

  protected:

    size_t  obuf_end;

    char * const obuf;

    const int Class_max,
              gen_max,
              itrain_max,
              * const genetic_train_const_ints = &Class_max;

          int Class_count,
              event_block_count,
              write_training_pool,
              write_pool_event_data,
              * const genetic_train_it_ints = &Class_count;

    const double  alpha_tol,
                  rs_full_factor,
                  train_tol,
                  * const genetic_train_const_dubs = &alpha_tol;

          double  prob_best,
                  prob_worst,
                  * const genetic_train_it_dubs = &prob_best,
                  ** const r2_pool_Framebead,
                  ** const alpha_pool_Framebead;

    MH_examiner ** const examiners;

    event_record  * bleader,
                  * wleader,
                  * bpool,
                  ** const records,
                  ** const leaders,
                  ** const pool,
                  ** const leader_board,
                  ** const candidates;


    // debugging
    void print_MH_genetic();

    // run
    inline void initialize_genetic_run()
    {
      initialize_basic_MHT_run();
      initialize_event_block();
      Class_count=event_block_count=0;
      prob_best=prob_worst=0.0;

      write_training_pool=(int) write_full_training_data;
      write_pool_event_data=(int) write_full_training_data;
    }

    void find_events(bool verbose_);
    void train_event_block(bool verbose_, bool &stable_convergence_, bool &unstable_convergence_);
    bool check_run_convergence(bool stable_conv_, bool unstable_conv_);

    // event

    void clear_genetic_event_data();
    inline void consolidate_genetic_event_data()
    {
      for (int i = 0; i < nbeads; i++)
      {
        comps_ordered[i]=i;
        stev_ordered[i]=stev_comp[i]=max(stev_comp[i],stev_comp_old[i]); // ensure we have strictly increasing events
      }
      event_block::consolidate_event_data(); // sort event states
    }

    inline void define_genetic_event_block()
      {event_block::define_event_block(sigma_scaled);}

    void synchronise_genetic_event_data();
    void report_genetic_event_data();
    inline void stage_event_search()
    {
      for (int i = 0; i < nbeads; i++) stev_comp_old[i] = stev_comp[i];
      stev_min=stev_ordered[nbeads-1]; // ensuring that in our next event search, we evaluate up to the last frame of previous event block
      take_records(leaders,pool,nlead);
    }

    // training
    double set_objective(bool verbose_, double &r2_scale_, bool stable_flag_);
    bool train_objective(bool verbose_,int &nit_,int &nit_objective_,double rho2_);
    inline void clear_genetic_training_data() {clear_basic_MHT_training_data();}
    double consolidate_genetic_training_data(double wsum_pool_,double *w_leaders_,double rho2_,int &nreplace_,double &r2_scale_);

    inline void report_genetic_training_data(int nreplace_,int &Class_count_,int &gen_count_)
    {
      if (nreplace_) write_Class_diagnostics(Class_count_++);
      write_generation_diagnostics(gen_count_++);
    }

    bool check_objective_convergence(int nit_, int nit_objective_, bool &training_success_);

    // sampling
    double set_leader_records();
    void respawn_pool(bool verbose_, double wsum_, double *w_leaders_, int nreload_=0);
    double compute_weights(double r2_min_, double rho2in_, event_record ** recs_, int n_);
    double compute_weights(double r2_min_, double rho2in_, double *w_leaders_);
    void compute_weighted_ustats(double wsum_, event_record ** recs_, int n_);

    inline void take_records(event_record ** rin_, event_record ** rout_, int ncap_)
      {for (int i = 0; i < ncap_; i++) int replace_status=rout_[i]->take_record(rin_[i]);}

    inline int take_records(event_record ** rin_, event_record ** rout_, int *repl_list_, int ncap_)
    {
      int nrepl=0;
      for (int i = 0; i < ncap_; i++) if (rout_[i]->take_record(rin_[i])) repl_list_[nrepl++]=i;
      return nrepl;
    }

    inline void duplicate_u(event_record *child_, event_record *parent_, MH_rng *rng_t_) {duplicate_u_basic(child_,parent_,rng_t_,0.1);}

    // io
    void stage_diagnostics();
    void close_diagnostics();
    void write_event_diagnostics(int event_block_count_);
    void write_Class_diagnostics(int Class_count_);
    void write_generation_diagnostics(int gen_count_);

    inline void write_genetic_it_ints(FILE * file_)
    {
      write_basic_MHT_it_ints(file_);
      fwrite(genetic_train_it_ints,sizeof(int),genetic_train_it_ilen,file_);
    }

    inline void write_genetic_it_dubs(FILE * file_)
    {
      write_basic_MHT_it_dubs(file_);
      fwrite(genetic_train_it_dubs,sizeof(double),genetic_train_it_dlen,file_);
    }

    inline int genetic_it_ilen_full()
      {return basic_MHT_it_ilen_full() + genetic_train_it_ilen;}

    inline int genetic_it_dlen_full()
      {return basic_MHT_it_dlen_full() + genetic_train_it_dlen;}

  private:

    // verbose
    void verbose_find_events_1();
    void verbose_find_events_2();
    void verbose_find_events_3();
    void verbose_set_objective_1();
    void verbose_set_objective_2();
    void verbose_train_objective_1(int nit_, int nsuccess_local_);
    void verbose_train_objective_2();
    void verbose_respawn_pool(int offset_);
};

class MH_doctor : public MH_genetic
{
  public:

    MH_doctor(MH_train_struct &mhts_, int test_id_, int test_relay_id_, int Frames_test_, double alpha_tol_);
    ~MH_doctor();

    void run(bool verbose_=true);

  protected:

    void initialize_doctor_run();
    void stage_doctor_diagnostics();
    inline void clear_doctor_event_data() {clear_genetic_event_data();}
    inline void consolidate_doctor_event_data() {consolidate_genetic_event_data();}
    inline void define_doctor_event_block() {define_genetic_event_block();}
    inline void synchronise_doctor_event_data() {synchronise_genetic_event_data();}
    inline void report_doctor_event_data() {report_genetic_event_data();}

  private:

    size_t test_buf_end;

    const int test_id,
              test_relay_id,
              Frames_test;

    char * const test_buffer;
    double  * const TEST_refp;

    MH_medic ** const medics;
};

#endif
