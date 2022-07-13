#ifndef MH_SOLVERS_HH
#define MH_SOLVERS_HH

#include "MH_workers.hh"

const int genetic_train_const_ilen=3;
const int genetic_train_const_dlen=3;
const int genetic_train_it_ilen=2;
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


    // debugging write outs

    inline void print_MH_genetic()
    {
      printf("(MH_genetic) Class_count: %d, event_block_count %d\n", Class_count,event_block_count);
      printf("(MH_genetic) prob_best: %e, prob_worst %e\n", prob_best,prob_worst);
      print_basic_MHT("     ");
      print_event_block(sigma_scaled,"     ");
    }

    // run
    inline void initialize_genetic_run()
    {
      initialize_basic_MHT_run();
      Class_count=event_block_count=0;
      prob_best=prob_worst=0.0;
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
        stev_ordered[i]=stev_comp[i];
      }
      event_block::consolidate_event_data();
    }

    inline void define_genetic_event_block()
      {event_block::define_event_block(sigma_scaled);}

    void synchronise_genetic_event_data();
    void report_genetic_event_data();
    inline void stage_event_search() {take_records(leaders,pool,nlead);}

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

    inline void verbose_find_events_1()
      {printf("(MH_genetic::find_events) Found event block %d with generation %d.", event_block_count,gen_count);}

    inline void verbose_find_events_2() {printf(" Earliest, latest frames: %d, %d\n", stev_earliest, stev_latest);}
    inline void verbose_find_events_3()
      {printf("(bead,frame), chronologically:\n");
      for (int i = 0; i < nbeads; i++) printf("(%d %d)\n",comps_ordered[i],stev_ordered[i]);}

    inline void verbose_set_objective_1()
    {
      printf("(MH_genetic::set_objective) %d %s frames, rho2 = %e. ",
      (stable_flag)?event_block::nst_stable():event_block::nst_unstable(),
      (stable_flag)?"STABLE":"UNSTABLE",
      rho2_objective);
    }

    inline void verbose_set_objective_2()
      {printf("nreplace: %d, r2 best (%d): %e, r2 worst (%d): %e.\n",nreplace,bleader_rid,br2,wleader_rid,wr2);}

    inline void verbose_train_objective_1(int nit_)
      {printf("(MH_genetic::train_objective) gen %d, Class %d, nit %d: %d candidates, best candidate (%d) r2 = %e. ", gen_count, Class_count, nit_, nsuccess, bpool->rid, bpool->get_r2());}

    inline void verbose_train_objective_2()
      {printf("%d replacements, r2 best (%d): %e, r2 worst (%d): %e. ", nreplace, bleader_rid, br2, wleader_rid, wr2);}

    inline void verbose_respawn_pool(int offset_)
    {
      // if (offset_) printf("(MH_genetic::respawn_pool) %d reloads, ", offset_);
      // else printf("(MH_genetic::respawn_pool) ");
      if (offset_) printf("RESPAWNING: %d reloads, ", offset_);
      else printf("RESPAWNING: ");
      printf("%d redraws, %d duplicates (%d unique).\n", nredraw, ndup, ndup_unique);
    }
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
