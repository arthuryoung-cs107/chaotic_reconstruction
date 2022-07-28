#include "MH_solvers.hh"

/*
  -----------------------------------------------------------------
  -----------------------------------------------------------------
  MH_tools
  -----------------------------------------------------------------
  -----------------------------------------------------------------
*/

// event_block -----------------------------------------------------------------
void event_block::print_event_block(double sigma_scaled_, const char indent_[])
{
  printf("%s(event_block) stev_earliest: %d, stev_latest: %d, nst_full: %d, nst_stable: %d, nst_unstable: %d\n",
  indent_,
  stev_earliest,
  stev_latest,
  nst_full(),
  nst_stable(),
  nst_unstable());
  printf("%s(event_block) rho2_objective: %e, netrho2: %e, rho2stable: %e, rho2unstable: %e\n",
  indent_,
  rho2_objective,
  compute_netrho2(),
  compute_netrho2stable(),
  compute_netrho2unstable());
  printf("%s(event_block) stev_comp: ",indent_); print_row_vec(stev_comp, ncomp);
  printf("%s(event_block) stev_ordered: ",indent_); print_row_vec(stev_ordered, ncomp);
  printf("%s(event_block) rho2stable_comp: ",indent_); print_row_vec(rho2stable_comp, ncomp);
  printf("%s(event_block) rho2unstable_comp: ",indent_); print_row_vec(rho2unstable_comp, ncomp);
}

// gaussian_likelihood -----------------------------------------------------------------
void gaussian_likelihood::print_gaussian_likelihood(const char indent_[])
  {printf("%s(gaussian_likelihood) sigma_noise = %e, noise_scale = %e, sigma_scaled = %e\n", indent_,sigma_noise,noise_scale,sigma_scaled);}


/*
  -----------------------------------------------------------------
  -----------------------------------------------------------------
  MH_learning
  -----------------------------------------------------------------
  -----------------------------------------------------------------
*/

// record -----------------------------------------------------------------
void record::print_record(const char indent_[])
  {printf("%s(record): ichunk_len=%d, dchunk_len=%d\n", indent_, ichunk_len,dchunk_len);}

// MH_trainer -----------------------------------------------------------------
void MH_trainer::print_MHT(const char indent_[])
{
  printf("%s(MH_trainer) leader_count: %d\n", indent_,leader_count);
  printf("%s(MH_trainer) gen_count: %d\n", indent_,gen_count);
  printf("%s(MH_trainer) nsuccess: %d\n", indent_,nsuccess);
  printf("%s(MH_trainer) ncandidates: %d\n", indent_,ncandidates);
  printf("%s(MH_trainer) bleader_rid: %d\n", indent_,bleader_rid);
  printf("%s(MH_trainer) wleader_rid: %d\n", indent_,wleader_rid);
  printf("%s(MH_trainer) nreplace: %d\n", indent_,nreplace);
  printf("%s(MH_trainer) ndup: %d\n", indent_,ndup);
  printf("%s(MH_trainer) ndup_unique: %d\n", indent_,ndup_unique);
  printf("%s(MH_trainer) nredraw: %d\n", indent_,nredraw);

  printf("%s(MH_trainer) rho2: %e\n",indent_,rho2);
  printf("%s(MH_trainer) br2: %e\n",indent_,br2);
  printf("%s(MH_trainer) wr2: %e\n",indent_,wr2);
}


/*
  -----------------------------------------------------------------
  -----------------------------------------------------------------
  FULL IMPLEMENTATIONS
  -----------------------------------------------------------------
  -----------------------------------------------------------------
*/

// records -----------------------------------------------------------------
void basic_record::print_basic_record(const char indent_[], bool print_all_)
{
  printf("%s(basic_record) int data:\n", indent_);
  printf("%s gen: %d\n",indent_,gen);
  printf("%s Class: %d\n",indent_,Class);
  printf("%s dup_count: %d\n",indent_,dup_count);
  printf("%s parent_count: %d\n",indent_,parent_count);
  printf("%s parent_rid: %d\n",indent_,parent_rid);
  printf("%s parent_gen: %d\n",indent_,parent_gen);
  printf("%s parent_Class: %d\n",indent_,parent_Class);

  printf("%s(basic_record) double data:\n", indent_);
  printf("%s w: %e\n",indent_,w);
  printf("%s r2: %e\n",indent_,get_r2());
  if (print_all_) print_record();
}
void event_record::print_event_record(int nlead_, int npool_)
{
  printf("(event_record) record %d of %d leaders, %d pool.\n", rid,nlead_,npool_);
  printf("(event_record) int data:\n");
  printf(" nfobs: %d\n",nfobs);
  printf(" nfstable: %d\n",nfstable);
  printf(" nfunstable: %d\n",nfunstable);

  printf("(event_record) double data:\n");
  printf(" r2net: %e\n",r2net);
  printf(" r2stable: %e\n",r2stable);
  printf(" r2unstable: %e\n",r2unstable);

  printf("(event_record) evframe_bead: "); print_row_vec(evframe_bead,nbeads);
  printf("(event_record) r2stable_bead: "); print_row_vec(r2stable_bead,nbeads);
  printf("(event_record) r2unstable_bead: "); print_row_vec(r2unstable_bead,nbeads);
  printf("(event_record) alpha_bead: "); print_row_vec(alpha_bead,nbeads);

  print_basic_record();
}

// workers -----------------------------------------------------------------
void basic_thread_worker::print_basic_tw(const char indent_[], double sigma_scaled_)
  {print_event_block(sigma_scaled_,indent_);}
void MH_examiner::print_MH_examiner(int thread_count_, double sigma_scaled_)
{
  printf("(MH_examiner) worker %d of %d:\n", thread_id, thread_count_);
  print_basic_tw("     ",sigma_scaled_);
}

// solvers -----------------------------------------------------------------
void basic_MH_trainer::print_basic_MHT(const char indent_[])
  {print_MHT(indent_); print_gaussian_likelihood(indent_);}

void MH_genetic::print_MH_genetic()
{
  printf("(MH_genetic) Class_count: %d, event_block_count %d\n", Class_count,event_block_count);
  printf("(MH_genetic) prob_best: %e, prob_worst %e\n", prob_best,prob_worst);
  print_basic_MHT("     ");
  print_event_block(sigma_scaled,"     ");
}
void MH_genetic::verbose_find_events_1()
  {printf("(MH_genetic::find_events) Found event block %d with generation %d.", event_block_count,gen_count);}
void MH_genetic::verbose_find_events_2()
  {printf(" Earliest, latest frames: %d, %d\n", stev_earliest, stev_latest);}
void MH_genetic::verbose_find_events_3()
  {printf("(bead,frame), chronologically:\n"); for (int i = 0; i < nbeads; i++) printf("(%d %d)\n",comps_ordered[i],stev_ordered[i]);}
void MH_genetic::verbose_set_objective_1()
{
  printf("(MH_genetic::set_objective) %d %s frames, rho2 = %e. ",
  (stable_flag)?event_block::nst_stable():event_block::nst_unstable(),
  (stable_flag)?"STABLE":"UNSTABLE",
  rho2_objective);
}
void MH_genetic::verbose_set_objective_2()
  {printf("nreplace: %d, r2 best (%d): %e, r2 worst (%d): %e. ",nreplace,bleader_rid,br2,wleader_rid,wr2);}
void MH_genetic::verbose_train_objective_1(int nit_, int nsuccess_local_)
  {printf("(MH_genetic::train_objective) gen %d, Class %d, nit %d: %d (%d) candidates, best candidate (%d) r2 = %e. ", gen_count, Class_count, nit_, nsuccess, nsuccess_local_, bpool->rid, bpool->get_r2());}
void MH_genetic::verbose_train_objective_2()
  {printf("%d replacements, r2 best (%d): %e, r2 worst (%d): %e. ", nreplace, bleader_rid, br2, wleader_rid, wr2);}
void MH_genetic::verbose_respawn_pool(int offset_)
{
  if (offset_) printf("RESPAWNING: %d reloads, ", offset_);
  else printf("RESPAWNING: ");
  printf("%d redraws, %d duplicates (%d unique).\n", nredraw, ndup, ndup_unique);
}
