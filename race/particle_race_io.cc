#include "particle_race.hh"

void record::print_record()
{
  printf("record ID:%d, fr:%d, l2:%e, params: ", global_index, frscore, l2score);
  for (int i = 0; i < 12; i++) printf("%e ", params[i]);
}




// RUNNER IO

void runner::print_raw_ics()
{
  printf("runner #%d: (t_phys, t0_raw, t0, ctheta0) = (%e %e %e %e) \n", thread_id, t_phys, t0_raw, t0, ctheta0);
  for (int j = 0; j < n; j++)
  {
    printf("  (bead %d): %e %e \n", j, x0[2*j], x0[2*j+1]);
  }
}

void runner::print_current_pos()
{
  printf("runner #%d: (t, ctheta) = (%e %e)\n", thread_id, time, ctheta);
  printf("(cx_im, cy_im, cl_im, cx, cy) = (%e %e %e %e %e)\n", cx_im, cy_im, cl_im, cx, cy);
  for (int j = 0; j < n; j++)
    printf("  (bead %d): %e %e %e\n", j, q[j].x, q[j].y, q[j].z);
}

void runner::print_params()
{
  printf("runner #%d: (%e %e) ", thread_id, rad, mass);
  double *par_vals = &Kn;
  for (int j = 0; j < param_len; j++) printf("%e ", par_vals[j]);
  printf("\n");
}





// RUNNER IO

void referee::print_referee_params()
{
  const double * double_params = &(dt_sim);
  const int * int_params = &(npool);

  printf("double params:\n");
  for (int i = 0; i < 7; i++) printf("%e ", double_params[i]);
  printf("\n\nint params:\n");
  for (int i = 0; i < 3; i++) printf("%d ", int_params[i]);
  printf("\n");
}




// RACE IO

void race::stage_diagnostics()
{
  odr->nlead = nlead;
  odr->npool = npool;
  odr->len = param_len;
  odr->dup_vec = dup_vec;
  odr->sample_weights = sample_weights;
  odr->leaders = leaders;
  odr->staged_flag = true;
}
void race::print_pool_leaders()
{
  for (int i = 0; i < pool_success_count; i++)
    printf("gen %d, pool particle %d: (ID, frame score, l2score) = (%d %d %e)\n", gen_count, i, pool_leaders[i]->global_index, pool_leaders[i]->frscore, pool_leaders[i]->l2score);
}
void race::print_leaders()
{
  printf("\n");
  for (int i = 0; i < leader_count; i++)
  {
    printf("(gen, i, ID, frame score, l2score) = (%d %d %d |%d| %e)\n", gen_count, i, leaders[i]->global_index, leaders[i]->frscore, leaders[i]->l2score);
  }
  print_best();
  print_worst_leader();
  printf("\n");
}
void race::print_leader_candidates()
{
  printf("\n");
  for (int i = 0; i < nlead; i++)
  {
    printf("gen %d, candidate particle %d: (ID, frame score, l2score) = (%d %d %e)\n", gen_count, i, leader_board[i]->global_index, leader_board[i]->frscore, leader_board[i]->l2score);
  }
  printf("\n");
}
void race::print_best()
{
 int best = find_best(leaders, leader_count);
 printf("BEST\n(gen, i, ID, frame score, l2score) = (%d %d %d |%d| %e)\n", gen_count, best, leaders[best]->global_index, leaders[best]->frscore, leaders[best]->l2score);
}
void race::print_worst_leader()
{
 int worst = find_worst(leaders, leader_count);
 printf("WORS\n(gen, i, ID, frame score, l2score) = (%d %d %d |%d| %e)\n", gen_count, worst, leaders[worst]->global_index, leaders[worst]->frscore, leaders[worst]->l2score);
}
void race::print_reference_positions(int len_)
{
  int k=0;
  printf("(beads, Frames) = (%d %d)\n", n, Frames);
  for (int i = 0; i < len_; i++)
  {
    printf("(Frame #%d): (t, theta) = (%e %e)\n", i, ts[i], d_ang[i]);
    for (int j = 0; j < n; j++)
    {
      printf("  (bead %d): %e %e \n", j, xs[k], xs[k+1]);
      k+=2;
    }
    printf("\n");
  }
}
void race::print_pool_params()
{
  for (int i = 0; i < npool; i++)
  {
    printf("pool params %d: ", i);
    for (int j = 0; j < param_len; j++) printf("%e ", pool_params[i][j]);
    printf("\n");
  }
}
void race::print_pool_records()
{
  for (int i = 0; i < npool; i++)
  {
    printf("(pool #%d) ", i); pool[i]->print_record();
    printf("\n");
  }
}
void race::make_best_swirl(char * name_)
{
  if (leader_count>0)
  {
    int best = find_best(leaders, leader_count);
    double * sbest = leaders[best]->params;

    runner * runner0 = runners[0];
    runner0->reset_sim(sbest);

    // Solve the system
    ODR_struct * odr_swbest = odr->spawn_swirlODR(name_);
    odr_swbest->set_vidspecs(t_phys, runner0->cx_im, runner0->cy_im, runner0->cl_im);
    runner0->setup_output_dir(odr_swbest, false, false);
    odr_swbest->write_sparam(runner0, ".sparam_best");

    runner0->solve(120,0.0005,1200);

    odr_swbest->end_writing();

    delete odr_swbest;
  }
  else printf("make_best_swirl: no leader candidates. Swirl not produced.\n");
}
