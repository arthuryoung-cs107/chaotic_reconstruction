#include "particle_race.hh"

void race::print_pool_leaders()
{
  for (int i = 0; i < pool_success_count; i++)
  {
    printf("gen %d, pool particle %d: (ID, frame score, l2score) = (%d %d %e)\n", gen_count, i, pool_leaders[i]->global_index, pool_leaders[i]->frscore, pool_leaders[i]->l2score);
  }
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

 void race::make_best_swirl()
 {
   if (leader_count>0)
   {
     int best = find_best(leaders, leader_count);
     double * sbest = leaders[best]->params;

     double g_phys=9.804,                 // Gravity (m/s^2)
            d_phys=0.00635,               // Diameter (m)
            t_phys=sqrt(d_phys/g_phys);   // Time unit (s)

            swirl_param sparam(sbest);

     // Create the hexagonal dish
     wall_list wl;
     const double r=5.72,fa=sqrt(0.75);
     wall_floor wf(0.);
     wall_par_planes wp0(0,1,0,r),wp1(fa,0.5,0,r),wp2(fa,-0.5,0,r);
     wl.add_wall(&wf);
     wl.add_wall(&wp0);
     wl.add_wall(&wp1);
     wl.add_wall(&wp2);

     // Set the initial positions of the splines
     proximity_grid pg;
     swirl sw(sparam,&pg,wl,3);
     sw.import("input_dir/input3.dat");

     // Solve the system
     ODR_struct odr("./dat_dir/", "race_3beads_swirlbest.odr/", "pts");
     odr.set_vidspecs(t_phys, sbest[11], sbest[12], sbest[13]);
     sw.setup_output_dir(&odr);

     sw.solve(120,0.0005,1200);

     odr.end_writing();
     odr.print_time_rotation();

     AYvec sparam_out(12);
     for (int i = 0; i < 12; i++) sparam_out.A_ptr[i] = sbest[i];
     sparam_out.fprintf_vec("./dat_dir/race_3beads_swirlbest.odr/sbest");

   }
   else
   {
     printf("make_best_swirl: no leader candidates. Swirl not produced.\n");
   }
 }
