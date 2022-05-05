#include "particle_relay.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

extern "C"
{
  #include "AYaux.h"
}

void relay::train_leg(int smooth_gen_max_, int stiff_gen_max_)
{
  bool smooth_training=true;
  int half_gen_max = gen_max_/2;

  do // train off of the smooth data
  {
    int success_local=0;
#pragma omp parallel
    {
      runner *rt = runners[thread_num()];
#pragma omp for reduction(+:success_local) nowait
      for (int i = 0; i < npool; i++)
      {
        success_local += rt->run_relay(pool[i], 0, event_frames, earliest_event, latest_event, residual_worst);
      }
    }
    gen_count++;
    printf("(gen %d): %d candidates. ", gen_count, success_local);
    if (check_pool_results()) smooth_training=false; // we win
    else if (gen_count == half_gen_max) smooth_training=false; // we give up
    else resample_pool(); // we try again
    if (debugging_flag) rep->write_gen_diagnostics(gen_count, leader_count);
  } while (smooth_training);

  printf("\nSmooth training terminated. Beginning post-event training.");

  do // train off of full data
  {
    int success_local=0;
#pragma omp parallel
    {
      runner *rt = runners[thread_num()];
#pragma omp for reduction(+:success_local) nowait
      for (int i = 0; i < npool; i++)
      {
        success_local += rt->run_relay(pool[i], 0, event_frames, earliest_event, latest_event, residual_worst);
      }
    }
    gen_count++;
    printf("(gen %d): %d candidates. ", gen_count, success_local);
    if (check_pool_results()) training_underway=false; // we win
    else if (gen_count == half_gen_max) training_underway=false; // we give up
    else resample_pool(); // we try again
    if (debugging_flag) rep->write_gen_diagnostics(gen_count, leader_count);
  } while (training_underway);

}


int relay::collect_candidates()
{
  pool_success_count = 0;
  for (int i = 0; i < npool; i++)
    if (pool[i]->success)
      candidates[pool_success_count++] = pool[i];
  if (pool_success_count > nlead) // if we have lots of good candidates
  {
    // we seek to pick the nlead best grades for comparison.
    int worst_best = find_worst_record(candidates, nlead);
    for (int i = nlead; i < pool_success_count; i++)
      if (candidates[worst_best]->isworse(candidates[i]))
      {
        candidates[worst_best] = candidates[i];
        worst_best = find_worst_record(candidates, nlead);
      }
    return nlead;
  }
  else return pool_success_count;
}

bool relay::check_pool_results()
{
  int pool_candidates_local = pool_candidates = collect_candidates();
  // if we now have a full leader roster
  repl_count=0;
  if (leader_count + pool_candidates_local >= nlead)
  {
    // fill up remainder of leaders
    for (int i = 0, gap = nlead-leader_count; i < gap; i++)
    {
      leader_board[leader_count] = leaders[leader_count];
      leaders[leader_count++]->take_vals(candidates[--pool_candidates_local]);
    }

    int worst_best = find_worst_record(leader_board, nlead);

    // consider the pool candidates, which are positioned adjacent to the current leaders on the leaderboard
    for (int i = nlead; i < nlead + pool_candidates_local; i++)
      if (leader_board[worst_best]->isworse(leader_board[i]))
      {
        leader_board[worst_best] = leader_board[i];
        worst_best = find_worst_record(leader_board, nlead);
      }

    worst_leader = worst_best; // this will be the worst leader
    residual_worst = leader_board[worst_best]->residual;

    double resacc_local = 0.0;
    for (int i = 0; i < nlead; i++)
    {
      resacc_local+=leader_board[i]->residual;
      if (leaders[i]->global_index != leader_board[i]->global_index)
      {
        repl_count++;
        leaders[i]->take_vals(leader_board[i]);
        leader_board[i] = leaders[i];
      }
    }
    pos_res_global=resacc_local/((double)(leader_count-1));
  }
  // otherwise, we can just fill in the leaderboard
  else for (int i = 0; i < pool_candidates; i++)
  {
    repl_count++;
    leader_board[leader_count] = leaders[leader_count];
    leaders[leader_count++]->take_vals(candidates[i]);
  }

  best_leader = find_best_record(leaders, leader_count);
  record * wl_rec = leaders[worst_leader], * bl_rec = leaders[best_leader];
  residual_best = bl_rec->residual;
  printf("Best/Worst: (ID, gen, parents, residual) = (%d/%d %d/%d %d/%d %e/%e), %d replacements. ", bl_rec->global_index, wl_rec->global_index, bl_rec->gen, wl_rec->gen, bl_rec->parent_count, wl_rec->parent_count, residual_best, residual_worst, repl_count);
  if (sqrt(residual_worst)<tau) return true;
  return false;
}
