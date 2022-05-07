#include "particle_relay.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

extern "C"
{
  #include "AYaux.h"
}

void relay::find_events(int min_frame_, int latest_frame_, bool verify_)
{
  bool first2finish = true;
  int ncheck;
  record ** rec_check;
  if (verify_)
  {
    ncheck=nlead;
    rec_check = leaders;
  }
  else
  {
    ncheck=npool;
    rec_check = pool;
  }
  // determine the next sequence of events in the data, requiring that that they must occur after the min_frame_ parameter
  #pragma omp parallel
    {
      runner *rt = runners[thread_num()];
      rt->clear_event_data();
  #pragma omp for nowait
      for (int i = 0; i < ncheck; i++)
      {
        rt->detect_events(rec_check[i], 0, min_frame_, latest_frame_);
      }
  #pragma omp critical
      {
        if (first2finish)
        {
          for (int i = 0; i < n*Frames; i++) global_event_frame_count[0][i] = rt->event_frame_count[0][i];
          first2finish=false;
        }
        else for (int i = 0; i < n*Frames; i++) global_event_frame_count[0][i] += rt->event_frame_count[0][i];
      }
    }
}


int relay::train_event_block(int event_block, int gen_max_, double tol_leeway_, bool train_full_)
{
  int end_point = (train_full_)?Frames:latest_event;
  bool training=true;
  do
  {
    int success_local=0;
#pragma omp parallel
    {
      runner *rt = runners[thread_num()];
      int * event_frames_ordered = ev->event_sorted;
      int * bead_order = ev->bead_order;
#pragma omp for reduction(+:success_local)
      for (int i = 0; i < npool; i++)
      {
        success_local += rt->run_relay(pool[i], 0, earliest_event, end_point, event_frames_ordered, bead_order, residual_worst);
      }
    }
    gen_count++;
    printf("(gen %d): %d candidates. ", gen_count, success_local);
    if (check_pool_results(tol_leeway_)) training=false; // we win
    else if (gen_count == gen_max_) training=false; // we give up
    else resample_pool(); // we try again
    if (debugging_flag) rep->write_gen_diagnostics(gen_count, leader_count);
  } while (training);

  return latest_event;
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

bool relay::check_pool_results(double tol_leeway_)
{
  int pool_candidates_local = pool_candidates = collect_candidates();
  // if we now have a full leader roster
  repl_count=0;
  if (leader_count + pool_candidates_local >= nlead)
  {
    // fill up remainder of leaders
    for (int i = 0, gap = nlead-leader_count; i < gap; i++)
    {
      repl_count++;
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
    residual_worst = leader_board[worst_best]->residual; root_res_worst = leader_board[worst_best]->root_residual;

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
  residual_best = bl_rec->residual; root_res_best = bl_rec->root_residual;
  residual_worst = wl_rec->residual; root_res_worst = wl_rec->root_residual;
  printf("Best/Worst: (ID, gen, parents, residual) = (%d/%d %d/%d %d/%d %e/%e), tau_sqr=%e, %d replacements. ", bl_rec->global_index, wl_rec->global_index, bl_rec->gen, wl_rec->gen, bl_rec->parent_count, wl_rec->parent_count, residual_best, residual_worst, tau_sqr, repl_count);
  if (root_res_worst<(1.0+tol_leeway_)*tau) return true;
  else if (repl_count==0) return true;
  return false;
}
