#include "particle_race.hh"

extern "C"
{
  #include "AYaux.h"
}

referee::~referee()
{
  if (alloc_flag)
  {
    delete [] leader_board;
    for (int i = 0; i < 2*nlead; i++) delete leaders[i];
    delete [] leaders;
    free_AYdmatrix(lead_params);
  }
}
void referee::alloc_records()
{
  leaders = new record*[nlead+npool];
  leader_board = new record*[nlead+npool];
  for (int i = 0; i < nlead+npool; i++) leaders[i] = new record(i);
  pool = leaders + nlead;
  pool_leaders = leader_board + nlead;

  lead_params = AYdmatrix(nlead+npool, param_len);
  pool_params = lead_params + nlead;

  alloc_flag = true;
}

int find_worst(record ** r, int ncap)
{
  int worst_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[worst_index]->isbetter(r[i]))
      worst_index = i;
  return worst_index;
}

int find_best(record ** r, int ncap)
{
  int best_index = 0;
  for (int i = 1; i < ncap; i++)
    if (r[best_index]->isworse(r[i]))
      best_index = i;
  return best_index;
}
