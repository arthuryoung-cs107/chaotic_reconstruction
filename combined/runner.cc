#include "particle_race.hh"

extern "C"
{
  #include "AYaux.h"
}

runner::runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_): swirl(sp_, pg_, wl_, n_), thread_id(thread_id_)
{}
runner::~runner()
{}
