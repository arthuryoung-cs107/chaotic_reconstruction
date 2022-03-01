#ifndef RUNNER_HH
#define RUNNER_HH

#include "swirl.hh"

class runner : public swirl
{
    public:
      int thread_id;
      runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_);
      ~runner();

    private:

};

#endif
