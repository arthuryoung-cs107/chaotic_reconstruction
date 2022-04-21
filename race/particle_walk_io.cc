#include "particle_walk.hh"

void walk::make_best_swirl(char * name_)
{
  if (leader_count>0)
  {
    int best = find_best_grade(leaders, leader_count);
    double * sbest = leaders[best]->params;

    walker * walker0 = walkers[0];
    walker0->reset_sim(sbest);

    // Solve the system
    ODR_struct * odr_swbest = ped->spawn_swirlODR(name_);
    odr_swbest->set_vidspecs(t_phys, walker0->cx_im, walker0->cy_im, walker0->cl_im);
    walker0->setup_output_dir(odr_swbest, false, false);
    odr_swbest->write_sparam(walker0, ".sparam_best");

    walker0->solve(120,0.0005,1200);

    odr_swbest->end_writing();

    delete odr_swbest;
  }
  else printf("make_best_swirl: no leader candidates. Swirl not produced.\n");
}
