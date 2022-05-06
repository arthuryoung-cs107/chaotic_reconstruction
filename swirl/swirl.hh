#ifndef SWIRL_HH
#define SWIRL_HH

#include "swirl_param.hh"
#include "particle.hh"
#include "p_grid.hh"
#include "dat_store.hh"
#include "wall.hh"
#include "RYdat2AYdat.hh"

#include <vector>

class swirl : public swirl_param {
    public:
        /** The total number of particles. */
        const int n;
        /** A pointer to the particle information. */
        particle* const q;
        /** The probability weight of this swirling simulation. */
        double wei;
        /** The logarithm of the factor to apply to the weight, stored
         * separately to avoid numerical overflow. */
        double logfac;
        /** The time. */
        double time;
        /** The dish angle. */
        double ctheta;
        /** The dish angular velocity. */
        double comega;
        /** The dish center x position. */
        double cx;
        /** The dish center y position. */
        double cy;
        /** The dish center y position. */
        double cz;
        /** The dish x velocity. */
        double cvx;
        /** The dish y velocity. */
        double cvy;
        /** The dish z velocity. */
        double cvz;
        /** A pointer to thedata structure for detecting inter-particle
         * contacts. */
        proximity_grid *pg;
        /** A pointer to the class for storing the particle positions,
         * for creating synthetic data. */
        dat_store *dstore;
        ODR_struct *odr;
        swirl(swirl_param &sp,proximity_grid *pg_,wall_list &wl,int n_);
        swirl(swirl **sw,int npar,double gamma);
        ~swirl();
        void copy(swirl &sw);
        void import(const char* filename);
        void import1(const char* filename, int max_beads_, int index_);
        void import_true(const char* filename);
        void solve(double dur,double dt,int frames);
        void advance(double dur,double ctheta_,double comega_,double dt);
        void step_forward(double dt);
        void calc_forces();
        void output(int k);
        void output(FILE *fp);
        void setup_output_dir(const char *odir_);
        void setup_output_dir(ODR_struct *odr_, bool make_filter_inputs_flag_=true, bool solve_verbose_=true);
        void init_positions(double time_,double ctheta_,double *f,gsl_rng *r);
        void update_weight(double *f,double dur,double t_wheels,double *lnorm);
        void output_state(FILE *fp);
        void accumulators(double *tdig);
        void jiggle(swirl_param &sp_min,swirl_param &sp_max,swirl_param &sp_rnd,gsl_rng *r,double v_pert,double omega_pert);

    protected:
      // we need this to be visible from runner.cc
      inline void set_swirl(double ctheta_,double comega_)
      {
          ctheta=ctheta_;
          comega=comega_;
          set_dish_coords();
      }

    private:
        /** A list of walls that make up the swirling dish. */
        wall_list &wl;
        /** The output directory filename. */
        char *odir;
        /** A temporary buffer for assemble output filenames. */
        char *obuf;
        inline void update_swirl(double dt) {
            ctheta+=dt*comega;
            set_dish_coords();
        }
        void set_dish_coords();
        void contact(int id1,int id2,double delx,double dely,double delz,double rsq,double dt);

        bool ODR_struct_flag = false;
        bool make_filter_inputs_flag;
        bool solve_verbose;
};

#endif
