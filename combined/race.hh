#ifndef RACE_HH
#define RACE_HH

#include <cstdio>
#include <gsl/gsl_rng.h>

#include "race_param.hh"
#include "swirl_param.hh"
#include "swirl.hh"
#include "RYdat2AYdat.hh"
#include "AYlinalg.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

/** The number of failed frames to allow before quitting. */
const int filter_nfail_thresh=4;

class race : public race_param {
    public:
        /** The minimum parameters. */
        swirl_param sp_min;
        /** The maximum parameters. */
        swirl_param sp_max;
        /** The random step parameters. */
        swirl_param sp_rnd;
        /** The simulation time in seconds. */
        double t_phys;
        /** The total number of beads. */
        int n;
        /** The total number of snapshots. */
        int nsnap;
        /** The total number of particles. */
        int npar;
        /** The failure count. */
        int nfail;
        /** The current frame. */
        int frame;
        /** The time points of the snapshots. */
        double* ts;
        /** The bead positions at the snapshots. */
        double* xs;
        /** The dish angle data at the snapshots. */
        double* d_ang;
        ODR_struct odr;
        race(race_param &rparam,swirl_param &sp_min_,swirl_param &sp_max_,swirl_param &sp_rnd_,wall_list &wl_,double t_phys_,ODR_struct &odr_,int offset=0);
        ~race();
        void init(int npar_);
        void run(int frames);
    private:
        /** The number of threads. */
        const int nt;
        /** The array AYuniform generators. */
        AYuniform * uni;
        /** A reference to the list of walls for the swirling simulation. */
        wall_list &wl;
        /** The array of proximity grids. */
        proximity_grid** const pg;
#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

#endif
