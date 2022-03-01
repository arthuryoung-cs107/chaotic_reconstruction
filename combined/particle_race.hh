#ifndef RACE_HH
#define RACE_HH

#include <cstdio>
#include <gsl/gsl_rng.h>
#include "swirl.hh"
#include "RYdat2AYdat.hh"
#include "AYlinalg.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

struct referee
{
    /** The maximum simulation timestep to use. */
    const double dt_sim;
    /** The length scale of the Gaussian (in pixels) for updating the particle
     * weights. */
    const double gau_scale;
    /** The coefficient that appears in the Gaussian for updating the particle
     * weights. */
    const double gau_coeff;


    /** number of particles being tested in each generation */
    const int npool;
    /** The number of best performing particles we store */
    const int nleaders;
    /** The length of the parameter vector corresponding to a particle */
    const int param_len;


    double ** pool = NULL;
    double ** leaders = NULL;

    double * l2_score = NULL;
    int * frame_score = NULL;

    referee(int npool_, int n_leaders_, int param_len_, double dt_sim_, double gau_scale_): npool(npool_), nleaders(n_leaders_), param_len(param_len_), dt_sim(dt_sim_), gau_scale(gau_scale_), gau_coeff(0.5/(gau_scale*gau_scale)) {}
    referee(referee &ref_): npool(ref_.npool), nleaders(ref_.nleaders), param_len(ref_.param_len), dt_sim(ref_.dt_sim), gau_scale(ref_.gau_scale), gau_coeff(0.5/(gau_scale*gau_scale)) {}
    ~referee()
    {
      if (pool!=NULL) free_AYdmatrix(pool);
      if (leaders!=NULL) free_AYdmatrix(leaders);
      if (l2_score!=NULL) delete l2_score;
      if (frame_score!=NULL) delete frame_score;
    }
    void alloc_records()
    {
      pool = AYdmatrix(npool, param_len);
      leaders = AYdmatrix(nleaders, param_len);
      l2_score = new double[nleaders];
      frame_score = new int[nleaders]; 
    }
};

class runner : public swirl
{
    public:
      int thread_id;
      runner(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int n_, int thread_id_);
      ~runner();

    private:

};

class race : public referee {
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
        /** The time points of the snapshots. */
        double* ts;
        /** The bead positions at the snapshots. */
        double* xs;
        /** The dish angle data at the snapshots. */
        double* d_ang;


        ODR_struct odr;


        race(referee &rparam,swirl_param &sp_min_,swirl_param &sp_max_,swirl_param &sp_rnd_,wall_list &wl_,double t_phys_,ODR_struct &odr_,int offset=0);
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
        /** The array of GSL random number generators. */
        gsl_rng** const rng;
        /** The array of proximity grids. */
        proximity_grid** const pg;

        /** Array of swirl simulations used by each thread to test particles */
        runner ** run;

#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

#endif
