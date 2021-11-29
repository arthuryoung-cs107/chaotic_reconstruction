#ifndef FILTER_HH
#define FILTER_HH

#include <cstdio>
#include <gsl/gsl_rng.h>

#include "fil_param.hh"
#include "swirl_param.hh"
#include "swirl.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

/** The number of failed frames to allow before quitting. */
const int filter_nfail_thresh=4;

class filter : public fil_param {
    public:
        /** The minimum parameters. */
        swirl_param sp_min;
        /** The maximum parameters. */
        swirl_param sp_max;
        /** The random step parameters. */
        swirl_param sp_rnd;
        /** The simulation time in seconds. */
        double t_phys;
        /** The minimum L^2 measure. */
        double min_l2;
        /** The minimum L^infty measure. */
        double min_linf;
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
        float* ts;
        /** The bead positions at the snapshots. */
        float* xs;
        /** The dish angle data at the snapshots. */
        float* d_ang;
        /** The array of pointers to swirling simulations. */
        swirl** sw;
        filter(fil_param &fparam,swirl_param &sp_min_,swirl_param &sp_max_,swirl_param &sp_rnd_,wall_list &wl_,double t_phys_,const char* filename,int offset=0);
        ~filter();
        void init(int npar_);
        void run(int frames);
        void resample();
        void step_frame();
        void set_snap(int k,short int xc,short int yc,double s);
        void setup_output_info(unsigned int fflags_,const char *odir_,int state_freq_=1);
        void output_particle(swirl *swp,const char* prefix);
        void output_data(swirl *swp,const char* prefix);
        void write_digest();
        void write_files();
        void output_states();
        int most_likely();
    private:
        /** The resampling buffer. */
        int* rsbuf;
        /** The types of file output to store. */
        unsigned int fflags;
        /** The frequency at which to store state information. */
        int state_freq;
        /** The number of threads. */
        const int nt;
        /** The thread particle allocation table. */
        int* const ttab;
        /** The thread resampling relocation table. */
        int* const rloc;
        /** The array of GSL random number generators. */
        gsl_rng** const rng;
        /** A reference to the list of walls for the swirling simulation. */
        wall_list &wl;
        /** The array of proximity grids. */
        proximity_grid** const pg;
        /** The output directory filename. */
        char *odir;
        /** A temporary buffer for assemble output filenames. */
        char *obuf;
        /** The file handle for the digest file, if used. */
        FILE *fdigest;
#ifdef _OPENMP
        inline int thread_num() {return omp_get_thread_num();}
#else
        inline int thread_num() {return 0;}
#endif
};

#endif
