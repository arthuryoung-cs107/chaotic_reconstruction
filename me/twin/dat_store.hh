#ifndef DAT_STORE_HH
#define DAT_STORE_HH

#include "particle.hh"

#include <vector>

class dat_store {
    public:
        const double t_phys;
        /** The dish x center in the images. */
        const double cx_im;
        /** The dish y center in the images. */
        const double cy_im;
        /** The bead diameter in the images. */
        const double cl_im;
        dat_store(double t_phys_,double cx_im_,double cy_im_,double cl_im_)
            : t_phys(t_phys_), cx_im(cx_im_), cy_im(cy_im_),
            cl_im(cl_im_), nsnap(0) {}
        void snapshot(double time,particle *q,int n_,double cx,double cy,double ctheta);
        void write(const char* filename);
        void print_theta_info();
    private:
        /** The number of particles. */
        int n;
        /** The number of snapshots. */
        int nsnap;
        /** The time points of the snapshots. */        
        std::vector<float> ts;
        /** The bead positions at the snapshots. */        
        std::vector<float> xs;
        /** The array of pointers to swirling simulations. */
        std::vector<float> d_ang;
};


#endif
