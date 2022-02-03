#ifndef SWIRL_PARAM_HH
#define SWIRL_PARAM_HH

#include <cstdio>
#include <cmath>

#include <gsl/gsl_rng.h>

/** The total number of varying parameters in swirling parameter class. */
const int sp_num_vparams=12;

struct swirl_param {
    /** The bead radius. */
    double rad;
    /** The bead mass. */
    double mass;
    /** The bead moment of inertia. */
    double momi;
    /** The bead normal spring constant. */
    double Kn; // likely to have more uncertainty relative to gamma
    /** The bead-bead normal damping. */
    double gamman_bead;
    /** The bead-base normal damping. */
    double gamman_base;
    /** The bead-wall normal damping. */
    double gamman_wall;
    /** The bead-bead friction coefficient. */
    double mu_bead;
    /** The bead-wall friction coefficient. */
    double mu_base;
    /** The bead-wall friction coefficient. */
    double mu_wall;
    /** The dish spinning amplitude. */
    double dsamp; // hold constant for all particles
    /** The dish x center in the images. */ // hold constant for all particles
    double cx_im; // hold constant for all particles
    /** The dish y center in the images. */ // hold constant for all particles
    double cy_im; // hold constant for all particles
    /** The bead diameter in the images. */ // hold constant for all particles
    double cl_im; // hold constant for all particles
    /** The wall scale factor. */
    double wall_sca;
    swirl_param(double rad_,double mass_,double Kn_,double gamman_bead_,double gamman_base_,double gamman_wall_,double mu_bead_,double mu_base_,double mu_wall_,double dsamp_,double cx_im_,double cy_im_,double cl_im_,double wall_sca_)
        : rad(rad_), mass(mass_), momi(0.4*mass*rad*rad), Kn(Kn_),
        gamman_bead(gamman_bead_), gamman_base(gamman_base_),
        gamman_wall(gamman_wall_), mu_bead(mu_bead_), mu_base(mu_base_),
        mu_wall(mu_wall_), dsamp(dsamp_), cx_im(cx_im_), cy_im(cy_im_),
        cl_im(cl_im_), wall_sca(wall_sca_) {}
    swirl_param(swirl_param &sp)
        : rad(sp.rad), mass(sp.mass), momi(sp.momi), Kn(sp.Kn),
        gamman_bead(sp.gamman_bead), gamman_base(sp.gamman_base),
        gamman_wall(sp.gamman_wall), mu_bead(sp.mu_bead), mu_base(sp.mu_base),
        mu_wall(sp.mu_wall), dsamp(sp.dsamp), cx_im(sp.cx_im), cy_im(sp.cy_im),
        cl_im(sp.cl_im), wall_sca(sp.wall_sca) {}
    swirl_param(swirl_param &sp,double l);
    swirl_param(swirl_param &sp_min,swirl_param &sp_max,double frac);
    swirl_param(swirl_param &sp_min,swirl_param &sp_max,gsl_rng *r);
    void add_param(swirl_param &sp,double gw);
    void scale_param(double cw);
    void copy_param(swirl_param &sp);
    void jiggle_param(swirl_param &sp_min,swirl_param &sp_max,swirl_param &swr,gsl_rng *r);
};

#endif
