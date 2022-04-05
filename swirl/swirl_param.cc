#include "swirl_param.hh"

#include <gsl/gsl_randist.h>

/** Initializes the swirling parameters as a fraction of the difference between
 * two existing swirling parameter classes. This is used to create a vector of random
 * displacements for applying during the filtering.
 * \param[in] sp_min the minimum parameters.
 * \param[in] sp_max the maximum parameters.
 * \param[in] frac the fraction of the difference to use. */
swirl_param::swirl_param(swirl_param &sp_min,swirl_param &sp_max,double frac) :
    rad(0), mass(0), momi(0) {
    double *dmin=&sp_min.Kn,*dmax=&sp_max.Kn,*d=&(this->Kn);
    for(int k=0;k<sp_num_vparams;k++) d[k]=frac*(dmax[k]-dmin[k]);
}

/** Initializes the swirling parameters, where each is chosen uniformly over a range.
 * \param[in] sp_min the minimum parameters.
 * \param[in] sp_max the maximum parameters.
 * \param[in] r a pointer to the GSL random number generator to use. */
swirl_param::swirl_param(swirl_param &sp_min,swirl_param &sp_max,gsl_rng *r)
    : rad(sp_min.rad), mass(sp_min.mass), momi(sp_min.momi) {
    double *dmin=&sp_min.Kn,*dmax=&sp_max.Kn,*d=&(this->Kn);
    for(int k=0;k<sp_num_vparams;k++) d[k]=dmin[k]+(dmax[k]-dmin[k])*gsl_rng_uniform(r);
}

/** Initializes the swirling parameters to be zero, in preparation for
 * computing them as an average from other sets of swirling parameters.
 * \param[in] sp a reference to the swirling parameter to get the basic bead info from.
 * \param[in] l a dummy floating point number. */
swirl_param::swirl_param(swirl_param &sp,double l) : rad(sp.rad), mass(sp.mass),
    momi(sp.momi) {
    double *d=&(this->Kn);
    for(int k=0;k<sp_num_vparams;k++) d[k]=0.;
}

/** Adds contributions to the swirling parameters from another set of swirling
 * parameters.
 * \param[in] sp the parameters to add.
 * \param[in] gw the scaling factor to apply. */
void swirl_param::add_param(swirl_param &sp,double gw) {
    double *dp=&sp.Kn,*d=&(this->Kn);
    for(int k=0;k<sp_num_vparams;k++) d[k]+=gw*dp[k];
}

/** Copies the swirling parameters from another set. The function assumes that
 * bead mass and bead radius are the same.
 * \param[in] sp the parameters to copy. */
void swirl_param::copy_param(swirl_param &sp) {
    double *dp=&sp.Kn,*d=&(this->Kn);
    for(int k=0;k<sp_num_vparams;k++) d[k]=dp[k];
}

/** Scales the swirling parameters.
 * \param[in] cw the scaling factor to apply. */
void swirl_param::scale_param(double cw) {
    double *d=&(this->Kn);
    for(int k=0;k<sp_num_vparams;k++) d[k]*=cw;
}

/** Applies random displacements, ensuring that they do not go outside minimum
 * and maximum limits.
 * \param[in] sp_min the minimum limits on the parameters.
 * \param[in] sp_max the maximum limits on the parameters.
 * \param[in] sp_rnd the random displacements to apply to the parameters.
 * \param[in] r a pointer to the GSL random number generator to use. */
void swirl_param::jiggle_param(swirl_param &sp_min,swirl_param &sp_max,swirl_param &sp_rnd,gsl_rng *r) {
    double *dr=&sp_rnd.Kn,*dmin=&sp_min.Kn,*dmax=&sp_max.Kn,*d=&(this->Kn);
    for(int k=0;k<sp_num_vparams;k++) {

        // Apply a random displacement to the parameter
        d[k]+=dr[k]*(2*gsl_rng_uniform(r)-1);

        // If the new parameter exceeds the limits, then move it back into the
        // range
        if(d[k]<dmin[k]) d[k]=dmin[k];
        else if(d[k]>dmax[k]) d[k]=dmax[k];
    }
}
