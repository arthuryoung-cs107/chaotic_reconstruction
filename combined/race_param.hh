#ifndef RACE_PARAM_HH
#define RACE_PARAM_HH

struct race_param
{
    /** The maximum simulation timestep to use. */
    const double dt_sim;
    /** The length scale of the Gaussian (in pixels) for updating the particle
     * weights. */
    const double gau_scale;
    /** The coefficient that appears in the Gaussian for updating the particle
     * weights. */
    const double gau_coeff;
    race_param(double dt_sim_, double gau_scale_): dt_sim(dt_sim_), gau_scale(gau_scale_), gau_coeff(0.5/(gau_scale*gau_scale)) {}
};

#endif
