#ifndef FIL_PARAM_HH
#define FIL_PARAM_HH

struct fil_param {
    /** The maximum simulation timestep to use. */
    const double dt_sim;
    /** The threshold fraction of effective particles at which a resampling
     * step is triggered. */
    const double rs_thresh;
    /** The size of the position perturbation (in pixels) to apply when
     * initializing the simulation from the experimental data. */
    const double x_pert;
    /** The size of the velocity perturbation to apply at each frame step. */
    const double v_pert;
    /** The size of the angular velocity perturbation to apply at each frame
     * step. */
    const double omega_pert;
    /** The training wheels parameter, setting how much drift there is of the
     * particles in the simulation toward true measurements. */
    const double t_wheels;
    /** The length scale of the Gaussian (in pixels) for updating the particle
     * weights. */
    const double gau_scale;
    /** The coefficient that appears in the Gaussian for updating the particle
     * weights. */
    const double gau_coeff;
    fil_param(double dt_sim_,double rs_thresh_,double x_pert_,double v_pert_,double omega_pert_,double t_wheels_,double gau_scale_)
        : dt_sim(dt_sim_), rs_thresh(rs_thresh_), x_pert(x_pert_),
        v_pert(v_pert_), omega_pert(omega_pert_), t_wheels(t_wheels_),
        gau_scale(gau_scale_), gau_coeff(0.5/(gau_scale*gau_scale)) {}
    fil_param(fil_param &f) : dt_sim(f.dt_sim), rs_thresh(f.rs_thresh),
        x_pert(f.x_pert), v_pert(f.v_pert), omega_pert(f.omega_pert),
        t_wheels(f.t_wheels), gau_scale(f.gau_scale), gau_coeff(f.gau_coeff) {}
};

#endif
