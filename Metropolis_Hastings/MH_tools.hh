#ifndef MH_TOOLS_HH
#define MH_TOOLS_HH

struct gaussian_likelihood
{
  gaussian_likelihood(double noise_tol_, double cl_): noise_tol(noise_tol_), cl(cl_), coeff(cl*noise_tol) {}
  ~gaussian_likelihood();

  const double  noise_tol,
                cl,
                coeff;

  inline double expected_r2(int * obs_, int nbeads_)
  {
    int net_obs=0;
    for (int i = 0; i < nbeads_; i++) net_obs+=obs_[i];
    return coeff*coeff*((double)net_obs);
  }
};

struct event_detector
{
  event_detector(double alpha_tol_);
  ~event_detector();

  const double alpha_tol;

  inline double alpha_comp(double *a_, double t_p1, double t_m1)
  {return log(a_[0]/a_[2])/log(t_p1/t_m1);}
  inline double alpha_comp(double a_p1, double a_m1, double t_p1, double t_m1)
  {return log(a_p1/a_m1)/log(t_p1/t_m1);}
};

#endif
