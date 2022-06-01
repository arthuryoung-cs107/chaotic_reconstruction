#include "MH_learning.hh"

extern "C"
{
  #include "AYaux.h"
}

MH_io::MH_io(char * proc_loc_, char * test_dir_, char * data_name_, int id_, bool noise_data_, double noise_sigma_): proc_loc(proc_loc_), test_dir(test_dir_), data_name(data_name_),
id(id_), fullbuf_len(strlen(proc_loc)+strlen(test_dir)+strlen(data_name_)),
noise_data(noise_data_), noise_sigma(noise_sigma_),
obuf(new char[fullbuf_len+100]), ibuf(new char[fullbuf_len+10])
{
  sprintf(obuf,"%s%s",proc_loc,test_dir); obuf_len=strlen(obuf);
  sprintf(ibuf,"%s%s%s",proc_loc,test_dir,data_name); ibuf_len=strlen(ibuf);
  read_fisml(ibuf+ibuf_len);
}

MH_io::~MH_io()
{delete [] obuf; delete [] ibuf;}

void MH_io::read_fisml(char *ibuf_)
{
  // there is certainly a better way of doing this
  int int_buf[3];
  strcpy(ibuf_, ".fisml"); FILE * ref_sml = fopen(ibuf, "r");
  fscanf(ref_sml,"%d %d %d\n",int_buf,int_buf+1,int_buf+2);
  fscanf(ref_sml,"%d %d %d\n",int_buf,int_buf+1,int_buf+2); Frames=int_buf[2];
  fscanf(ref_sml,"%d %d %d\n",int_buf,int_buf+1,int_buf+2); nbeads=int_buf[1];
  fclose(ref_sml);
}

void MH_io::load_reference(double *ts_, double *xs_, double *d_ang_, double * comega_s_, double t_phys_)
{
  sprintf(ibuf+ibuf_len, ".filin");
  FILE * inputs = fopen(ibuf, "r");
  fread_safe(ts_, sizeof(double), Frames, inputs);
  fread_safe(xs_, sizeof(double), 2*Frames*nbeads, inputs);
  fread_safe(d_ang_, sizeof(double), Frames, inputs);
  fclose(inputs);

  if (noise_data)
  {
    AYrng ran;
    ran.rng_init_gsl(1);
    for (int i = 2*nbeads; i < 2*Frames*nbeads; i++) xs_[i]+=ran.rand_gau_gsl(0.0,noise_sigma);
  }

  // initialize the dish velocity data. Will save time on calculations later
  comega_s_[0]=0.0;
  for (int i = 1; i < Frames; i++)
  {
    double comega=d_ang_[i]-d_ang_[i-1];
    if(comega>M_PI) comega-=2*M_PI; else if(comega<-M_PI) comega+=2*M_PI;
    comega_s_[i] = t_phys_*comega/(ts_[i]-ts_[i-1]);
  }

}

thread_worker::thread_worker(swirl_param &sp_, proximity_grid * pg_, wall_list &wl_, int thread_id_, int ichunk_len_, int dchunk_len_, int param_len_, int nbeads_, int Frames_, int nlead_, int npool_, double dt_sim_, double t_phys_, double *ts_, double *xs_, double *d_ang_, double *comega_s_): swirl(sp_, pg_, wl_, nbeads_),
thread_id(thread_id_), ichunk_len(ichunk_len_), dchunk_len(dchunk_len_),
param_len(param_len_), nbeads(nbeads_), Frames(Frames_), nlead(nlead_), npool(npool_),
dt_sim(dt_sim_), t_phys(t_phys_),
ts(ts_), xs(xs_), d_ang(d_ang_), comega_s(comega_s_),
ichunk(new int[ichunk_len]), dchunk(new double[dchunk_len])
{pvals = &Kn;}

thread_worker::~thread_worker()
{delete [] ichunk; delete [] dchunk;}

MH_trainer::MH_trainer(MH_params &par_, swirl_param &sp_min_, swirl_param &sp_max_, wall_list &wl_, double t_wheels_): MH_params(par_),
sp_min(sp_min_), sp_max(sp_max_), wl(wl_),
t_wheels(t_wheels_),
ts(new double[par_.Frames]), xs(new double[2*par_.nbeads*par_.Frames]), d_ang(new double[par_.Frames]), comega_s(new double[par_.Frames]),
#ifdef _OPENMP
nt(omp_get_max_threads()), // each thread is a runner
#else
nt(1), // only one runner
#endif
pg(new proximity_grid*[nt]), rng(new AYrng*[nt]),
u_chunk(AYdmatrix(par_.nlead+par_.npool, par_.param_len))
{
  io->load_reference(ts, xs, d_ang, comega_s, t_phys);

  // Set up each runner's personal data
#pragma omp parallel
  {
    int t=thread_num();
    pg[t]=new proximity_grid();
    rng[t]=new AYrng();
    rng[t]->rng_init_gsl(t+1);
  }
}

MH_trainer::~MH_trainer()
{
  delete [] ts;
  delete [] xs;
  delete [] d_ang;
  delete [] comega_s;
  free_AYdmatrix(u_chunk);

  for (int i = 0; i < nt; i++)
  {
    delete pg[i];
    delete rng[i];
  }
  delete [] pg;
  delete [] rng;
}
