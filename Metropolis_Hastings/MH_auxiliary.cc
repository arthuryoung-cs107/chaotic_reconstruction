#include "MH_auxiliary.hh"
#include "MH_tools.hh"

int set_special_u(int id_, double *vec_)
{
  int id_out=0;
  if ((id_>=0)&&(id_<special_u_count)) id_out=id_;
  else printf("WARNING: invalid parameter ID (%d). Using default parameters\n", id_);

  for (int i = 0; i < full_ulen; i++) vec_[i] = special_u[id_out][i];
  return id_out;
}

int set_special_u(const char *id_, double *vec_)
{
  if (strcmp(id_, "true")==0) return set_special_u(0, vec_);
  else if (strcmp(id_, "min")==0) return set_special_u(1, vec_);
  else if (strcmp(id_, "max")==0) return set_special_u(2, vec_);
  else if (strcmp(id_, "pert")==0) return set_special_u(3, vec_);
  else if (strcmp(id_, "b3i0r5")==0) return set_special_u(4, vec_);
  else return set_special_u(-1, vec_);
}

// MH_io

MH_io::MH_io(char * proc_loc_, char * test_dir_, char * data_name_, int id_, bool noise_data_, double noise_sigma_): fullbuf_len(strlen(proc_loc_)+strlen(test_dir_)+strlen(data_name_)),
proc_loc(proc_loc_), test_dir(test_dir_), data_name(data_name_),
id(id_),
noise_data(noise_data_), noise_sigma(noise_sigma_),
obuf(new char[fullbuf_len+100]), ibuf(new char[fullbuf_len+100])
{
  sprintf(obuf,"%s%s",proc_loc,test_dir); obuf_len=strlen(obuf);
  sprintf(ibuf,"%s%s%s",proc_loc,test_dir,data_name); ibuf_len=strlen(ibuf);
  read_fisml(ibuf+ibuf_len);
}

MH_io::~MH_io()
{delete [] obuf; delete [] ibuf;}

void MH_io::read_fisml(char *ibuf_)
{
  int int_buf[3], return_check;
  strcpy(ibuf_, ".fisml"); FILE * ref_sml = fopen(ibuf, "r");
  return_check=fscanf(ref_sml,"%d %d %d\n",int_buf,int_buf+1,int_buf+2);
  return_check=fscanf(ref_sml,"%d %d %d\n",int_buf,int_buf+1,int_buf+2); Frames=int_buf[2];
  return_check=fscanf(ref_sml,"%d %d %d\n",int_buf,int_buf+1,int_buf+2); nbeads=int_buf[1];
  fclose(ref_sml);
}

void MH_io::load_reference(double *ts_, double *xs_, double *d_ang_, double * comega_s_, double t_phys_)
{
  sprintf(ibuf+ibuf_len, ".filin");
  FILE * inputs = fopen(ibuf, "r");
  fread_SAFE(ts_, sizeof(double), Frames, inputs);
  fread_SAFE(xs_, sizeof(double), 2*Frames*nbeads, inputs);
  fread_SAFE(d_ang_, sizeof(double), Frames, inputs);
  fclose(inputs);

  if (noise_data)
  {
    MH_rng ran(1);
    for (int i = 2*nbeads; i < 2*Frames*nbeads; i++) xs_[i]+=ran.rand_gau(0.0,noise_sigma);
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

// swirl_system_struct

swirl_system_struct::swirl_system_struct(int ulen_, int nbeads_, int Frames_): ulen(ulen_), nbeads(nbeads_), Frames(Frames_) {}

// record_struct

record_struct::record_struct(int ulen_, int nbeads_, int Frames_, int ichunk_len_, int dchunk_len_) : swirl_system_struct(ulen_, nbeads_, Frames_), ichunk_len(ichunk_len_), dchunk_len(dchunk_len_) {}

// thread_worker_struct

thread_worker_struct::thread_worker_struct(int ulen_, int nbeads_, int Frames_, int nlead_, int npool_, double dt_sim_, double t_phys_, double *ts_, double *xs_, double *d_ang_, double *comega_s_): swirl_system_struct(ulen_, nbeads_, Frames_),
nlead(nlead_), npool(npool_), dt_sim(dt_sim_), t_phys(t_phys_), ts(ts_), xs(xs_), d_ang(d_ang_), comega_s(comega_s_) {}

// MH_params

MH_params::MH_params(MH_io *io_, int ulen_, int nlead_, int npool_, double dt_sim_, double t_phys_, double sigma_): swirl_system_struct(ulen_, io_->nbeads, io_->Frames),
io(io_),
nlead(nlead_), npool(npool_), dt_sim(dt_sim_), t_phys(t_phys_), sigma(sigma_) {}

// MH_train_struct

MH_train_struct::MH_train_struct(MH_params *par_, swirl_param *sp_min_, swirl_param *sp_max_, wall_list *wl_): par(par_), sp_min(sp_min_), sp_max(sp_max_), wl(wl_) {}
