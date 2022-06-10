#include "MH_auxiliary.hh"
#include "MH_tools.hh"

int set_special_params(int id_, double *vec_)
{
  int id_out=0;
  if ((id_>=0)&&(id_<special_param_count)) id_out=id_;
  else printf("WARNING: invalid parameter ID (%d). Using default parameters\n", id_);

  for (int i = 0; i < full_param_len; i++) vec_[i] = special_parameters[id_out][i];
  return id_out;
}

int set_special_params(const char *id_, double*vec_)
{
  if (strcmp(id_, "true")==0) return set_special_params(0, dvec_);
  else if (strcmp(id_, "min")==0) return set_special_params(1, dvec_);
  else if (strcmp(id_, "max")==0) return set_special_params(2, dvec_);
  else if (strcmp(id_, "pert")==0) return set_special_params(3, dvec_);
  else if (strcmp(id_, "b3i0r5")==0) return set_special_params(4, dvec_);
  else return set_special_params(-1, dvec_);
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
