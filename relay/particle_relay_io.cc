#include <sys/stat.h>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>

#include "particle_relay.hh"

extern "C"
{
  #include "AYaux.h"
}

void reporter::init_relay( char * proc_loc_, char * rydat_dir_, char * file_name_, int relay_id_, bool noise_data_, double noise_tol_in_)
{
  relay_id = relay_id_; noise_data = noise_data_; noise_tol_in = noise_tol_in_;
  reading_flag=true;
  rydat_dir=string_gen_pruned(rydat_dir_); file_name=string_gen_pruned(file_name_);
  size_t len_loc = strlen(proc_loc_) + strlen(rydat_dir);
  out_buf = new char[len_loc+100]; filin_dir= new char[len_loc+strlen(file_name)+1];

  sprintf(out_buf, "%s%s", proc_loc_, rydat_dir_); obuf_end = strlen(out_buf);
  sprintf(filin_dir, "%s%s%s", proc_loc_, rydat_dir_, file_name_);

  size_t filin_dir_len = (size_t)(strlen(filin_dir) + 50);
  in_buf = new char[filin_dir_len]; strcpy(in_buf, filin_dir);
  ibuf_end = strlen(in_buf);

  char * file_it = in_buf + ibuf_end;
  sprintf(file_it, ".fisml");
  int val1, val2, val3;
  std::ifstream sml_stream; std::string line; sml_stream.open(in_buf);
  std::getline(sml_stream, line); std::istringstream dimline0(line);
  dimline0 >> val1 >> depth >> val3;
  std::getline(sml_stream, line); std::istringstream dimline1(line);
  dimline1 >> val1 >> val2 >> Frames;
  std::getline(sml_stream, line); std::istringstream dimline2(line);
  dimline2 >> val1 >> P >> val3;

  sml_stream.close();
}

void reporter::load_relay(double *ts_, double *xs_, double *d_ang_)
{
  if (reading_flag)
  {
    int nsnaps = Frames;
    char * file_name = in_buf + ibuf_end;
    sprintf(file_name, ".filin");
    FILE * inputs = fopen(in_buf, "r"); // or "rb"?
    fread_safe(ts_, sizeof(double), nsnaps, inputs);
    fread_safe(xs_, sizeof(double), 2*nsnaps*P, inputs);
    fread_safe(d_ang_, sizeof(double), nsnaps, inputs);
    fclose(inputs);

    if (noise_data)
    {
      AYrng ran;
      ran.rng_init_gsl(1);
      for (int i = 2*P; i < 2*nsnaps*P; i++) xs_[i]+=ran.rand_gau_gsl(0.0,noise_tol_in);
    }
  }
  else printf("ODR_struct: not staged for reading binary filter inputs\n");
}

void reporter::write_startup_diagnostics(int *header_, int *int_params_, double *double_params_)
{
  char * buf_it = out_buf+obuf_end; sprintf(buf_it, "%s.re%d_startup.redat", file_name, relay_id);
  FILE * out_dat = fopen(out_buf, "wb");
  fwrite(header_, sizeof(int), 2, out_dat);
  fwrite(int_params_, sizeof(int), header_[0], out_dat);
  fwrite(double_params_, sizeof(double), header_[1], out_dat);
  fclose(out_dat);
}

void reporter::write_event_diagnostics(int event_)
{
  char * buf_it = out_buf+obuf_end; sprintf(buf_it, "%s.re%d_event%d.redat", file_name, relay_id, event_);
  FILE * out_dat = fopen(out_buf, "wb");

  fwrite(global_event_frame_count, sizeof(int), beads*Frames, out_dat);
  for (int i = 0; i < npool; i++)
  {
    // write the integer records
    fwrite(&(pool[i]->global_index), sizeof(int), record_int_len, out_dat);
    // write the integer records
    fwrite(&(pool[i]->residual), sizeof(double), record_double_len, out_dat);

    fwrite(pool[i]->event_positions, sizeof(int), beads*record_int_chunk_count, out_dat);
    fwrite(pool[i]->residual_data, sizeof(double), beads*record_double_chunk_count, out_dat);

    // write the parameter
    fwrite(pool[i]->params, sizeof(double), param_len, out_dat);
  }
  fclose(out_dat);
}

void reporter::write_postevent_diagnostics(int event_)
{
  int header[] = {event_int_len,event_double_len}; // {int_len, double_len}

  char * buf_it = out_buf+obuf_end; sprintf(buf_it, "%s.re%d_postevent%d.redat", file_name, relay_id, event_);
  FILE * out_dat = fopen(out_buf, "wb");

  fwrite(header, sizeof(int), 2, out_dat);
  fwrite(postevent_int_vec, sizeof(int), header[0], out_dat);
  fwrite(postevent_double_vec, sizeof(double), header[1], out_dat);
  fwrite(event_frames, sizeof(int), beads, out_dat);
  for (int i = 0; i < npool; i++)
  {
    // write the parameter
    fwrite(pool[i]->params, sizeof(double), param_len, out_dat);
  }
  fclose(out_dat);
}

void reporter::write_gen_diagnostics(int gen_count_, int leader_count_)
{
  int header[] = {gen_int_len,gen_double_len}; // {int_len, double_len}

  char * buf_it = out_buf+obuf_end; sprintf(buf_it, "%s.re%d_gen%d.redat", file_name, relay_id, gen_count_);
  FILE * out_dat = fopen(out_buf, "wb");

  fwrite(header, sizeof(int), 2, out_dat);
  fwrite(gen_int_vec, sizeof(int), header[0], out_dat);
  fwrite(gen_double_vec, sizeof(double), header[1], out_dat);
  fwrite(lead_dup_count, sizeof(int), leader_count_, out_dat);
  fwrite(sample_weights, sizeof(double), leader_count_, out_dat);
  fwrite(lead_par_w_mean, sizeof(double), param_len, out_dat);
  fwrite(lead_par_w_var, sizeof(double), param_len, out_dat);
  for (int i = 0; i < leader_count_; i++)
  {
    // write the integer records
    fwrite(&(leaders[i]->global_index), sizeof(int), record_int_len, out_dat);
    // write the integer records
    fwrite(&(leaders[i]->residual), sizeof(double), record_double_len, out_dat);
    // write the parameter

    fwrite(leaders[i]->residual_data, sizeof(double), beads, out_dat);

    fwrite(leaders[i]->params, sizeof(double), param_len, out_dat);
  }
  fclose(out_dat);
}

void reporter::close_diagnostics(int gen_count_)
{
  int header[] = {gen_int_len,gen_double_len}; // {int_len, double_len}

  char * buf_it = out_buf+obuf_end; sprintf(buf_it, "%s.re%d_end.redat", file_name, relay_id);
  FILE * out_dat = fopen(out_buf, "wb");

  fwrite(header, sizeof(int), 2, out_dat);
  fwrite(gen_int_vec, sizeof(int), header[0], out_dat);
  fwrite(gen_double_vec, sizeof(double), header[1], out_dat);
  fclose(out_dat);
}

ODR_struct * reporter::spawn_swirlODR(char * name_)
{
  strcpy(out_buf + obuf_end, "");
  ODR_struct * odr_out = new ODR_struct();
  odr_out->init_swirl(out_buf, name_, file_name, false);
  return odr_out;
}
