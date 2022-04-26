#include <sys/stat.h>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>

#include "particle_relay.hh"

void reporter::init_relay( char * proc_loc_, char * rydat_dir_, char * file_name_, int relay_id_)
{
  relay_id = relay_id_;
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

void reporter::write_startup_diagnostics()
{
  int header[] = {9,5}; // {int_len, double_len}
  int int_params[] = {nlead, npool, param_len, beads, Frames, record_int_len, record_double_len, record_int_chunk_count, record_double_chunk_count};
  double double_params[] = {dt_sim, noise_tol, alpha_tol, max_weight_ceiling, t_phys};

  char * buf_it = out_buf+obuf_end; sprintf(buf_it, "%s.re%d_startup.redat", file_name, relay_id, gen_count_);
  FILE * out_dat = fopen(out_buf, "wb");
  fwrite(header, sizeof(int), 2, out_dat);
  fwrite(int_params, sizeof(int), header[0], out_dat);
  fwrite(double_params, sizeof(double), header[1], out_dat);
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
  fwrite(event_end, sizeof(int), beads, out_dat);
  for (int i = 0; i < leader_count_; i++)
  {
    // write the parameter
    fwrite(leaders[i]->params, sizeof(double), param_len, out_dat);
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
    fwrite(&(leaders[i]->residual), sizeof(double), 1, out_dat);
    // write the parameter

    fwrite(leaders[i]->residual_data, sizeof(double), beads, out_dat);

    fwrite(leaders[i]->params, sizeof(double), param_len, out_dat);
  }
  fclose(out_dat);
}

void reporter::close_diagnostics(int gen_count_)
{
  int header[] = {gen_int_len,gen_double_len}; // {int_len, double_len}

  char * buf_it = out_buf+obuf_end; sprintf(buf_it, "%s.re%d_end.redat", file_name, relay_id, gen_count_);
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
