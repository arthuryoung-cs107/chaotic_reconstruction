#include <sys/stat.h>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>

#include "particle_walk.hh"

void pedestrian::init_walk( char * proc_loc_, char * rydat_dir_, char * file_name_, int walk_id_)
{
  walk_id = walk_id_;
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

void pedestrian::write_gen_diagnostics(int gen_count_, int leader_count_, int worst_leader_, int best_leader_, double t_wheels_, double min_res_)
{
  int header[] = {8,2}; // {int_len, double_len}
  int int_params[] = {gen_count_, Frames, leader_count_, grade_int_len, grade_double_len, param_len, worst_leader_, best_leader_};
  double double_params[] = {t_wheels_, min_res_};

  char * buf_it = out_buf+obuf_end; sprintf(buf_it, "%s.wa%d_gen%d.wadat", file_name, walk_id, gen_count_, min_res_);
  FILE * out_dat = fopen(out_buf, "wb");

  fwrite(header, sizeof(int), 2, out_dat);
  fwrite(int_params, sizeof(int), header[0], out_dat);
  fwrite(double_params, sizeof(double), header[1], out_dat);
  fwrite(lead_dup_count, sizeof(int), leader_count_, out_dat);
  fwrite(sample_weights, sizeof(double), leader_count_, out_dat);
  fwrite(frame_kill_count, sizeof(int), Frames, out_dat);
  fwrite(gen_frame_res_data, sizeof(double), 4*Frames, out_dat);
  fwrite(gen_param_mean, sizeof(double), param_len, out_dat);
  fwrite(gen_param_var, sizeof(double), param_len, out_dat);
  for (int i = 0; i < leader_count_; i++)
  {
    // write the integer grades
    fwrite(&(leaders[i]->global_index), sizeof(int), grade_int_len, out_dat);
    // write the integer grades
    fwrite(&(leaders[i]->l2score), sizeof(double), grade_double_len, out_dat);
    // write the parameter
    fwrite(leaders[i]->params, sizeof(double), param_len, out_dat);
  }
  fclose(out_dat);
}

void pedestrian::close_diagnostics(int gen_count_, int worst_leader_, int best_leader_, double t_wheels_, double min_res_)
{
  int header[] = {6,2}; // {int_len, double_len}
  int int_params[] = {gen_count_, nlead, npool, nA, worst_leader_, best_leader_};
  double double_params [] = {t_wheels_, min_res_};

  char * buf_it = out_buf+obuf_end;
  sprintf(buf_it, "%s.wa%d_end.wadat", file_name, walk_id);
  FILE * out_dat = fopen(out_buf, "wb");
  fwrite(header, sizeof(int), 2, out_dat);
  fwrite(int_params, sizeof(int), header[0], out_dat);
  fwrite(double_params, sizeof(double), header[1], out_dat);
  fclose(out_dat);
}

ODR_struct * pedestrian::spawn_swirlODR(char * name_)
{
  strcpy(out_buf + obuf_end, "");
  ODR_struct * odr_out = new ODR_struct();
  odr_out->init_swirl(out_buf, name_, file_name, false);
  return odr_out;
}
