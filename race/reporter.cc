#include <sys/stat.h>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>

#include "particle_race.hh"

void reporter::init_race( char * proc_loc_, char * rydat_dir_, char * file_name_, int race_id_): race_id(race_id_)
{
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

void reporter::write_gen_diagnostics(int gen_count_, int leader_count_, int worst_leader_, int best_leader_)
{
  int int_params[] = {gen_count_, leader_count_, record_int_len, record_double_len, len, worst_leader_, best_leader_};

  char * buf_it = out_buf+obuf_end;
  sprintf(buf_it, "%s.rc%d_gen%d.rcdat", file_name, race_id, gen_count_);
  FILE * out_dat = fopen(out_buf, "wb");
  fwrite(int_params, sizeof(int), 7, out_dat);
  fwrite(sample_weights, sizeof(double), leader_count_, out_dat);
  for (int i = 0; i < leader_count_; i++)
  {
    // write the integer records
    fwrite(&(leaders[i]->global_index), sizeof(int), record_int_len, out_dat);
    // write the integer records
    fwrite(&(leaders[i]->l2score), sizeof(double), record_double_len, out_dat);
    // write the parameter
    fwrite(leaders[i]->params, sizeof(double), len, out_dat);
  }
  fclose(out_dat);
}

void reporter::close_diagnostics(int gen_count_, int leader_count_, int worst_leader_, int best_leader_)
{
  int int_params[] = {gen_count_, leader_count_, record_int_len, record_double_len, len, worst_leader_, best_leader_};

  char * buf_it = out_buf+obuf_end;
  sprintf(buf_it, "%s.rc%d_end.rcdat", file_name, race_id);
  FILE * out_dat = fopen(out_buf, "wb");
  fwrite(int_params, sizeof(int), 7, out_dat);
  for (int i = 0; i < leader_count_; i++)
  {
    // write the integer records
    fwrite(&(leaders[i]->global_index), sizeof(int), record_int_len, out_dat);
    // write the integer records
    fwrite(&(leaders[i]->l2score), sizeof(double), record_double_len, out_dat);
    // write the parameter
    fwrite(leaders[i]->params, sizeof(double), len, out_dat);
  }
  fclose(out_dat);
}

ODR_struct * reporter::spawn_swirlODR(char * name_)
{
  strcpy(out_buf + obuf_end, "");
  ODR_struct * odr_out = new ODR_struct();
  odr_out->init_swirl(out_buf, name_, file_name, false);
  return odr_out;
}
