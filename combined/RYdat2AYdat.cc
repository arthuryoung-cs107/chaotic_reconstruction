#include <sys/stat.h>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>

#include "RYdat2AYdat.hh"

extern "C"
{
  #include "AYaux.h"
}

ODR_struct::ODR_struct(char *rydat_loc_, char *rydat_dir_, char *file_name_, int Frames_): AYdata(Frames_, 2), rydat_loc(string_gen_pruned(rydat_loc_)), rydat_dir(string_gen_pruned(rydat_dir_)), file_name(string_gen_pruned(file_name_)), rydat_loc_alloc_flag(true),
rydat_dir_alloc_flag(true),
file_name_alloc_flag(true)
{
    size_t in_buf_len = (size_t)(strlen(rydat_loc)+strlen(rydat_dir)+strlen(file_name) + 50);
    in_buf = new char[in_buf_len]; sprintf(in_buf, "%s%s%s", rydat_loc, rydat_dir, file_name);

    char *file_it = in_buf + strlen(in_buf); //address at the end of the input string
    sprintf(file_it, ".%d", 0);

    std::ifstream ry_dat;
    ry_dat.open(in_buf);
    int lines = 0;
    for(std::string line; std::getline(ry_dat, line);) ++lines; // determining the number of line in the file
    P=lines-1;
    ry_dat.close();
    set_dims();

    data = AYd3tensor(Frames, P, len_dat); data_alloc_flag = true;
    specs = AYdmatrix(Frames, len_specs); specs_alloc_flag = true;

    double t, cx, cy, wall_sca;
    double x, y, z, q1, q2, q3, q4;

    for ( int i = 0; i < Frames; i++)
    {
      sprintf(file_it, ".%d", i);
      std::ifstream ry_dat;
      ry_dat.open(in_buf);
      int count = 0;
      for(std::string line; std::getline(ry_dat, line);)
      {
        std::istringstream in(line);
        if (count == 0)
        {
          in >> t >> cx >> cy >> wall_sca;
          specs[i][0] = t;
          specs[i][1] = cx;
          specs[i][2] = cy;
          specs[i][3] = wall_sca;
        }
        else
        {
          in >> x >> y >> z >> q1 >> q2 >> q3 >> q4;
          data[i][count-1][0] = x;
          data[i][count-1][1] = y;
          data[i][count-1][2] = z;
          data[i][count-1][3] = q1;
          data[i][count-1][4] = q2;
          data[i][count-1][5] = q3;
          data[i][count-1][6] = q4;
        }
        count++;
      }
      ry_dat.close();
    }
}

ODR_struct::ODR_struct(char * proc_loc_, char * rydat_dir_, char * file_name_): AYdata(), rydat_dir(string_gen_pruned(rydat_dir_)), file_name(string_gen_pruned(file_name_)), rydat_dir_alloc_flag(true), file_name_alloc_flag(true)
{
  size_t proc_len = (size_t)(strlen(proc_loc_)+strlen(rydat_dir)+1);
  proc_dir = new char[proc_len]; sprintf(proc_dir, "%s%s", proc_loc_, rydat_dir); proc_dir_alloc_flag=true; mkdir(proc_dir, S_IRWXU);

  size_t out_buf_len = (size_t)(strlen(proc_dir)+strlen(file_name)+50);
  out_buf = new char[out_buf_len]; sprintf(out_buf, "%s%s", proc_dir, file_name);

  Frames = 0; depth = 2; dims = AYimatrix(depth, 3); dims_alloc_flag = true;
}

ODR_struct::~ODR_struct()
{
  if (data_alloc_flag) free_AYd3tensor(data);
  if (specs_alloc_flag) free_AYdmatrix(specs);
  if (ibuf_alloc_flag) delete in_buf;
  if (obuf_alloc_flag) delete out_buf;
  if (rydat_loc_alloc_flag) delete rydat_loc;
  if (rydat_dir_alloc_flag) delete rydat_dir;
  if (file_name_alloc_flag) delete file_name;
  if (proc_dir_alloc_flag) delete proc_dir;
}

void ODR_struct::set_dims()
{
  dims[0][0] = 1; dims[0][1] = 1; dims[0][2] = len_specs;
  dims[1][0] = 1; dims[1][1] = P; dims[1][2] = len_dat;
}

void ODR_struct::fprintf_split(bool verbose_)
{
  char * file_it = out_buf + strlen(out_buf);

  for (int i = 0; i < Frames; i++)
  {
    sprintf(file_it, ".%d.aydat", i);
    FILE * data_file = fopen(out_buf, "wb");
    fwrite(specs[i], sizeof(double), dims[0][1]*dims[0][2], data_file);
    fwrite(data[i][0], sizeof(double), dims[1][1]*dims[1][2], data_file);
    fclose(data_file);
  }

  sprintf(file_it, ".rysml");
  FILE * aysml_file = fopen(out_buf, "w");
  fprintf(aysml_file, "%d %d %d\n", 1, Frames, depth);
  for (int i = 0; i < depth; i++) fprintf(aysml_file, "%d %d %d\n", dims[i][0], dims[i][1], dims[i][2]);
  fclose(aysml_file);

  if (verbose_) printf("ODR_struct: wrote file(s) %s%s.aydat/rysml\n", proc_dir, file_name);
}

void ODR_struct::prepare_datdir(char *proc_loc_)
{
  size_t proc_len = (size_t)(strlen(proc_loc_)+strlen(rydat_dir)+1);
  proc_dir = new char[proc_len]; sprintf(proc_dir, "%s%s", proc_loc_, rydat_dir); proc_dir_alloc_flag=true; mkdir(proc_dir, S_IRWXU);
  printf("made directory  %s\n", proc_dir);
}

FILE * ODR_struct::write_split(int k_, particle * q_)
{
  char * file_it = out_buf + strlen(out_buf);

  if (writing_flag)
  {
    sprintf(file_it, ".%d.aydat", k_);
    FILE * data_out = fopen(out_buf, "wb");
    Frames++;
    return data_out;

    // fwrite(specs[i], sizeof(double), dims[0][1]*dims[0][2], data_file);
    // fwrite(data[i][0], sizeof(double), dims[1][1]*dims[1][2], data_file);
    // fclose(data_file);


  }
}
