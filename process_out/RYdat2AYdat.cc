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

ODR_struct::ODR_struct(char name_[], int Frames_): AYdata(Frames_, 2)
{
    set_dims();
    int len_specs = dims[0][1];
    int len_dat = dims[1][1];
    memset(name, 0, 49); snprintf(name, 50, "%s", name_);
    memset(rydat_dir, 0, 74); snprintf(rydat_dir, 75, "../twin/%s", name);
    char file_it[100]; memset(file_it, 0, 99); snprintf(file_it, 100, "%spts.0", rydat_dir);

    double t, cx, cy, wall_sca;
    double x, y, z, q1, q2, q3, q4;

    std::ifstream ry_dat;
    ry_dat.open(file_it);
    int lines = 0;
    for(std::string line; std::getline(ry_dat, line);) ++lines; // determining the number of line in the file
    P=lines-1;
    ry_dat.close();

    data = new AYtens(Frames, len_dat, P);
    specs = new AYmat(len_specs, Frames);

    for ( int i = 0; i < Frames; i++)
    {
      memset(file_it, 0, 99); snprintf(file_it, 100, "%spts.%d", rydat_dir, i);
      std::ifstream ry_dat;
      ry_dat.open(file_it);
      int count = 0;
      for(std::string line; std::getline(ry_dat, line);)
      {
        std::istringstream in(line);
        if (count == 0)
        {
          in >> t >> cx >> cy >> wall_sca;
          specs->set(0, i, t);
          specs->set(1, i, cx);
          specs->set(2, i, cy);
          specs->set(3, i, wall_sca);
        }
        else
        {
          in >> x >> y >> z >> q1 >> q2 >> q3 >> q4;
          data->set(0, count-1, i, x);
          data->set(1, count-1, i, y);
          data->set(2, count-1, i, z);
          data->set(3, count-1, i, q1);
          data->set(4, count-1, i, q2);
          data->set(5, count-1, i, q3);
          data->set(6, count-1, i, q4);
        }
        count++;
      }
      ry_dat.close();
    }
}

ODR_struct::~ODR_struct()
{
  delete data; delete specs;
}

void ODR_struct::set_dims()
{
  dims[0][0] = 1; dims[0][1] = 4; dims[0][2] = 1;
  dims[1][0] = 1; dims[1][1] = 7; dims[1][2] = P;
}


void ODR_struct::fprintf_split(char name_[], bool verbose_)
{
  char aysml_name[300]; memset(aysml_name, 0, 299); snprintf(aysml_name, 300, "%s%s", directory, name_);
  char specfile[300];

  for (int i = 0; i < Frames; i++)
  {
    memset(specfile, 0, 299); snprintf(specfile, 300, "%s%s_%d.aydat", directory, name_, i);
    FILE * data_file = fopen(specfile, "wb");
    fwrite(specs->AT[i], sizeof(double), dims[0][1], data_file);
    fwrite(data->T_AT[i], sizeof(double), dims[1][1]*dims[1][2], data_file);
    fclose(data_file);
  }
  AYdata_aysml_gen(aysml_name, 1);
  if (verbose_) printf("AYdata: wrote file(s) %s.aydat/aysml\n", name);
}

void ODR_struct::prepare_datdir(char name_[])
{memset(directory, 0, 99); snprintf(directory, 100, "%s%s", name_, name); mkdir(directory, S_IRWXU);
printf("made directory  %s\n", directory);}
