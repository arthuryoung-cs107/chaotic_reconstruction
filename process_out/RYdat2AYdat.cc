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
    int len_dat = dims[1][2];
    int len_specs = dims[0][1];
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

    data = new AYtens(Frames, P, len_dat);
    specs = new AYmat(Frames, len_specs);

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
          specs->set(i, 0, t);
          specs->set(i, 1, cx);
          specs->set(i, 2, cy);
          specs->set(i, 3, wall_sca);
        }
        else
        {
          in >> x >> y >> z >> q1 >> q2 >> q3 >> q4;
          data->set(count-1, 0, i, x);
          data->set(count-1, 1, i, y);
          data->set(count-1, 2, i, z);
          data->set(count-1, 3, i, q1);
          data->set(count-1, 4, i, q2);
          data->set(count-1, 5, i, q3);
          data->set(count-1, 6, i, q4);
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
  dims[1][0] = 1; dims[1][1] = P; dims[1][2] = 7;
}
