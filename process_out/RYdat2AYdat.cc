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

ODR_struct::ODR_struct(char name_[], int len_)
{
    len = len_;

    memset(name, 0, 49); snprintf(name, 50, "%s", name_);
    memset(rydat_dir, 0, 74); snprintf(rydat_dir, 75, "../twin/%s", name);
    char file_it[100]; memset(file_it, 0, 99); snprintf(file_it, 100, "%spts.0", rydat_dir);

    double t, cx, cy, wall_sca;
    double x, y, z, q1, q2, q3, q4;

    std::ifstream ry_dat;
    ry_dat.open(file_it);
    int lines = 0;
    for(std::string line; std::getline(ry_dat, line);) ++lines; // determining the number of line in the file
    int dat_lines=lines-1;
    ry_dat.close();

    printf("len, dat_lines: %d %d \n", len, dat_lines);

    data = new AYtens(len, dat_lines, 7);
    specs = new AYmat(len, 4);

    for ( int i = 0; i < len; i++)
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
          specs->set(i, 2, cy );
          specs->set(i, 3, wall_sca);
        }
        else
        {
          in >> x >> y >> z >> q1 >> q2 >> q3 >> q4;
          data->mat[i].set(count-1, 0, x);
          data->mat[i].set(count-1, 1, y);
          data->mat[i].set(count-1, 2, z);
          data->mat[i].set(count-1, 3, q1);
          data->mat[i].set(count-1, 4, q2);
          data->mat[i].set(count-1, 5, q3);
          data->mat[i].set(count-1, 6, q4);
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
