#include <sys/stat.h>
#include <sstream>
#include <string>
#include <fstream>

#include "blas.h"
#include "omp.h"
#include "AYlinalg.hh"
// #include "AYtens.hh"

extern "C"
{
  #include "nrutil.h"
  #include "knuth_lcg.h"
  #include "auxiliary_functions.h"
}

void preliminary_test1()
{
  AYmat mat1(3, 5);
  mat1.init_123();
  mat1.print_mat();
  char name1[50]; name_gen(name1, 50, "mat1");
  char name2[50]; name_gen(name2, 50, "mat2");
  mat1.fprintf_mat(name1);

  AYmat * mat2 = aysml_read(name2);
  mat2->print_mat();
}

void convert_RYdata()
{
  char * check;
  int len = 1200 + 1;
  char rydat_dir[50]; name_gen(rydat_dir, 50, "../twin/circ6.odr/");
  char specs_name[100]; name_gen(specs_name, 100, "dat_dir/cir6_specs");
  char tens_name[100]; name_gen(tens_name, 100, "dat_dir/cir6_data");

  char file_it[100];
  double t, cx, cy, wall_sca;
  double x, y, z, q1, q2, q3, q4;

  AYmat specs_mat(len, 4);

  std::ifstream ry_dat;
  ry_dat.open(file_it);
  int lines = 0;
  for(std::string line; std::getline(ry_dat, line);) ++lines;
  int dat_lines=lines-1;
  ry_dat.close();

  AYtens rytensor(len, dat_lines, 7);

  for ( int i = 0; i < 2; i++)
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
        specs_mat.set(i, 0, t);
        specs_mat.set(i, 1, cx);
        specs_mat.set(i, 2, cy );
        specs_mat.set(i, 3, wall_sca);
      }
      else
      {
        in >> x >> y >> z >> q1 >> q2 >> q3 >> q4;
        rytensor.M_ptr[i]->set(count-1, 0, x);
        rytensor.M_ptr[i]->set(count-1, 1, y);
        rytensor.M_ptr[i]->set(count-1, 2, z);
        rytensor.M_ptr[i]->set(count-1, 3, q1);
        rytensor.M_ptr[i]->set(count-1, 4, q2);
        rytensor.M_ptr[i]->set(count-1, 5, q3);
        rytensor.M_ptr[i]->set(count-1, 6, q4);
      }
      count++;
    }
    ry_dat.close();
  }
  specs_mat.fprintf_mat(specs_name);
  rytensor.fprintf_tens(tens_name);
}

void process_circ6()
{

}

int main()
{
  // preliminary_test1();
  // process_circ6();
  convert_RYdata();
  return 0;
}
