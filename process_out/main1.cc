#include "RYdat2AYdat.hh"

extern "C"
{
  #include "AYaux.h"
}
void process_circ6()
{
  char rydat_dir[50]; name_gen(rydat_dir, 50, "circ6.odr/");
  char dat_dir[50]; name_gen(dat_dir, 50, "dat_dir/");
  char file_name[50]; name_gen(file_name, 50, "pts.");
  ODR_struct circ6odr(rydat_dir);

  circ6odr.specs->print_mat();
  circ6odr.prepare_datdir(dat_dir);
  circ6odr.fprintf_split(file_name);
}

int main()
{
  process_circ6();

  return 0;
}
