#include "RYdat2AYdat.hh"

extern "C"
{
  #include "AYaux.h"
}
void process_circ6()
{
  char rydat_loc[] = "../twin/";
  char rydat_dir[] = "circ6.odr/";
  char proc_loc[] = "dat_dir/";
  char file_name[] = "pts";
  ODR_struct circ6odr(rydat_loc, rydat_dir, file_name, 1201);

  circ6odr.prepare_datdir(proc_loc);
  circ6odr.fprintf_split();
}

int main()
{
  process_circ6();

  return 0;
}
