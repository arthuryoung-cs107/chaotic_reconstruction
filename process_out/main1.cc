#include "RYdat2AYdat.hh"

extern "C"
{
  #include "AYaux.h"
}
void process_circ6()
{
  char rydat_dir[50]; name_gen(rydat_dir, 50, "circ6.odr/");
  ODR_struct circ6odr(rydat_dir);

  circ6odr.specs->print_mat();
}

int main()
{
  process_circ6();

  return 0;
}
