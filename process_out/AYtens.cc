#include <cstdio>

#include "AYlinalg.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
}

AYtens::AYtens(int W_, int N_, int M_): W(W_), M(M_), N(N_)
{
  M_ptr = (AYmat**)malloc((W_)*sizeof(AYmat*));

}
AYtens::~AYtens()
{
  for (int i = 0; i < W; i++) delete *(M_ptr+i);
  free(M_ptr);
}
void AYtens::fprintf_tens(char name[], bool verbose)
{
  int i;
  char specfile[300]; memset(specfile, 0, 299); snprintf(specfile, 300, "%s.aytens", name);
  char smlfile[300]; memset(smlfile, 0, 299); snprintf(smlfile, 300, "%s.aysml", name);
  FILE * data_file = fopen(specfile, "wb");
  FILE * aysml_file = fopen(smlfile, "w");
  fprintf(aysml_file, "%d %d %d");
  fclose(aysml_file);
  for ( i = 0; i < W; i++) fwrite(M_ptr[i]->A_ptr, sizeof(double), M*N, data_file);
  fclose(data_file);
  if (verbose) printf("wrote file: %s.aytens/aysml\n", name);
}
