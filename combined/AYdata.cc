#include "AYlinalg.hh"

extern "C"
{
  #include "AYaux.h"
}

AYdata::AYdata(int Frames_, int depth_ ): Frames(Frames_), depth(depth_), dims(AYimatrix(depth_, 3))
{}

AYdata::~AYdata()
{free_AYimatrix(dims);}

void AYdata::set_dims()
{}

void AYdata::AYdata_aysml_gen(char name_[], int split_)
{
  char smlfile[300]; memset(smlfile, 0, 299); snprintf(smlfile, 300, "%s.aysml", name_);
  FILE * aysml_file = fopen(smlfile, "w");
  fprintf(aysml_file, "%d %d %d \n", split_, Frames, depth);
  for (int i = 0; i < depth; i++)
  {
    fprintf(aysml_file, "%d %d %d\n", dims[i][0], dims[i][1], dims[i][2]);
  }
  fclose(aysml_file);
}

void AYdata::fprintf_split(char name_[], bool verbose_)
{}
