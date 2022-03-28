#include "AYlinalg.hh"

extern "C"
{
  #include "AYaux.h"
}

AYdata::AYdata(int Frames_, int depth_ ): Frames(Frames_), depth(depth_), dims(AYimatrix(depth_, 3))
{}

AYdata::AYdata()
{}

AYdata::~AYdata()
{if (dims!=NULL) free_AYimatrix(dims);}

void AYdata::set_dims()
{}

void AYdata::fprintf_split(char name_[], bool verbose_)
{}
