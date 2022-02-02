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
