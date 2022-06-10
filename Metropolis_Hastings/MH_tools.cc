#include "MH_tools.hh"

template <class T> T ** Tmatrix(int M_, int N_)
{
  T * chunk = new T[M_*N_],
    ** rows = new T*[M_];
  for (int i = 0,j=0; i < M_; i++,j+=N_)
    rows[i] = chunk+j;
  return rows;
}

template <class T> void free_Tmatrix(T ** Tmat_)
{
  delete [] Tmat_[0];
  delete [] Tmat_;
}
