#include <cstdio>

#include "AYlinalg.hh"

extern "C"
{
  #include "AYaux.h"
}


AYvec::AYvec(int M_): AYmat(M_, 1) {}
AYvec::~AYvec() {}

void AYvec::set(int i, double val)
{
  if ((i >= 0)&&(i<M)) A_ptr[i] = val;
  else printf("AYvec: set failed, bad dimensions\n");
}
double AYvec::get(int i)
{
  if ((i >= 0)&&(i<M)) return A_ptr[i];
  else {printf("AYvec: get failed, bad dimensions\n"); return 0;}
}

void AYvec::GSL_2_AYvec_copy(gsl_matrix * mat_in) {GSL_2_AYmat_copy(mat_in);}
void AYvec::GSL_2_AYvec_copy(gsl_vector * vec_in) {GSL_2_AYmat_copy(vec_in);}
void AYvec::AYvec_2_GSL_copy(gsl_matrix * mat_in) {AYmat_2_GSL_copy(mat_in);}
void AYvec::AYvec_2_GSL_copy(gsl_vector * vec_in) {AYmat_2_GSL_copy(vec_in);}

void AYvec::print_vec(bool space_)
{print_mat(space_);}

AYvec * AYvec::copy_gen()
{
  AYvec * x_out = new AYvec(M);
  memcpy(x_out->A_ptr, A_ptr, M*sizeof(double));
  return x_out;
}

AYvec * AYvec::row_slice_gen(int i_)
{
  printf("AYvec: row_slice_gen warning, requesting a row slice of a vector. Will return a 1x1 AYvec\n");
  AYvec * v_out = new AYvec(1);
  *(v_out->A_ptr) = A_ptr[i_];
  return v_out;
}
AYvec * AYvec::col_slice_gen(int j_)
{
  if (j_ == 0)
  {
    printf("AYvec: col_slice_gen warning, requesting a column slice of a vector. Will return the same vector");
    return copy_gen();
  }
  else
  {
    printf("AYvec: col_slice_gen error, requesting j = %d column slice of a vector.  Will return the same vector", j_);
    return copy_gen();
  }
}
AYmat * AYvec::transpose_gen()
{
  printf("AYvec: transpose_gen warning, requesting the transpose of a vector. Will return 1x%d AYmatrix\n", M);
  int i, j;
  AYmat * X_out = new AYmat(N, M);
  for ( i = 0; i < M; i++) {for ( j = 0; j < N; j++) { X_out->A_ptr[i*N + j] = AT[j][i];}}
  return X_out;
}

double AYvec::norm_2()
{
  double out = 0.0;
  for (int i = 0; i < M; i++) out += (A_ptr[i]*A_ptr[i]);
  return sqrt(out);
}
double AYvec::dot(AYmat *B_)
{
  double out = 0.0;
  if ((B_->M ==  M)&&(B_->N == 1))
  {
    for (int i = 0; i < M; i++) out += B_->A_ptr[i]*A_ptr[i];
  }
  else printf("AYvec: dot failed, dim(A) = (%d %d) vs. dim(B) = (%d %d)\n", M, N, B_->M, B_->N);

  return out;
}
void AYvec::svd(gsl_vector *S, gsl_matrix *V, gsl_vector *work) {printf("AYvec: svd failed, cannot take SVD of a vector\n");}
void AYvec::svd_check() {printf("AYvec: svd_check failed, cannot take SVD of a vector\n");}
