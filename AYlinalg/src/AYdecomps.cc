#include "AYlinalg.hh"

extern "C"
{
  #include "AYaux.h"
}

AY_SVDspace::AY_SVDspace(AYmat * mat_): M_in(mat_->M), N_in(mat_->N), U(gsl_matrix_alloc(mat_->M, mat_->N)), V(gsl_matrix_alloc(mat_->N, mat_->N)), s(gsl_vector_alloc(mat_->N)), work(gsl_vector_alloc(mat_->N)) // assume long thin matrix for now
{}
AY_SVDspace::AY_SVDspace(int M_, int N_): M_in(M_), N_in(N_), U(gsl_matrix_alloc(M_, N_)), V(gsl_matrix_alloc(N_, N_)), s(gsl_vector_alloc(N_)), work(gsl_vector_alloc(N_)) // assume long thin matrix for now
{}
AY_SVDspace::~AY_SVDspace() { gsl_matrix_free(U); gsl_matrix_free(V); gsl_vector_free(s); gsl_vector_free(work); }

void AY_SVDspace::load_U(AYmat * X_) {for (int i = 0; i < M_in; i++) {for (int j = 0; j < N_in; j++) {gsl_matrix_set(U, i, j, X_->AT[j][i]);}}}
void AY_SVDspace::svd() {gsl_linalg_SV_decomp(U, V, s, work);}
void AY_SVDspace::unpack(AYmat * U_, AYmat * S_, AYmat * V_)
{
  for (int j = 0; j < N_in; j++)
  {
    for (int i = 0; i < M_in; i++)
    {
      U_->set(i, j, gsl_matrix_get(U, i, j));
      if (i < N_in)
      {
        V_->set(i, j, gsl_matrix_get(V, i, j));
        S_->set(i, j, 0.0);
      }
    }
    S_->set(j, j, gsl_vector_get(s, j));
  }
}

AY_Choleskyspace::AY_Choleskyspace(AYsym * mat_): N_in(mat_->N), mat_gsl(gsl_matrix_alloc(mat_->N, mat_->N))
{load_mat(mat_);}

AY_Choleskyspace::AY_Choleskyspace(int N_): N_in(N_), mat_gsl(gsl_matrix_alloc(N_, N_))
{}

AY_Choleskyspace::~AY_Choleskyspace()
{
  gsl_matrix_free(mat_gsl);
  if (x_gsl!=NULL) gsl_vector_free(x_gsl);
}

void AY_Choleskyspace::load_mat(AYsym * mat_)
{
  int i,j;
  for ( i = 0; i < N_in; i++)
  {
    //lower diagonal components and diagonal only, exploiting gsl technique
    for ( j = 0; j < i; j++) gsl_matrix_set(mat_gsl, i, j, mat_->A[j][i-j]);
    for ( j = i; j < N_in; j++) gsl_matrix_set(mat_gsl, i, j, mat_->A[i][j-i]);
  }
}
void AY_Choleskyspace::load_mat(AYsym * mat_, double scal_)
{
  for (int i = 0; i < N_in; i++)
  {
    //lower diagonal components and diagonal only, exploiting gsl technique
    for (int j = 0; j < i; j++) gsl_matrix_set(mat_gsl, i, j, scal_*mat_->A[j][i-j]);
    for (int j = i; j < N_in; j++) gsl_matrix_set(mat_gsl, i, j, scal_*mat_->A[i][j-i]);
  }
}

void AY_Choleskyspace::Cholesky_decomp()
{gsl_linalg_cholesky_decomp1(mat_gsl);}

void AY_Choleskyspace::Cholesky_decomp(AYsym * mat_, AYsym * L_)
{
  int i,j;
  for ( i = 0; i < N_in; i++)
  {
    //lower diagonal components and diagonal only, exploiting gsl technique
    for ( j = 0; j < i; j++) gsl_matrix_set(mat_gsl, i, j, mat_->A[j][i-j]);
    for ( j = i; j < N_in; j++) gsl_matrix_set(mat_gsl, i, j, mat_->A[i][j-i]);
  }
  gsl_linalg_cholesky_decomp1(mat_gsl);
  for ( i = 0; i < N_in; i++) for ( j = i; j < N_in; j++) L_->A[i][j-i] = gsl_matrix_get(mat_gsl, j, i);
}

void AY_Choleskyspace::iCholesky_decomp(AYsym * mat_, AYsym * L_, double threshold_)
{
  int i,j;
  double * l1_thresh = new double[N_in];
  for ( i = 0; i < N_in; i++)
  {
    l1_thresh[i] = 0.0;
    for ( j = 0; j < i; j++)
    {
      gsl_matrix_set(mat_gsl, i, j, mat_->A[j][i-j]);
      l1_thresh[i]+=fabs(mat_->A[j][i-j]);
    }
    j = i;
    gsl_matrix_set(mat_gsl, i, j, mat_->A[i][j-i]);
    for ( j = i+1; j < N_in; j++)
    {
      gsl_matrix_set(mat_gsl, i, j, mat_->A[i][j-i]);
      l1_thresh[i]+=fabs(mat_->A[i][j-i]);
    }
    l1_thresh[i] /= (double)(N_in-1);
    l1_thresh[i] *= threshold_;
  }
  gsl_linalg_cholesky_decomp1(mat_gsl);
  for ( i = 0; i < N_in; i++) // going through columns of gsl matrix
  {
    j = i;
    L_->A[i][j-i] = gsl_matrix_get(mat_gsl, j, i);
    for ( j = i+1; j < N_in; j++) L_->A[i][j-i] = ((fabs(mat_->A[i][j-i])) < l1_thresh[i]) ? 0.0 : gsl_matrix_get(mat_gsl, j, i);
  }
  delete l1_thresh;
}

void AY_Choleskyspace::alloc_workspace()
{x_gsl = gsl_vector_alloc(N_in);}

void AY_Choleskyspace::solve_system(AYvec * x_in)
{
  if (x_gsl==NULL) alloc_workspace();
  gsl_linalg_cholesky_decomp1(mat_gsl);
  gsl_linalg_cholesky_svx(mat_gsl, x_gsl);
  x_in->GSL_2_AYvec_copy(x_gsl);
}
void AY_Choleskyspace::solve_system(AYvec * x_in, AYvec * b_in)
{
  if (x_gsl==NULL) alloc_workspace();
  gsl_linalg_cholesky_decomp1(mat_gsl);
  b_in->AYvec_2_GSL_copy(x_gsl);
  gsl_linalg_cholesky_svx(mat_gsl, x_gsl);
  x_in->GSL_2_AYvec_copy(x_gsl);
}
void AY_Choleskyspace::solve_system(AYsym * A_, AYvec * x_in, AYvec * b_in)
{
  if (x_gsl==NULL) alloc_workspace();
  load_mat(A_);
  gsl_linalg_cholesky_decomp1(mat_gsl);
  b_in->AYvec_2_GSL_copy(x_gsl);
  gsl_linalg_cholesky_svx(mat_gsl, x_gsl);
  x_in->GSL_2_AYvec_copy(x_gsl);
}
void AYlinalg_Cholesky_solve(AYsym * L_, AYvec * z_, AYvec *r_) // assumes L_ is already the decomposition
{
  int i, j, N = L_->N;
  for ( i = 0; i < N; i++)
  {
    double sum = r_->A_ptr[i];
    for ( j = 0; j < i; j++) sum -= z_->A_ptr[j]*L_->A[j][i-j];
    z_->A_ptr[i] = sum/(L_->A[i][0]);
  }
  for ( i = N-1; i >= 0; i--)
  {
    double sum = z_->A_ptr[i];
    for ( j = N-1; j > i; j--) sum -= z_->A_ptr[j]*L_->A[i][j-i];
    z_->A_ptr[i] = sum/(L_->A[i][0]);
  }
}

void AYlinalg_svd(AYmat * mat_, AY_SVDspace * space_) // assume long thin matrix for now
{
  mat_->AYmat_2_GSL_copy(space_->U);
  gsl_linalg_SV_decomp(space_->U, space_->V, space_->s, space_->work);
}
