#include <cstdio>

#include "AYlinalg.hh"

extern "C"
{
  #include "nrutil.h"
  #include "auxiliary_functions.h"
  #include "knuth_lcg.h"
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

void AYlinalg_svd(AYmat * mat_, AY_SVDspace * space_) // assume long thin matrix for now
{
  mat_->AYmat_2_GSL_copy(space_->U);
  gsl_linalg_SV_decomp(space_->U, space_->V, space_->s, space_->work);
}


AYvec * AYmat_2_AYvec_gen(AYmat * X_in)
{
  AYvec * x_out = new AYvec((X_in->M)*(X_in->N));
  memcpy(x_out->A_ptr, X_in->A_ptr, (X_in->M)*(X_in->N)*sizeof(double));
  return x_out;
}

void AYmat_2_AYvec_copy(AYmat * X_in, AYvec * x_in)
{
  if ((X_in->M)*(X_in->N) == x_in->M) memcpy(x_in->A_ptr, X_in->A_ptr, (X_in->M)*(X_in->N)*sizeof(double));
  else printf("AYmat: AYmat_2_AYvec_copy failed, inequal dimensions\n");
}

AYmat * AYvec_2_AYmat_gen(AYvec * x_in)
{
  AYmat * X_out = new AYmat((x_in->M), 1);
  memcpy(X_out->A_ptr, x_in->A_ptr, (x_in->M)*sizeof(double));
  return X_out;
}
void AYvec_2_AYmat_copy(AYvec * x_in, AYmat * X_in)
{
  if ((x_in->M) == (X_in->M)*(X_in->N)) memcpy(X_in->A_ptr, x_in->A_ptr, (X_in->M)*sizeof(double));
  else printf("AYvec: AYvec_2_AYmat_copy failed, inequal dimensions\n");
}

AYmat * aysml_read(char name[])
{
  int i, j, M, N;
  double data;
  size_t success;
  char extract[1000];
  char aysml_specfile[200]; memset(aysml_specfile, 0, 199); snprintf(aysml_specfile, 200, "%s.aysml", name);
  char aydat_specfile[200]; memset(aydat_specfile, 0, 199); snprintf(aydat_specfile, 200, "%s.aydat", name);

  FILE * aysml_stream = fopen(aysml_specfile, "r");

  char * check = fgets(extract, sizeof(extract), aysml_stream);
  M = atoi(extract);
  int next = (int) log10((double) M) + 1;
  N = atoi(extract+next);

  AYmat * out = new AYmat(M, N);
  fclose(aysml_stream);

  FILE * aydat_stream = fopen(aydat_specfile, "r");
  success = fread(out->A_ptr, sizeof(double), M*N, aydat_stream );
  fclose(aydat_stream);

  return out;
}

AYvec * aysml_read_vec(char name[])
{
  int i, j, M, N;
  double data;
  size_t success;
  char extract[1000];
  char aysml_specfile[200];
  char aydat_specfile[200];
  memset(aysml_specfile, 0, 199);
  memset(aydat_specfile, 0, 199);
  snprintf(aysml_specfile, 200, "%s.aysml", name);
  snprintf(aydat_specfile, 200, "%s.aydat", name);

  FILE * aysml_stream = fopen(aysml_specfile, "r");

  char * check = fgets(extract, sizeof(extract), aysml_stream);
  M = atoi(extract);
  int next = (int) log10((double) M) + 1;
  N = atoi(extract+next);

  if (N != 1) printf("aysml_read_vec warning: aysml has N = %d. Reading first M = %d values\n", N, M);

  AYvec * out = new AYvec(M);
  fclose(aysml_stream);

  FILE * aydat_stream = fopen(aydat_specfile, "r");

  for ( i = 0; i < M; i++)
  {
    success = fread(&data, sizeof(double), 1, aydat_stream );
    out->set(i, data);
  }
  fclose(aydat_stream);

  return out;
}

AYmat * GSL_2_AYmat_gen(gsl_matrix * mat_in)
{
  AYmat * mat_out = new AYmat(mat_in->size1, mat_in->size2);
  for (int i = 0; i < mat_out->M; i++)
  {
    for (int j = 0; j < mat_out->N; j++) mat_out->set(i, j, gsl_matrix_get(mat_in, i, j));
  }
  return mat_out;
}

AYmat * GSL_2_AYmat_gen(gsl_vector * vec_in)
{
  AYmat * vec_out = new AYmat(vec_in->size, 1);
  for (int i = 0; i < vec_out->M; i++) vec_out->set(i, 0, gsl_vector_get(vec_in, i));
  return vec_out;
}

AYmat * GSL_2_diagAYmat_gen(gsl_vector * vec_in)
{
  AYmat * mat_out = new AYmat(vec_in->size, vec_in->size);
  mat_out->init_0();
  for (int i = 0; i < mat_out->M; i++) mat_out->set(i, i, gsl_vector_get(vec_in, i));
  return mat_out;
}
AYvec * GSL_2_AYvec_gen(gsl_matrix * mat_in)
{
  AYvec * vec_out = new AYvec((mat_in->size1)*(mat_in->size2));
  for (int i = 0; i < mat_in->size1; i++) {for (int j = 0; j < mat_in->size2; j++) vec_out->set((i*(mat_in->size2))+j, gsl_matrix_get(mat_in, i, j));}
  return vec_out;
}
AYvec * GSL_2_AYvec_gen(gsl_vector * vec_in)
{
  AYvec * vec_out = new AYvec(vec_in->size);
  for (int i = 0; i < vec_out->M; i++) vec_out->set(i, gsl_vector_get(vec_in, i));
  return vec_out;
}
