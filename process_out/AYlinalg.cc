#include <cstdio>
#include <sstream>
#include <string>
#include <fstream>

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

AYsym::AYsym(int N_): N(N_)
{
  double Nd = (double) N;
  double Np1 = Nd + 1.0;
  double Nd2 = Nd/2.0;

  int len = (int) (Np1*Nd2); // Baby Gauss theorem

  A = (double**)malloc((size_t)(N)*sizeof(double*));
  *A = (double*)malloc((size_t)(len)*sizeof(double));
  int j = N;
  for (int i = 1; i < N; i++) {A[i] = A[i-1] + j; j--;}
}

AYsym::~AYsym()
{
  for (int i = 0; i < N; i++) free(A[i]);
  free(A);
}

void AYsym::mult_vec(AYvec * in_, AYvec * out_) // specialized for vector multiplication.
{
  if ((in_->M==N)&&(out_->M==N))
  {
    for (int i = 0; i < N; i++)
    {
      out_->A_ptr[i] += A[i][0]*in_->A_ptr[i];
      for (int j = i+1; j < N-i; j++)
      {
        out_->A_ptr[i] += A[i][j]*in_->A_ptr[i];
      }


    }
  }
  else printf("AYsym: mult_set error, dimensions are (%d %d)(%d %d) = (%d %d)\n", N, N, in_->M, in_->N, in_->M, in_->N);
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
  char aysml_specs[300]; memset(aysml_specs, 0, 299); snprintf(aysml_specs, 300, "%s.aysml", name);
  char aydat_name[300]; memset(aydat_name, 0, 299); snprintf(aydat_name, 300, "%s.aydat", name);
  int M, N;
  std::ifstream tens_file;
  tens_file.open(aysml_specs);
  std::string line;
  std::getline(tens_file, line);
  std::istringstream in(line);
  in >> M >> N;
  AYmat * out = new AYmat(M, N);
  FILE * aydat_stream = fopen(aydat_name, "r");
  size_t success = fread(out->A_ptr, sizeof(double), M*N, aydat_stream);
  fclose(aydat_stream);
  return out;
}

AYvec * aysml_read_vec(char name[])
{
  int i, j, M, N;
  double data;
  size_t success;
  char extract[1000];
  char aysml_specfile[300];
  char aydat_specfile[300];
  memset(aysml_specfile, 0, 299);
  memset(aydat_specfile, 0, 299);
  snprintf(aysml_specfile, 300, "%s.aysml", name);
  snprintf(aydat_specfile, 300, "%s.aydat", name);

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

AYtens * aysml_read_tens(char name[])
{
  char aysml_specs[300]; memset(aysml_specs, 0, 299); snprintf(aysml_specs, 300, "%s.aysml", name);
  char aytens_name[300]; memset(aytens_name, 0, 299); snprintf(aytens_name, 300, "%s.aytens", name);
  int type, M, N, W;
  std::ifstream tens_file;
  tens_file.open(aysml_specs);
  std::string line;
  std::getline(tens_file, line);
  std::istringstream in(line);
  in >> type >> M >> N >> W;

  AYtens * T_out = new AYtens(W, M, N);
  if (type == 1)
  {
    FILE * aydat_stream = fopen(aytens_name, "r");
    size_t success = fread(**(T_out->T_AT), sizeof(double), M*N*W, aydat_stream);
    fclose(aydat_stream);
  }
  return T_out;
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

double *** AYd3tensor(int W_, int M_, int N_)
{
  double * M00 = (double*)malloc((size_t)(W_*N_*M_)*sizeof(double)); // allocating space for all the values of each matrix
  double ** M0 = (double**)malloc((size_t)(W_*M_)*sizeof(double*)); // allocating space for all the pointers to each COLUMN
  double *** T = (double***)malloc((size_t)(W_)*sizeof(double**));

  for (int j = 0; j < W_; j++)
  {
    for (int i = 0; i < M_; i++)
    {
      M0[i+(j*M_)] = M00 + (i*N_) + (j*(N_*M_));
    }
    T[j] = M0 + j*(M_);
  }

  return T;
}
void free_AYd3tensor(double *** t_) {free(t_[0][0]); free(t_[0]); free(t_);}
