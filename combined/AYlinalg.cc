#include <cstdio>
#include <sstream>
#include <string>
#include <fstream>

#include "AYlinalg.hh"

extern "C"
{
  #include "AYaux.h"
}

AYsym::AYsym(int N_): N(N_)
{
  double Nd = (double) N;
  double Np1 = Nd + 1.0;
  double Nd2 = Nd/2.0;

  len = (int) (Np1*Nd2); // Baby Gauss theorem

  double *row00 = (double*)malloc((size_t)(len)*sizeof(double));
  A =(double**)malloc((size_t)(N)*sizeof(double*));
  int k = N;
  A[0] = row00;
  for (int i = 1; i < N; i++, k--) A[i] = A[i-1] + k;
}

AYsym::~AYsym()
{
  free(A[0]);
  free(A);
}


void AYsym::init_eye()
{
  for (int i = 0; i < N; i++)
  {
    A[i][0] = 1.0;
    for (int j = 1; j < (N-i); j++) A[i][j] = 0.0;
  }
}

void AYsym::print_mat(bool space_)
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < i; j++)
    {
      printf("%f ", A[j][i-j]);
    }
    for (int j = i; j < N; j++)
    {
      printf("%f ", A[i][j-i]);
    }
    printf("\n");
  }
  if (space_) printf("\n");
}

void AYsym::fprintf_sym(char name[], bool verbose_)
{
  size_t n_end = strlen(name);
  char * buf = new char[n_end + 10]; strcpy(buf, name); char * buf_it = buf + n_end;

  strcpy(buf_it, ".aydat"); FILE * data_file = fopen(buf, "wb");
  fwrite(*A, sizeof(double), len, data_file); fclose(data_file);

  strcpy(buf_it, ".aysml"); FILE * aysml_file = fopen(buf, "w");
  fprintf(aysml_file, "1 %d %d", N, len); fclose(aysml_file);

  if (verbose_) printf("AYsym: wrote file  %s.aydat/aysml\n", name);
  delete buf;
}

void AYsym::init_123()
{for (int i = 0; i < len; i++) A[0][i] = (double) (i+1);}

void AYsym::init_sqrmat(AYmat * m_)
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < (N - i); j++)
    {
      A[i][j] = 0.0; // this is being set. Adjust other dims
      for (int k = 0; k < m_->M; k++)
      {
        A[i][j] += (m_->AT[i][k])*(m_->AT[j+i][k]);
      }
    }
  }
}

void AYsym::mult_vec(AYvec * in_, AYvec * out_, bool diff_) // specialized for vector multiplication. Adds on top of output vector to keep it quick and simple
{
  if ((in_->M==N)&&(out_->M==N))
  {
    if (diff_)
    {
      for (int i = 0; i < N; i++)
      {
        out_->A_ptr[i] -= A[i][0]*in_->A_ptr[i]; // diagonal element
        for (int j = 1; j < (N - i); j++)
        {
          out_->A_ptr[i] -= A[i][j]*in_->A_ptr[i+j];
          out_->A_ptr[i+j] -= A[i][j]*in_->A_ptr[i];
        }
      }
    }
    else
    {
      for (int i = 0; i < N; i++)
      {
        out_->A_ptr[i] += A[i][0]*in_->A_ptr[i]; // diagonal element
        for (int j = 1; j < (N - i); j++)
        {
          out_->A_ptr[i] += A[i][j]*in_->A_ptr[i+j];
          out_->A_ptr[i+j] += A[i][j]*in_->A_ptr[i];
        }
      }
    }
  }
  else printf("AYsym: mult_set error, dimensions are (%d %d)(%d %d) = (%d %d)\n", N, N, in_->M, in_->N, out_->M, out_->N);
}

double AYsym::vT_A_v(AYvec * v, AYvec * w)
{
  double out=0.0;
  if ((v->M==N)&&(w->M==N))
  {
    for (int i = 0; i < N; i++)
    {
      w->A_ptr[i] += A[i][0]*v->A_ptr[i]; // diagonal element
      for (int j = 1; j < (N - i); j++)
      {
        w->A_ptr[i] += A[i][j]*v->A_ptr[i+j];
        w->A_ptr[i+j] += A[i][j]*v->A_ptr[i];
      }
      out += (w->A_ptr[i])*(v->A_ptr[i]); // actively updating the inner product
    }
  }
  else printf("AYsym: vT_A_v error, dimensions are (%d %d)(%d %d) = (%d %d)\n", N, N, v->M, v->N, w->M, w->N);
  return out;
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

char * string_gen_pruned(char * in_)
{
  size_t len = (size_t)(strlen(in_) + 1);
  char * out = new char[len];
  strcpy(out, in_);
  return out;
}

char * string_gen_pruned(const char * in_)
{
  size_t len = (size_t)(strlen(in_) + 1);
  char * out = new char[len];
  strcpy(out, in_);
  return out;
}
