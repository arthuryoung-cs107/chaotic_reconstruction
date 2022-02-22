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

void AYsym::init_XTWX(AYmat * m_, AYdiag * w_)
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < (N - i); j++)
    {
      A[i][j] = 0.0; // this is being set. Adjust other dims
      for (int k = 0; k < m_->M; k++) // multiply along the full column
      {
        A[i][j] += (w_->A[k])*(m_->AT[i][k])*(m_->AT[j+i][k]);
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
        w->A_ptr[i+j] += A[i][j]*v->A_ptr[i]; // only updates the vector DOWNSTREAM
      }
      out += (w->A_ptr[i])*(v->A_ptr[i]); // actively updating the inner product. w is already fully calculated by this point
    }
  }
  else printf("AYsym: vT_A_v error, dimensions are (%d %d)(%d %d) = (%d %d)\n", N, N, v->M, v->N, w->M, w->N);
  return out;
}

AYdiag::AYdiag(int N_) : N(N_), A(new double[N_])
{}

AYdiag::AYdiag(AYmat * m_) : N((m_->M > m_->N) ? m_->N : m_->M), A(new double[N])
{
  for (int i = 0; i < N; i++) A[i] = m_->AT[i][i];
}

AYdiag::AYdiag(AYvec * v_) : N(v_->M), A(new double[v_->M])
{
  memcpy(A, v_->A_ptr, N*sizeof(double));
}

AYdiag::~AYdiag()
{delete A;}

void AYdiag::print_mat(bool space_)
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < i; j++)
    {
      printf("%f ", 0.0);
    }
    printf("%f ", A[i]);
    for (int j = i+1; j < N; j++)
    {
      printf("%f ", 0.0);
    }
    printf("\n");
  }
  if (space_) printf("\n");
}

void AYdiag::fprintf_diag(char * name, bool verbose) // just interpret as a vector
{
  size_t n_end = strlen(name);
  char * buf = new char[n_end + 10]; strcpy(buf, name); char * buf_it = buf + n_end;

  strcpy(buf_it, ".aydat"); FILE * data_file = fopen(buf, "wb");
  fwrite(A, sizeof(double), N, data_file); fclose(data_file);

  strcpy(buf_it, ".aysml"); FILE * aysml_file = fopen(buf, "w");
  fprintf(aysml_file, "%d 1", N); fclose(aysml_file);

  if (verbose) printf("AYdiag: wrote file  %s.aydat/aysml\n", name);
  delete buf;
}

void AYdiag::init_eye()
{
  for (int i = 0; i < N; i++) A[i] = 1.0;
}

void AYdiag::init_123()
{
  for (int i = 0; i < N; i++) A[i] = (double)(i+1);
}

void AYdiag::init_mat(AYmat * m_)
{
  for (int i = 0; i < N; i++) A[i] = m_->AT[i][i];
}

void AYdiag::init_vec(AYvec * v_)
{
  memcpy(A, v_->A_ptr, N*sizeof(double));
}

void AYdiag::mult_vec(AYvec * in_, AYvec * out_, bool diff_)
{
  if ((N==in_->M)&&(N==out_->M))
  {
    if (diff_) for (int i = 0; i < N; i++) out_->A_ptr[i] -= (A[i])*in_->A_ptr[i];
    else for (int i = 0; i < N; i++) out_->A_ptr[i] += (A[i])*in_->A_ptr[i];
  }
  else
  {
    char msg[200]; sprintf(msg, "AYdiag: mult_vec fatal error, dimension mismatch (AYdiag: %d x %d, in: %d, out: %d)", N, N, in_->M, out_->M); AYfatalerror(msg);
  }
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

void AYfatalerror(char * msg, int exit_code)
{printf("%s\n", msg); exit(exit_code);}
void AYfatalerror(const char *msg, int exit_code){printf("%s\n", msg); exit(exit_code);}
