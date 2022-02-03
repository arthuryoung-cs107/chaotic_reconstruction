#include <sys/stat.h>
#include <cstdio>

#include "AYlinalg.hh"

extern "C"
{
  #include "AYaux.h"
}

AYtens::AYtens(int W_, int M_, int N_): W(W_), M(M_), N(N_), T_AT(AYd3tensor(W_, N_, M_)), mat((AYmat*)malloc((size_t)(W_*sizeof(AYmat))))
{
  for (int i = 0; i < W; i++)
  {
    mat[i].M = M;
    mat[i].N = N;
    mat[i].AT = T_AT[i];
    mat[i].A_ptr = mat[i].AT[0];
  }
}
AYtens::~AYtens()
{
  free_AYd3tensor(T_AT);
  free(mat);
}

double AYtens::get(int i, int j, int k) // the more intuitive way of indexing, I guess
{return T_AT[k][j][i];}
void AYtens::set(int i, int j, int k, double val)
{T_AT[k][j][i] = val;}

void AYtens::print_dims()
{printf("AYtens: W = %d long set of matrices of dimension (M, N) = (%d, %d)\n", W, M, N);}

void AYtens::print_tens(bool space_)
{
  int i, j, k;

  for (k = 0; k < W; k++)
  {
    for ( i = 0; i < M; i++)
    {
      for ( j = 0; j < N; j++) printf("%f ", T_AT[k][j][i]);
      printf("\n");
    }
    printf("\n");
  }
  if (space_) printf("\n");
}
void AYtens::fprintf_tens(char name[], bool verbose_)
{
  char specfile[300]; memset(specfile, 0, 299); snprintf(specfile, 300, "%s.aytens", name);
  char smlfile[300]; memset(smlfile, 0, 299); snprintf(smlfile, 300, "%s.aysml", name);
  FILE * data_file = fopen(specfile, "wb");
  FILE * aysml_file = fopen(smlfile, "w");
  fprintf(aysml_file, "1 %d %d %d", M, N, W);
  fclose(aysml_file);
  fwrite(**T_AT, sizeof(double), W*M*N, data_file);
  fclose(data_file);
  if (verbose_) printf("wrote file: %s.aytens/aysml\n", name);
}
void AYtens::init_0()
{for (int i = 0; i < W*M*N; i++) *(**T_AT + i) = 0.0;}

void AYtens::init_123()
{int count = 1; for (int k = 0; k < W; k++) for (int i = 0; i < M; i++) for (int j = 0; j < N; j++, count++) T_AT[k][j][i] = (double) count;}

void AYtens::init_mats123()
{for (int k = 0; k < W; k++) mat[k].init_123();}
