#include <cstdio>
#include <sstream>
#include <string>
#include <fstream>

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

AYtens::AYtens(char * name)
{
  size_t n_end = strlen(name);
  char * buf = new char[n_end + 50]; strcpy(buf, name); char * buf_it = buf + n_end;

  int num_lines=0;
  int split;
  std::ifstream aysml; strcpy(buf_it, ".aysml"); aysml.open(buf);
  for (std::string line; std::getline(aysml, line);)
  {
    if (num_lines++==0)
    {
      std::istringstream in(line);
      in >> split >> M >> N >> W;
    }
  }
  aysml.close();

  if (num_lines==1)
  {
    T_AT = AYd3tensor(W, N, M);
    mat=(AYmat*)malloc((size_t)(W*sizeof(AYmat)));
    for (int i = 0; i < W; i++)
    {
      mat[i].M = M;
      mat[i].N = N;
      mat[i].AT = T_AT[i];
      mat[i].A_ptr = mat[i].AT[0];
    }
    switch (split)
    {
      case 0: // divided across many files
        for (int i = 0; i < W; i++)
        {
          sprintf(buf_it, ".%d.aydat", i); FILE * data_file = fopen(buf, "r");
          fread_safe(mat[i].A_ptr, sizeof(double), M*N, data_file);
          fclose(data_file);
        }
        break;

      case 1: // sourced from one file
      {
        strcpy(buf_it, ".aydat"); FILE * data_file = fopen(buf, "r");
        fread_safe(T_AT[0][0], sizeof(double), W*M*N, data_file);
        fclose(data_file);
      }
        break;

      default:
        printf("AYtens: read failed. Invalid aysml (split = %d)\n", split);
        exit (EXIT_FAILURE);
    }
  }
  else
  {
    printf("AYtens: read failed. Invalid aysml (lines = %d)\n", num_lines);
    exit (EXIT_FAILURE);
  }
  delete buf;

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
void AYtens::fprintf_tens(char name[], int split_, bool verbose_)
{
  size_t n_end = strlen(name);
  char * buf = new char[n_end + 50]; strcpy(buf, name); char * buf_it = buf + n_end;

  switch(split_)
  {
    case 1: // write as one dat file
    {
      strcpy(buf_it, ".aydat"); FILE * data_file = fopen(buf, "wb");
      fwrite(**T_AT, sizeof(double), W*M*N, data_file); fclose(data_file);
    }
      break;
    case 0: // write as many dat files

      for (int i = 0; i < W; i++)
      {
        sprintf(buf_it, ".%d.aydat", i); FILE * data_file = fopen(buf, "wb");
        fwrite(T_AT[i][0], sizeof(double), M*N, data_file); fclose(data_file);
      }
      break;

    default:
    {
      strcpy(buf_it, ".aydat"); FILE * data_file = fopen(buf, "wb");
      fwrite(**T_AT, sizeof(double), W*M*N, data_file); fclose(data_file);
    }
  }

  strcpy(buf_it, ".aysml"); FILE * aysml_file = fopen(buf, "w");
  fprintf(aysml_file, "%d %d %d %d", split_, M, N, W); fclose(aysml_file);

  if (verbose_) printf("AYtens: wrote file(s)  %s.aydat/aysml\n", name);
  delete buf;
}
void AYtens::init_0()
{for (int i = 0; i < W*M*N; i++) *(**T_AT + i) = 0.0;}

void AYtens::init_123()
{int count = 1; for (int k = 0; k < W; k++) for (int i = 0; i < M; i++) for (int j = 0; j < N; j++, count++) T_AT[k][j][i] = (double) count;}

void AYtens::init_mats123()
{for (int k = 0; k < W; k++) mat[k].init_123();}
