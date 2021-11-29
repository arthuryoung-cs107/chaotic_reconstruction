#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "auxiliary_functions.h"
#include "nrutil.h"

void aysml_gen(char name[], int m, int n)
{
  char specfile[300];
  memset(specfile, 0, 299);
  snprintf(specfile, 300, "%s.aysml", name);
  FILE * aysml_file = fopen(specfile, "w");
  fprintf(aysml_file, "%d %d", m, n);
  fclose(aysml_file);
}

void fprintf_matrix(double ** matrix, int M, int N, char name[])
{
  int i, j;
  char specfile[200];
  memset(specfile, 0, 199);
  snprintf(specfile, 200, "%s.aydat", name);
  FILE * data_file = fopen(specfile, "wb");

  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      fwrite(&((*matrix)[i*N + j]), sizeof(double), 1, data_file); // write by rows
    }
  }
  fclose(data_file);

  aysml_gen( name, M, N);
}

void name_gen(char ptr[], int length, const char * name)
{
  memset(ptr, 0, length-1);
  // snprintf(ptr, length, name);
  memcpy(ptr, name, length*sizeof(char));
}

// math functions
double dmax_element(double * vector, int vector_low, int vector_high )
{
  double initial = 0;
  for (int i = vector_low; i <= vector_high; i++) //assumes zero index
  {
    if (fabs(vector[i]) > initial)
    {
      initial = fabs(vector[i]);
    }
  }
  return initial;
}
double dmin_element(double * vector, int vector_low, int vector_high ) //assumes index zero
{
  double initial = dmax_element(vector, vector_low, vector_high);
  for (int i = vector_low; i <= vector_high; i++) //assumes zero index
  {
    if (fabs(vector[i]) < initial && vector[i] != 0)
    {
      initial = fabs(vector[i]);
    }
  }
  return initial;
}
int imax_element(int * ivector, int vector_low, int vector_high )
{
  double dinitial = 0;
  double vector;
  for (int i = vector_low; i <= vector_high; i++) //assumes zero index
  {
    vector = ivector[i];
    if (fabs(vector) > dinitial)
    {
      dinitial = fabs(vector);
    }
  }
  int initial = dinitial;
  return initial;
}
int imin_element(int * ivector, int vector_low, int vector_high )
{
  int initial = imax_element(ivector, vector_low, vector_high);
  double dinitial = initial;
  double vector;
  for (int i = vector_low; i <= vector_high; i++) //assumes zero index
  {
    vector = ivector[i];
    if (fabs(vector) < dinitial && vector != 0)
    {
      dinitial = fabs(vector);
    }
  }
  int initial2 = dinitial;
  return initial2;
}
void zerom_init(double **a, int arl, int arh, int acl, int ach)
{
  int i, j;
  for ( i = arl; i <= arh; i++)
  {
    for ( j = acl; j <= ach ; j++)
    {
      a[i][j] = 0;
    }
  }
}
void zeromint_init(int **a, int arl, int arh, int acl, int ach)
{
  int i, j;
  for ( i = arl; i <= arh; i++)
  {
    for ( j = acl; j <= ach ; j++)
    {
      a[i][j] = 0;
    }
  }
}
void onev_init(double *a, int arl, int arh)
{
  int i;
  for ( i = arl; i <= arh; i++)
  {
    a[i] = 1;
  }
}
void zerov_init(double *a, int arl, int arh)
{
  int i;
  for ( i = arl; i <= arh; i++)
  {
    a[i] = 0;
  }
}
void zerovint_init(int *a, int arl, int arh)
{
  int i;
  for ( i = arl; i <= arh; i++)
  {
    a[i] = 0;
  }
}
void dmatrix_mult(double ** a, int arl, int arh, int acl, int ach, double ** b, int brl, int brh, int bcl, int bch, double ** c) // a seg fault happens in here when dmv mult is called
{
  int i, j, k;
  double sum;
  int am = arh - arl + 1; //producing dimensions of matrices
  int an = ach - acl + 1;
  int bm = brh - brl + 1;
  int bn = bch - bcl + 1;

  if (an != bm)
  {
      printf(" Error: inner matrix dimensions inequivalent. Matrix multiplication failed. \n");
  }

  for ( i = arl; i <= arh; i++)
  {
    for ( j = bcl; j <= bch; j++) // outer loop for resolving appropriate dimensions of output
    {
      sum = 0;
      for(k = acl; k <= ach; k++)
      {
        sum = sum + a[i][k]*b[k][j];
      }
      c[i][j] = sum;
    }
  }
}
void dmv_mult(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c)
{
  int i, j;
  double sum;

  int an = ach - acl + 1;
  int bm = brh - brl + 1;

  if (an != bm)
  {
    printf(" Error: matrix-vector pair of inequivalent inner-dimensions. Matrix multiplication failed \n");
  }
  else
  {
    for ( i = arl; i <= arh; i++ ) // for each row of output vector
    {
      sum = 0;
      for(j = acl; j <= ach; j++) // each column component of i'th row of a
      {
        sum = sum + a[i][j]*b[j];
      }
      c[i] = sum;
    }
  }
}
void dmatrix_add(double ** a, int arl, int arh, int acl, int ach, double **b, int brl, int brh, int bcl, int bch, int sign, double ** result) //assumes equivalent matrix dimensions, not neccessarily iterating regime
{
  int i, j;
  int am = arh - arl + 1;
  int an = ach - acl + 1;
  int bm = brh - brl + 1;
  int bn = bch - bcl + 1;
  double sign_d = (double) sign;
  if (am != bm || an != bn)
  {
    printf(" Error: matrix of inequivalent dimensions. Matrix addition failed \n");
  }
  for ( i = arl; i <= arh; i++)
  {
    for ( j = acl; j <= ach; j++)
    {
      result[i][j] = a[i][j] + (sign_d * b[i][j]);
    }
  }
}
void dvector_add(double * a, int arl, int arh, double * b, int sign, double * result) //assumes equivalent indexing regime
{
  int i;
  double sign_d = (double)sign;
  for ( i = arl; i <= arh; i++)
  {
    result[i] = a[i] +  sign_d*b[i];
  }
}
void dv_scalmult(double * a, int arl, int arh, double scalar, double * result) // potentially more fitting as a scalar?
{
  for (int i = arl; i <= arh; i++)
  {
    result[i] = a[i]*scalar;
  }
}
void dm_scalmult(double ** a, int arl, int arh, int acl, int ach, double scalar, double ** result)
{
  for (int i = arl; i <= arh; i++)
  {
    for (int j = acl; j <= ach; j++)
    {
      result[i][j] = scalar*a[i][j];
    }
  }
}
void mat2vec(double ** G, int nrl, int nrh, int ncl, int nch, double * g, int mrl, int mrh)
{
  int nrows = nrh - nrl + 1;
  int ncols = nch - ncl + 1;
  int mrows = mrh - mrl + 1;
  if ( mrows != nrows*ncols)
  {
    printf("mat2vec failed: Number of elements in inputted vector inequal to number of elements in inputted matrix \n" );
  }
  if ( nrl != 1 || ncl != 1 || mrl != 1)
  {
    printf("mat2vec failed: ONLY USE MAT2VEC WITH INDEX 1 SCHEMES \n" );
  }
  else
  {
    int vec_iterator = mrl;
    for ( int j = ncl; j <= nch; j++)
    {
      for (int i = nrl; i <= nrh; i++)
      {
        g[vec_iterator] = G[i][j];
        vec_iterator++;
      }
    }
  }
}
void vec2mat(double *g, int nrl, int nrh, double ** G, int mrl, int mrh, int mcl, int mch)
{
  int nrows = nrh - nrl + 1;
  int mrows = mrh - mrl + 1;
  int mcols = mch - mcl + 1;
  if ( nrows != mrows*mcols)
  {
    printf("vec2mat failed: Number of elements in inputted vector inequal to number of elements in inputted matrix \n" );
  }
  else
  {
    int vec_iterator = nrl;
    for (int j = mcl; j <= mch; j++)
    {
      for (int i = mrl; i <= mrh; i++)
      {
        G[i][j] = g[vec_iterator];
        vec_iterator++;
      }
    }
  }
}
double norm_l2 (double * input_vector, int nrl, int nrh) // standard definition, no squaring after square root
{
  double sum = 0;
  double return_sum, local_element, local_element_2;
  for (int i = nrl; i <= nrh; i++)
  {
    local_element = input_vector[i];
    local_element_2 = local_element*local_element;
    sum = sum + local_element_2;
  }
  return_sum = sqrt(sum);
  return return_sum;
}
double norm_frob (double ** X, int nrl, int nrh, int ncl, int nch)
{
  double sum = 0;
  double return_sum, local_element, local_element_2;
  for (int i = nrl; i <= nrh; i++)
  {
    for (int j = ncl; j <= nch; j++)
    {
      local_element = X[i][j];
      local_element_2 = local_element*local_element;
      sum = sum + local_element_2;
    }
  }
  return_sum = sqrt(sum);
  return return_sum;
}

double up_order(double a)
{
  double result = log(a)/log(10);
  result = ceil(result);
  return pow(10, result);
}

void AY_GSLmatrix_add(gsl_matrix * W, gsl_matrix * K, gsl_matrix * work, double scalar)
{
  gsl_matrix_memcpy(work, K);
  gsl_matrix_scale(work, scalar);
  gsl_matrix_add(W, work);
}
