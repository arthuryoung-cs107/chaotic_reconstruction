#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "AYaux.h"

void aysml_gen(char name[], int m, int n)
{
  char specfile[300];
  memset(specfile, 0, 299);
  snprintf(specfile, 300, "%s.aysml", name);
  FILE * aysml_file = fopen(specfile, "w");
  fprintf(aysml_file, "%d %d", m, n);
  fclose(aysml_file);
}

void name_gen(char ptr[], int length, const char * name)
{
  memset(ptr, 0, length-1);
  // snprintf(ptr, length, name);
  memcpy(ptr, name, length*sizeof(char));
}

double knuth_random_uni(double low, double high, uint64_t * carry) // knuth rng
{
  double rand_uni = ((double) lcg_uni(carry))/(lcg_sze());
  rand_uni = (rand_uni - 0.5)*(high - low) + ((high + low)/(2.0));
  return rand_uni;
}

uint64_t lcg_uni(uint64_t *lcg_carry) // call this evertime you want a random integer, or rand()
{
  // **********************************************
  // * KNUTH's 64-BIT LCG RANDOM NUMBER GENERATOR *
  // **********************************************

  #ifdef KNUTH
    //ax+c mod m: m=2^64
    *lcg_carry= (__uint128_t) *lcg_carry*6364136223846793005ULL+1442695040888963407ULL; //knuth
  #else
    *lcg_carry=rand();
  #endif

  return *lcg_carry;
}

uint64_t lcg_fwd(uint64_t seed,uint64_t jump) // eq to srand(seed)
{
  // *********************
  // * LCG JUMP FUNCTION *
  // *********************

  //being lazy here, should use the function  a^k x + (a^k-1)*c/(a-1)

  #ifdef KNUTH
    uint64_t ix;
    uint64_t ret;
    ret=0;
    for(ix=0;ix<jump;ix++)
    {
      ret=lcg_uni(&seed);
    }
    return ret;
  #else
    srand(seed);
    return seed;
  #endif
}

double lcg_sze() // replaces randmax
{
  // *******************************
  // * THE DOMAIN OF RANDOM NUMERS *
  // *******************************/

  #ifdef KNUTH
    return pow(2,64);
  #else
    return (double) RAND_MAX + 1.0;
  #endif
}
double boxmuller_knuth(double mean, double variance, uint64_t * carry)
{
  double sze = lcg_sze();
  double rand_uni1 = ((double) lcg_uni(carry))/sze; //numerator is double, forces double arithmatic
  double rand_uni2 = ((double) lcg_uni(carry))/sze;
  double box1 = -2.0 * log(rand_uni1);
  double gauss = pow(box1, 0.5)* cos(2.0*M_PI*rand_uni2);
  double absg = fabs(gauss);

  if (absg > 4)
  {
    gauss = 4.0 * gauss/absg;
  }

  gauss = (gauss * pow(variance, 0.5)) + mean;
  return gauss;
}

int ** AYimatrix(int M_, int N_)
{
  int * M00 = (int*)malloc((size_t)(N_*M_)*sizeof(int)); // allocating space for all the values of each matrix
  int ** M0 = (int**)malloc((size_t)(M_)*sizeof(int*)); // allocating space for all the pointers to each COLUMN

  for (int i = 0; i < M_; i++) M0[i] = M00 + (i*N_);

  return M0;
}
void free_AYimatrix(int ** in_) {free(in_[0]); free(in_);}


double ** AYdmatrix(int M_, int N_)
{
  double * M00 = (double*)malloc((size_t)(N_*M_)*sizeof(double)); // allocating space for all the values of each matrix
  double ** M0 = (double**)malloc((size_t)(M_)*sizeof(double*)); // allocating space for all the pointers to each COLUMN

  for (int i = 0; i < M_; i++) M0[i] = M00 + (i*N_);

  return M0;
}
void free_AYdmatrix(double ** in_) {free(in_[0]); free(in_);}


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
void free_AYd3tensor(double *** in_) {free(in_[0][0]); free(in_[0]); free(in_);}
