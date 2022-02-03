#ifndef AYAUX_H  /* Include guard */
#define AYAUX_H

#include <stdint.h>

void aysml_gen(char name[], int m, int n);
void name_gen(char ptr[], int length, const char * name);

double knuth_random_uni(double low, double high, uint64_t * carry);
uint64_t lcg_uni(uint64_t *lcg_carry); // call this evertime you want a random integer, or rand()
uint64_t lcg_fwd(uint64_t seed,uint64_t jump); // eq to srand(seed)
double lcg_sze(); // replaces randmax
double boxmuller_knuth(double mean, double variance, uint64_t * carry);

int ** AYimatrix(int M_, int N_);
void free_AYimatrix(int ** m_);

double ** AYdmatrix(int M_, int N_);
void free_AYdmatrix(double ** m_);

double *** AYd3tensor(int W_, int M_, int N_);
void free_AYd3tensor(double *** t_);

#endif
