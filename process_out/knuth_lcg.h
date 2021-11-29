#ifndef knuth_lcg_H  /* Include guard */
#define knuth_lcg_H

#include <stdint.h>

double knuth_random_uni(double low, double high, uint64_t * carry);
uint64_t lcg_uni(uint64_t *lcg_carry); // call this evertime you want a random integer, or rand()
uint64_t lcg_fwd(uint64_t seed,uint64_t jump); // eq to srand(seed)
double lcg_sze(); // replaces randmax
double boxmuller_knuth(double mean, double variance, uint64_t * carry);

#endif // knuth_lcg_H
