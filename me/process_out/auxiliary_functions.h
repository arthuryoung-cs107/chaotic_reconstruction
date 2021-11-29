#ifndef auxiliary_functions_H  /* Include guard */
#define auxiliary_functions_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


void aysml_gen(char name[], int m, int n);
void fprintf_matrix(double ** matrix, int M, int N, char prefix[]);
void name_gen(char ptr[], int length, const char * name);

// math functions
double dmax_element(double * vector, int vector_low, int vector_high );
double dmin_element(double * vector, int vector_low, int vector_high ); //assumes index zero
int imax_element(int * ivector, int vector_low, int vector_high );
int imin_element(int * ivector, int vector_low, int vector_high );
void zerom_init(double **a, int arl, int arh, int acl, int ach);
void zeromint_init(int **a, int arl, int arh, int acl, int ach);
void onev_init(double *a, int arl, int arh);
void zerov_init(double *a, int arl, int arh);
void zerovint_init(int *a, int arl, int arh);
void dmatrix_mult(double ** a, int arl, int arh, int acl, int ach, double ** b, int brl, int brh, int bcl, int bch, double ** c );
void dmv_mult(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c);
void dmatrix_add(double ** a, int arl, int arh, int acl, int ach, double **b, int brl, int brh, int bcl, int bch, int sign, double ** result);
void dvector_add(double * a, int arl, int arh, double * b, int sign, double * result);
void dv_scalmult(double * a, int arl, int arh, double scalar, double * result);
void dm_scalmult(double ** a, int arl, int arh, int acl, int ach, double scalar, double ** result);
void mat2vec(double ** G, int nrl, int nrh, int ncl, int nch, double * g, int mrl, int mrh);
void vec2mat(double *g, int nrl, int nrh, double ** G, int mrl, int mrh, int mcl, int mch);
double norm_l2(double * input_vector, int nrl, int nrh);
double norm_frob(double ** X, int nrl, int nrh, int ncl, int nch);
double up_order(double a);
void transpose_dmatrix(double ** a, int arl, int arh, int acl, int ach, double ** b);

void AY_GSLmatrix_add(gsl_matrix * W, gsl_matrix * K, gsl_matrix * work, double scalar);


#endif // auxiliary_functions_H
