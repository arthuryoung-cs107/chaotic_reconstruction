#include <cstdio>
#include <sstream>
#include <string>
#include <fstream>

#include "AYlinalg.hh"
#include "AYblas.h"

extern "C"
{
  #include "AYaux.h"
}

// note: this is designed to play nicely with BLAS, hence we store column major
AYmat::AYmat(int M_, int N_): M(M_), N(N_), AT(AYdmatrix(N, M)), dmatrix_alloc_flag(true)
{
  A_ptr = AT[0];
}

AYmat::AYmat(char * name)
{
  size_t n_end = strlen(name);
  char * buf = new char[n_end + 10]; strcpy(buf, name); char * buf_it = buf + n_end;

  int num_lines=0;

  std::ifstream aysml; strcpy(buf_it, ".aysml"); aysml.open(buf);
  for (std::string line; std::getline(aysml, line);)
  {
    if (num_lines++==0)
    {
      std::istringstream in(line);
      in >> M >> N;
    }
  }
  aysml.close();

  if (num_lines==1)
  {
    if (M>1)
    {
      AT = AYdmatrix(N, M); dmatrix_alloc_flag = true; A_ptr = AT[0];

      strcpy(buf_it, ".aydat"); FILE * data_file = fopen(buf, "r");
      fread_safe(A_ptr, sizeof(double), M*N, data_file);
      fclose(data_file);
    }
    else
    {
      printf("AYmat: read failed. Invalid aysml (M = %d). Did you mean AYsym?\n", M);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    printf("AYmat: read failed. Invalid aysml (lines = %d)\n", num_lines);
    exit (EXIT_FAILURE);
  }
  delete buf;
}

AYmat::~AYmat()
{
  if  (dmatrix_alloc_flag) free_AYdmatrix(AT);
  if (GSL_flag != 0)
  {
    if (GSL_flag == 1)
    {
      gsl_vector_free(v_gsl);
    }
    else if (GSL_flag == 2)
    {
      gsl_matrix_free(A_gsl);
    }
  }
}

void AYmat::set(int i, int j, double val)
{
  if ( (i >= 0)&&(i<M)&&(j>=0)&&(j<N) )
  {
    AT[j][i] = val;
  }
  else
  {
    printf("AYmat: set failed. Setting (%d, %d) of %d by %d matrix \n", i, j, M, N);
  }
}

double AYmat::get(int i, int j)
{
  if ( (i >= 0)&&(i<M)&&(j>=0)&&(j<N) )
  {
    return AT[j][i];
  }
  else
  {
    printf("AYmat: get failed. Getting (%d, %d) of %d by %d matrix \n", i, j, M, N);
    return 0;
  }
}

void AYmat::GSL_init()
{
  int i, j;
  if (GSL_flag == 0)
  {
    if (N == 1)
    { v_gsl = gsl_vector_alloc(M); GSL_flag = 1; }
    else
    { A_gsl = gsl_matrix_alloc(M, N); GSL_flag = 2; }
  }
  if (N == 1) { for ( i = 0; i < M; i++) gsl_vector_set(v_gsl, i, AT[0][i]); }
  else
  {
    for ( i = 0; i < M; i++)
    {
      for ( j = 0; j < N; j++)
      {
        gsl_matrix_set(A_gsl, i, j, AT[j][i]);
      }
    }
  }
}

void AYmat::GSL_send()
{
  int i, j;
  if (GSL_flag == 0) printf("GSL_send failed: GSL matrix not allocated.\n");
  else if (GSL_flag == 1) for ( i = 0; i < M; i++) AT[0][i] = gsl_vector_get(v_gsl, i);
  else if (GSL_flag == 2)
  {
    for ( i = 0; i < M; i++)
    {
      for ( j = 0; j < N; j++)
      {
        AT[j][i] = gsl_matrix_get(A_gsl, i, j);
      }
    }
  }
}
void AYmat::GSL_2_AYmat_copy(gsl_matrix * mat_in)
{
  if ((mat_in->size1 == M) && (mat_in->size2 == N))
  {
    for (int i = 0; i < M; i++)
    {
      for (int j = 0; j < N; j++) AT[j][i] = gsl_matrix_get(mat_in, i, j);
    }
  }
  else printf("GSL_2_AYmat_copy (gslmat) failed: dimension mismatch\n");
}
void AYmat::GSL_2_AYmat_copy(gsl_vector * vec_in)
{
  if ((vec_in->size == M) && (N == 1))
  {
    for (int i = 0; i < M; i++) A_ptr[i] = gsl_vector_get(vec_in, i);
  }
  else printf("GSL_2_AYmat_copy (gslvec) failed: dimension mismatch\n");
}
void AYmat::AYmat_2_GSL_copy(gsl_matrix * mat_in)
{
  if ((mat_in->size1 == M) && (mat_in->size2 == N))
  {
    for (int i = 0; i < M; i++)
    {
      for (int j = 0; j < N; j++) gsl_matrix_set(mat_in, i, j, AT[j][i]);
    }
  }
  else printf("AYmat_2_GSL_copy (gslmat) failed: dimension mismatch\n");
}
void AYmat::AYmat_2_GSL_copy(gsl_vector * vec_in)
{
  if ((vec_in->size == M) && (N == 1))
  {
    for (int i = 0; i < M; i++) gsl_vector_set(vec_in, i, A_ptr[i]);
  }
  else printf("AYmat_2_GSL_copy (gslvec) failed: dimension mismatch\n");
}

void AYmat::print_mat(bool space_)
{
  int i, j;
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      printf("%f ", AT[j][i]);
    }
    printf("\n");
  }
  if (space_) printf("\n");
}

void AYmat::print_dims()
{printf("AYmat: matrix has dimensions (M, N) = (%d, %d)\n", M, N);}

void AYmat::fprintf_mat(char name[], bool verbose)
{
  size_t n_end = strlen(name);
  char * buf = new char[n_end + 10]; strcpy(buf, name); char * buf_it = buf + n_end;

  strcpy(buf_it, ".aydat"); FILE * data_file = fopen(buf, "wb");
  fwrite(A_ptr, sizeof(double), M*N, data_file); fclose(data_file);

  strcpy(buf_it, ".aysml"); FILE * aysml_file = fopen(buf, "w");
  fprintf(aysml_file, "%d %d", M, N); fclose(aysml_file);

  if (verbose) printf("AYmat: wrote file  %s.aydat/aysml\n", name);
  delete buf;
}

void AYmat::init_123()
{
  for (int i = 0; i < M; i++) for (int j = 0; j < N; j++) AT[j][i] = (double)((i*N + j)+1);
}

void AYmat::init_0()
{for (int i = 0; i < N*M; i++) A_ptr[i] = 0.0;}

void AYmat::init_randuni(double low_, double high_)
{
  AYuniform gen(low_, high_);
  for (int i = 0; i < N*M; i++) A_ptr[i] = gen.rand_gen();
}

void AYmat::init_randgen(AYrng * rng_in)
{for ( int i = 0; i < N*M; i++) AT[0][i] = rng_in->rand_gen();}

AYmat * AYmat::copy_gen()
{
  AYmat * X_out = new AYmat(M, N);
  memcpy(X_out->A_ptr, A_ptr, M*N*sizeof(double));
  return X_out;
}

AYmat * AYmat::row_slice_gen(int i_)
{
  AYmat * v_out = new AYmat(N, 1);
  for (int j = 0; j < N; j++)
  {
    v_out->A_ptr[j] = AT[j][i_];
  }
  return v_out;
}

AYmat * AYmat::col_slice_gen(int j_)
{
  AYmat * v_out = new AYmat(M, 1);
  for (int i = 0; i < M; i++) v_out->A_ptr[i] = AT[j_][i];
  return v_out;
}

void AYmat::copy_set(AYmat * X_in)
{
  if ((X_in->M == M)&&(X_in->N == N) ) memcpy(X_in->A_ptr, A_ptr, M*N*sizeof(double));
  else printf("AYmat: copy_set failed, inequal dimensions\n");
}
void AYmat::copy_set(AYmat * X_in, double scal_)
{
  if ((X_in->M == M)&&(X_in->N == N) ) for (int i = 0; i < M*N; i++) X_in->A_ptr[i] = scal_*(A_ptr[i]);
  else printf("AYmat: copy_set failed, inequal dimensions\n");
}
void AYmat::copy_set_col(AYmat * mat_out, int local_index, int out_index)
{
  if (M == mat_out->M) memcpy(mat_out->AT[out_index], AT[local_index], M*sizeof(double));
  else printf("AYmat: copy_set_col failed, inequal first dimension (%d vs. %d)\n", M, mat_out->M);
}

AYmat * AYmat::transpose_gen()
{
  int i, j;
  AYmat * X_out = new AYmat(N, M);
  for ( i = 0; i < M; i++) {for ( j = 0; j < N; j++) { X_out->A_ptr[i*N + j] = AT[j][i];}}
  return X_out;
}

void AYmat::transpose_set(AYmat * X_in)
{
  if ((X_in->M == N)&&(X_in->N == M) )
  {
    int i, j;
    for ( i = 0; i < M; i++) { for ( j = 0; j < N; j++) { X_in->A_ptr[i*N + j] = AT[j][i];}}
  }
  else printf("AYmat: transpose_set failed, bad dimensions\n");
}

void AYmat::add(AYmat * B_in, AYmat * C_in, double alpha, double beta )
{
  if ( (B_in->M==M)&&(B_in->N==N)&&(C_in->M==M)&&(C_in->N==N)  )
  {
    for (int i = 0; i < M*N; i++) C_in->A_ptr[i] = A_ptr[i] + alpha*B_in->A_ptr[i] + beta*C_in->A_ptr[i];
  }
  else printf("AYmat: add failed, inequal dimensions\n");
}
void AYmat::add(double alpha_, AYmat * out_)
{
  if ((out_->M == M)&&(out_->M == M))
  {
    for (int i = 0; i < M*N; i++) out_->A_ptr[i] = A_ptr[i] + alpha_;
  }
  else printf("AYmat: add failed, inequal dimensions\n");
}

void AYmat::scal_mult(double c)
{for (int i = 0; i < M*N; i++) A_ptr[i] *= c;}

// take in one matrix, post multiply this matrix with input, return pointer to new matrix
AYmat * AYmat::mult_gen(AYmat * B_in, double alpha, bool trans1_, bool trans2_)
{
  int dim1_A = M; int dim2_A = N;
  int dim1_B = B_in->M; int dim2_B = B_in->N;
  int LDA = M; int LDB = B_in->M;
  char trans1 = 'n';
  char trans2 = 'n';
  if (trans1_)
  {
    dim1_A = N; dim2_A = M;
    trans1 = 't';
  }
  if (trans2_)
  {
    dim1_B = B_in->N; dim2_B = B_in->M;
    trans2 = 't';
  }
  if (dim2_A == dim1_B)
  {
    double beta = 0.0;
    int inc = 1;
    AYmat * C_out = new AYmat(dim1_A, dim2_B);
    if (B_in->N == 1) // matrix vector multiplication
    {
      dgemv_(&trans1, &(M), &(N), &alpha, (A_ptr), &(M), (B_in->A_ptr), &inc, &beta, (C_out->A_ptr), &inc);
    }
    else // matrix-matrix multiplication
    {
      dgemm_(&trans1,&trans2,&(dim1_A),&(dim2_B),&(dim2_A),&alpha,(A_ptr),&(LDA), (B_in->A_ptr),&(LDB),&beta,(C_out->A_ptr),&(C_out->M));
    }
    return C_out;
  }
  else
  {
    printf("AYmat: mult_gen failed, bad dimensions\n");
    return new AYmat(M, B_in->N);
  }

}

// take in two matrices, post multiply this one with B, set contents of C_in
void AYmat::mult_set(AYmat * B_in, AYmat * C_in, double alpha, double beta, bool trans1_, bool trans2_)
{
  int dim1_A = M; int dim2_A = N;
  int dim1_B = B_in->M; int dim2_B = B_in->N;
  int LDA = M; int LDB = B_in->M;
  char trans1 = 'n';
  char trans2 = 'n';

  if (trans1_)
  {
    dim1_A = N; dim2_A = M;
    trans1 = 't';
  }
  if (trans2_)
  {
    dim1_B = B_in->N; dim2_B = B_in->M;
    trans2 = 't';
  }
  if ((dim1_B==dim2_A)&&(C_in->M==dim1_A)&&(C_in->N==dim2_B))
  {
    int inc = 1;
    if (B_in->N == 1) // matrix vector multiplication
    {
      dgemv_(&trans1, &(M), &(N), &alpha, (A_ptr), &(M), (B_in->A_ptr), &inc, &beta, C_in->A_ptr, &inc);
    }
    else // matrix-matrix multiplication
    {
      dgemm_(&trans1,&trans2,&(dim1_A),&(dim2_B),&(dim2_A),&alpha,(A_ptr),&(LDA), (B_in->A_ptr),&(LDB),&beta,(C_in->A_ptr),&(C_in->M));
    }
  }
  else
  {
    printf("AYmat: mult_set failed, bad dimension. Attempting (%d x %d) = (%d x %d).(%d x %d)\n", C_in->M, C_in->N, dim1_A, dim2_A, dim1_B, dim2_B );
  }
}

double AYmat::mean()
{
  double sum = 0.0;
  for (int j = 0; j < N; j++) for (int i = 0; i < M; i++) sum += AT[j][i];
  return sum/((double) M*N);
}

double AYmat::variance()
{
  double diff;
  double sum = 0.0;
  double mu = mean();

  for (int j = 0; j < N; j++)
  {
    for (int i = 0; i < M; i++)
    {
      diff = AT[j][i] - mu;
      sum += diff*diff;
    }
  }
  return sum/((double) M*N);
}

double AYmat::inner(AYmat * B_in)
{
  if ((B_in->M==M)&&(B_in->N==N))
  {
    int i, j;
    double sum = 0;
    for ( j = 0; j < N; j++) for ( i = 0; i < M; i++) sum += (AT[j][i] * B_in->AT[j][i]);
    return sum;
  }
  else
  {
    printf("AYmat: inner failed, inequal dimensions\n");
    return 0;
  }
}

double AYmat::norm_frob()
{
  double out = 0.0;
  for (int j = 0; j < N; j++) for (int i = 0; i < M; i++) out += (AT[j][i])*(AT[j][i]);
  return sqrt(out);
}

double AYmat::norm_frob(AYmat * X_)
{
  double out = 0.0;
  for (int j = 0; j < N; j++) for (int i = 0; i < M; i++) out += (AT[j][i]-X_->AT[j][i])*(AT[j][i]-X_->AT[j][i]);
  return sqrt(out);
}

double AYmat::norm_1()
{
  int i, j;
  double out = 0.0;

  for ( j = 0; j < N; j++) for ( i = 0; i < M; i++) out += fabs(AT[j][i]);

  return out;
}

double AYmat::max_val_mag()
{
  int i;
  double out = 0.0;
  for ( i = 0; i < M*N; i++) { if (fabs(A_ptr[i]) > out) out = fabs(A_ptr[i]); }
  return out;
}

double AYmat::max_val_mag(int * index, int start_index)
{
  if (start_index<(M*N-1))
  {
    int i;
    double out = fabs(A_ptr[start_index]); * index = start_index;
    for ( i = start_index+1; i < M*N; i++)
    {
      if (fabs(A_ptr[i]) > out)
      {
        * index = i;
        out = fabs(A_ptr[i]);
      }
    }
    return out;
  }
  else
  {
    * index = M*N-1; // just return last value
    return fabs(A_ptr[M*N-1]);
  }
}
double AYmat::min_val_mag()
{
  int i;
  double out = fabs(A_ptr[0]);
  for ( i = 1; i < M*N; i++) { if (fabs(A_ptr[i])<out) out = fabs(A_ptr[i]); }
  return out;
}

double AYmat::min_val_mag(int * index, int start_index)
{
  if (start_index<(M*N-1))
  {
    int i;
    double out = fabs(A_ptr[start_index]); * index = start_index;
    for ( i = start_index+1; i < M*N; i++)
    {
      if (fabs(A_ptr[i])<out)
      {
        * index = i;
        out = fabs(A_ptr[i]);
      }
    }
    return out;
  }
  else
  {
    * index = M*N-1; // just return the last element
    return fabs(A_ptr[M*N-1]);
  }
}

int AYmat::norm_0(double tol)
{
  int i, j;
  int out = 0;
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      if (fabs(AT[j][i]) > tol) out++;
    }
  }
  return out;
}

int AYmat::max_mag_elements(AYmat *top_vec_, int *index_array_)
{
  if (top_vec_->M <= M*N) // ensure we are only looking for a subset of magnitudes
  {
    int i, i_it;
    double min_large;

    for ( i = 0; i < top_vec_->M; i++)
    {
      top_vec_->A_ptr[i] = A_ptr[i];  // initializing the magnitudes to be the first few
      index_array_[i] = i; // setting first few indices
    }
    min_large = top_vec_->min_val_mag(&i_it); // determine the smallest magnitude on this list
    for ( i = top_vec_->M; i < N*M; i++) // go through the rest of them
    {
      if (fabs(A_ptr[i]) > min_large) // if we find somethng bigger than the current smallest element on the top list
      {
        top_vec_->A_ptr[i_it] = A_ptr[i]; // kick off the smallest 'large' element
        index_array_[i_it] = i;
        min_large = top_vec_->min_val_mag(&i_it); // determine the new smallest 'large' element
      }
    }
    return i_it;
  }
  else
  {
    printf("AYmat: max_mag_elements failed, searching for %d values within M.N = %d elements\n", top_vec_->N, M*N);
    return -1;
  }
}

int AYmat::min_mag_elements(AYmat *top_vec_, int *index_array_)
{
  if (top_vec_->M <= M*N) // ensure we are only looking for a subset of magnitudes
  {
    int i, i_it;
    double max_small;

    for ( i = 0; i < top_vec_->M; i++)
    {
      top_vec_->A_ptr[i] = A_ptr[i];  // initializing the magnitudes to be the first few
      index_array_[i] = i; // setting first few indices
    }
    max_small = top_vec_->max_val_mag(&i_it); // determine the largest magnitude on this list
    for ( i = top_vec_->M; i < N*M; i++) // go through the rest of them
    {
      if (fabs(A_ptr[i]) < max_small) // if we find somethng smaller than the current largest element on the top list
      {
        top_vec_->A_ptr[i_it] = A_ptr[i]; // kick off the largest 'small' element
        index_array_[i_it] = i;
        max_small = top_vec_->max_val_mag(&i_it); // determine the new largest 'small' element
      }
    }
    return i_it;
  }
  else
  {
    printf("AYmat: min_mag_elements failed, searching for %d values within M.N = %d elements\n", top_vec_->N, M*N);
    return -1;
  }
}

void AYmat::max_mag_elements_ordered(AYmat *top_vec_, int *index_array_)
{
  if (top_vec_->M <= M*N) // ensure we are only looking for a subset of magnitudes
  {
    int i_it, index_temp;
    double val, val_temp;
    i_it = max_mag_elements(top_vec_, index_array_);
    val = top_vec_->max_val_mag(&i_it);
    val_temp = top_vec_->A_ptr[0]; index_temp = index_array_[0];
    top_vec_->A_ptr[0] = top_vec_->A_ptr[i_it]; index_array_[0] = index_array_[i_it];
    top_vec_->A_ptr[i_it] = val_temp; index_array_[i_it] = index_temp;
    max_mag_elements_recursive(top_vec_, index_array_, 1);
  }
  else printf("AYmat: max_mag_elements_ordered failed, searching for %d values within M.N = %d elements\n", top_vec_->N, M*N);
}
void AYmat::min_mag_elements_ordered(AYmat *top_vec_, int *index_array_)
{
  if (top_vec_->M <= M*N) // ensure we are only looking for a subset of magnitudes
  {
    int i_it, index_temp;
    double val, val_temp;
    i_it = min_mag_elements(top_vec_, index_array_);
    val = top_vec_->min_val_mag(&i_it);
    val_temp = top_vec_->A_ptr[0]; index_temp = index_array_[0];
    top_vec_->A_ptr[0] = top_vec_->A_ptr[i_it]; index_array_[0] = index_array_[i_it];
    top_vec_->A_ptr[i_it] = val_temp; index_array_[i_it] = index_temp;
    min_mag_elements_recursive(top_vec_, index_array_, 1);
  }
  else printf("AYmat: min_mag_elements_ordered failed, searching for %d values within M.N = %d elements\n", top_vec_->N, M*N);
}
void AYmat::max_mag_elements_recursive(AYmat *top_vec_, int * index_array_, int i_next)
{
  int i_it, index_temp;
  double val, val_temp;
  val = top_vec_->max_val_mag(&i_it, i_next);
  val_temp = top_vec_->A_ptr[i_next]; index_temp = index_array_[i_next];
  top_vec_->A_ptr[i_next] = top_vec_->A_ptr[i_it]; index_array_[i_next] = index_array_[i_it];
  top_vec_->A_ptr[i_it] = val_temp; index_array_[i_it] = index_temp;
  if (i_next < (top_vec_->M-1)) max_mag_elements_recursive(top_vec_, index_array_, i_next+1);
}

void AYmat::min_mag_elements_recursive(AYmat *top_vec_, int * index_array_, int i_next)
{
  int i_it, index_temp;
  double val, val_temp;
  val = top_vec_->min_val_mag(&i_it, i_next);
  val_temp = top_vec_->A_ptr[i_next]; index_temp = index_array_[i_next];
  top_vec_->A_ptr[i_next] = top_vec_->A_ptr[i_it]; index_array_[i_next] = index_array_[i_it];
  top_vec_->A_ptr[i_it] = val_temp; index_array_[i_it] = index_temp;
  if (i_next < (top_vec_->M-1)) min_mag_elements_recursive(top_vec_, index_array_, i_next+1);
}

void AYmat::Proj_1(AYmat *z_, int * ind_vec_, double R_)
{
  if ((z_->M*z_->N) == (M*N))
  {
    double tau, tau_test, acc=0.0, val_it, y_it;
    int K=0, k;
    max_mag_elements_ordered(z_, ind_vec_);
    do
    {
      val_it = fabs(z_->A_ptr[K]);
      acc += val_it;
      tau_test = (acc-R_)/((double)(K+1));
      tau = (tau_test < val_it) ? tau_test : tau;
      K++;
    } while ((tau_test < val_it) && (K<M*N));
    if (tau < 0.0) // if already in the unit ball
    {for (int i = 0; i < N*M; i++) z_->A_ptr[i] = A_ptr[i];}
    else
    {
      for (k = 0; k < K-1; k++) // the previous iteration gives us the index where this stops
      {
        val_it =  (A_ptr[ind_vec_[k]]) - tau;
        z_->A_ptr[ind_vec_[k]] = ((A_ptr[ind_vec_[k]] > 0.0) ? val_it : -1.0*val_it);
      }
      for ( k = K-1; k < M*N; k++) z_->A_ptr[ind_vec_[k]] = 0.0;
    }
  }
  else printf("AYmat: Proj_1 failed, dimension mismatch\n");
}

void AYmat::svd(gsl_vector * S_in, gsl_matrix * V_in, gsl_vector * work) {GSL_init(); gsl_linalg_SV_decomp(A_gsl, V_in, S_in, work);}

void AYmat::svd_check()
{
  int i, j;
  AY_SVDspace space(M, N);
  for ( i = 0; i < M; i++) { for ( j = 0; j < N; j++) { gsl_matrix_set(space.U, i, j, AT[j][i]); }}
  space.svd();
  printf("\nU:\n");
  for ( i = 0; i < M; i++)
  {
    for ( j = 0; j < N; j++)
    {
      printf("%e ", gsl_matrix_get(space.U, i, j));
    }
    printf("\n");
  }
  printf("\nV:\n");
  for ( i = 0; i < N; i++)
  {
    for ( j = 0; j < N; j++)
    {
      printf("%e ", gsl_matrix_get(space.V, i, j));
    }
    printf("\n");
  }
  printf("\ns:\n");
  for ( i = 0; i < N; i++) printf("%e\n", gsl_vector_get(space.s, i));
}
