#ifndef AYLINALG_HH
#define AYLINALG_HH

#include <cstring>
#include <cstdint>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

class AYrng
{
    public:
        AYrng();
        AYrng(uint64_t seed_);
        AYrng(uint64_t seed_, uint64_t jump_);
        ~AYrng();

        uint64_t seed=5.555e18;
        uint64_t jump=100;
        uint64_t carry;

        gsl_rng * gsl_gen=NULL;

        void rng_true_init(uint64_t seed_=8.888e18, uint64_t jump_=100);
        virtual void rng_init(uint64_t seed_=0);
        void rng_init_gsl(uint64_t seed_=1);
        virtual double rand_gen() {return 0.0;}
        virtual double rand_gen_gsl() {return 0.0;}
        double rand_uni_gsl(double low_, double high_);
        double rand_gau_gsl(double mu_, double var_);
        double dseed();
        double dcarry();
};

class AYuniform : public AYrng
{
  public:
    double low, high;
    AYuniform(double low_=0.0, double high_=1.0);
    AYuniform(double low_, double high_, uint64_t seed_);
    ~AYuniform();
    double rand_gen();
    inline double rand_gen_gsl() {return rand_uni_gsl(low, high);};
};

class AYnormal : public AYrng
{
  public:
    double mu, var;
    AYnormal(double mu_=0.0, double var_=1.0);
    AYnormal(double mu_, double var_, uint64_t seed_);
    ~AYnormal();
    double rand_gen();
    inline double rand_gen_gsl() {return rand_gau_gsl(mu, var);}
};

class AYmat
{
    public:
      AYmat(int M_, int N_);
      AYmat(char * name);
      AYmat(const char * name): AYmat((char*) name) {}
      AYmat();
      ~AYmat();

      int M, N;
      double * A_ptr;
      double ** AT;
      bool dmatrix_alloc_flag = false;
      gsl_matrix * A_gsl;
      gsl_vector * v_gsl;

      virtual void set(int i, int j, double val);
      virtual double get(int i, int j);

      virtual void GSL_init();
      void GSL_send();

      void GSL_2_AYmat_copy(gsl_matrix * mat_in);
      void GSL_2_AYmat_copy(gsl_vector * vec_in);
      void AYmat_2_GSL_copy(gsl_matrix * mat_in);
      void AYmat_2_GSL_copy(gsl_vector * vec_in);

      void print_mat(bool space_ = true);
      void print_mat(char * name_, bool space_ = true) {printf("%s\n", name_); print_mat(space_);};
      void print_mat(const char * name_, bool space_ = true) {print_mat((char *) name_, space_);};
      void print_dims();
      void fprintf_mat(char name[], bool verbose = false);
      void fprintf_mat(const char name[], bool verbose = false) {fprintf_mat((char*) name, verbose);}
      void init_123();
      void init_0();
      void init_randuni(double low_=0.0, double high_=1.0);
      void init_randgen(AYrng * rng_in);
      void copy_set_col(AYmat * mat_out, int local_index, int out_index);

      virtual AYmat * copy_gen();
      virtual void copy_set(AYmat * X);
      virtual void copy_set(AYmat * X, double scal_);

      virtual AYmat * row_slice_gen(int i_);

      virtual AYmat * col_slice_gen(int j_);

      virtual AYmat * transpose_gen();
      void transpose_set(AYmat * XT);

      void add(AYmat * B, AYmat * C, double alpha=1.0, double beta=0.0);
      void add(double alpha_, AYmat * out_);
      void scal_mult(double c);

      AYmat * mult_gen(AYmat * B, double alpha=1.0, bool trans1_=false, bool trans2_=false);
      void mult_set(AYmat * B, AYmat * C, double alpha=1.0, double beta=0.0, bool trans1_=false, bool trans2_=false);

      double mean();
      double variance();

      double inner(AYmat *B);
      double norm_frob();
      double norm_frob(AYmat * X_);
      double norm_1();
      double max_val_mag();
      double max_val_mag(int * index, int start_index=0);
      double min_val_mag();
      double min_val_mag(int * index, int start_index=0);
      int norm_0(double tol = 1e-10);

      int max_mag_elements(AYmat *top_vec_, int *index_array_);
      void max_mag_elements_ordered(AYmat *top_vec_, int *index_array_);

      int min_mag_elements(AYmat *top_vec_, int *index_array_);
      void min_mag_elements_ordered(AYmat *top_vec_, int *index_array_);

      void Proj_1(AYmat *z_, int * ind_vec_, double R_=1.0);

      virtual void svd(gsl_vector *S, gsl_matrix *V, gsl_vector *work);
      virtual void svd_check();

    protected:
      int GSL_flag = 0;
      friend class AYcolstack;
      friend class DCT_mapping;
      friend class AY_SVDspace;
      friend class AYtens;
      friend class AYsym;
      friend class data_cloud;

    private:
      void max_mag_elements_recursive(AYmat *top_vec_, int * index_array_, int i_next);
      void min_mag_elements_recursive(AYmat *top_vec_, int * index_array_, int i_next);
};

class AYvec : public AYmat
{
  public:
    AYvec(int M_);
    AYvec(char * name_);
    AYvec(const char * name): AYvec((char*)name) {}

    ~AYvec();

    void set(int i, double val);
    double get(int i);

    void GSL_2_AYvec_copy(gsl_matrix * mat_in);
    void GSL_2_AYvec_copy(gsl_vector * vec_in);
    void AYvec_2_GSL_copy(gsl_matrix * mat_in);
    void AYvec_2_GSL_copy(gsl_vector * vec_in);

    void print_vec(bool space_ = true) {print_mat(space_);};
    void print_vec(char * name_ ,bool space_ = true) {printf("%s\n", name_); print_vec(space_);}
    void print_vec(const char * name_ ,bool space_ = true) {print_vec((char*)name_, space_);}
    void fprintf_vec(char * name_, bool space_ = true);
    void fprintf_vec(const char * name_, bool space_ = true) {fprintf_vec((char*) name_, space_);}

    AYvec * copy_gen();
    AYvec * row_slice_gen(int i_);
    AYvec * col_slice_gen(int j_);
    AYmat * transpose_gen();

    double norm_2() {return norm_frob();}
    double norm_2(AYvec * x_) {return norm_frob(x_);}
    double dot(AYmat *B_);
    void svd(gsl_vector *S, gsl_matrix *V, gsl_vector *work);
    void svd_check();
};

class AYcolstack : public AYmat
{
  public:
    AYcolstack(AYmat * A_in_): AYmat()
    {
      M = (A_in_->M)*(A_in_->N);
      N = 1;
      AT = A_in_->AT;
      A_ptr = A_in_->A_ptr;
    }
    ~AYcolstack()
    {}
};

class AYtens
{
  public:
    int W, M, N;
    double *** T_AT;
    AYmat * mat;

    AYtens(int W_, int M_, int N_);
    AYtens(char * name);
    AYtens(const char * name): AYtens((char*)name) {}
    ~AYtens();
    double get(int i_, int j_, int k_);
    void set(int i_, int j_, int k_, double val_);
    void print_dims();
    void print_tens(bool space_ = true);
    void fprintf_tens(char name_[], int split_=1, bool verbose_=false);
    void fprintf_tens(const char name_[], int split_=1, bool verbose_=false) {fprintf_tens((char*)name_, split_, verbose_);}
    void init_0();
    void init_123();
    void init_mats123();

  private:
};

class AYdiag
{
  public:
    int N;

    double * A;

    AYdiag(int N_);
    AYdiag(AYmat * m_);
    AYdiag(AYvec * v_);
    ~AYdiag();

    void print_mat(bool space_ = true);
    void fprintf_diag(char * name_, bool verbose_ = false);
    void fprintf_diag(const char * name_, bool verbose_= false) {fprintf_diag((char*)name_, verbose_);}
    void init_eye();
    void init_123();
    void init_mat(AYmat * m_ );
    void init_vec(AYvec * v_);
    void mult_vec(AYvec * in_, AYvec * out_, bool diff_ = false);

};

class AYsym
{
  public:
    int N, len;

    double ** A;

    AYsym(int N_);
    ~AYsym();

    void print_mat(bool space_ = true);
    void print_mat(char * name_, bool space_ = true) {printf("%s\n", name_); print_mat(space_);}
    void print_mat(const char * name_, bool space_ = true) {print_mat((char*)name_, space_);}
    void fprintf_sym(char name[], bool verbose_=false);
    void fprintf_sym(const char name[], bool verbose_=false) {fprintf_sym((char*)name, verbose_);}
    void init_eye();
    void init_123();
    void init_sqrmat(AYmat * m_ );
    void init_XTWX(AYmat * m_, AYdiag * w_);
    void mult_vec(AYvec * in_, AYvec * out_, bool diff_ = false);
    double vT_A_v(AYvec *v, AYvec * w);
};

class AYdata
{
  public:
    int Frames;
    int depth;
    int ** dims=NULL;
    AYdata(int Frames_, int depth_);
    AYdata();
    ~AYdata();
    virtual void set_dims();
    virtual void fprintf_split(char name_[], bool verbose_=false);

};

class AY_SVDspace // for now, assuming that we are working with a long thin matrix. A relatively simple conditional can be implemented in order to handle short fat case
{
    public:
      AY_SVDspace(AYmat * mat_);
      AY_SVDspace(int M_, int N_);
      ~AY_SVDspace();

      int M_in, N_in;
      gsl_matrix * U;
      gsl_matrix * V;
      gsl_vector * s;
      gsl_vector * work;

      void svd();
      void load_U(AYmat * X_);
      void unpack(AYmat * U_, AYmat * S_, AYmat * V_);
};

class AY_Choleskyspace
{
    public:
      AY_Choleskyspace(AYsym * mat_);
      AY_Choleskyspace(int N_);
      ~AY_Choleskyspace();

      gsl_matrix * mat_gsl;
      gsl_vector * x_gsl=NULL;

      int N_in;
      void load_mat(AYsym * mat_);
      void load_mat(AYsym * mat_, double scal_);
      void Cholesky_decomp();
      void Cholesky_decomp(AYsym * mat_, AYsym * L_);
      void iCholesky_decomp(AYsym * mat_, AYsym * L_, double threshold_ = 1e-1);
      void unpack(AYsym * L_);
      void alloc_workspace();
      void solve_system(AYvec* x_in, AYvec * b_in);
      void solve_system(AYsym * A_, AYvec* x_in, AYvec * b_in);
      void solve_system(AYvec* x_in);
};

void AYlinalg_svd(AYmat * mat_, AY_SVDspace * space_);
void AYlinalg_Cholesky_solve(AYsym * L_, AYvec *z_, AYvec * r_ );

AYvec * AYmat_2_AYvec_gen(AYmat * X_in);
void AYmat_2_AYvec_copy(AYmat * X_in, AYvec * x_in);
AYmat * AYvec_2_AYmat_gen(AYvec * x_in);
void AYvec_2_AYmat_copy(AYvec * x_in, AYmat * X_in);

AYvec * GSL_2_AYvec_gen(gsl_matrix * mat_in);
AYvec * GSL_2_AYvec_gen(gsl_vector * vec_in);
AYmat * GSL_2_AYmat_gen(gsl_matrix * mat_in);
AYmat * GSL_2_AYmat_gen(gsl_vector * vec_in);
AYmat * GSL_2_diagAYmat_gen(gsl_vector * vec_in);

char * string_gen_pruned(char * in_);
char * string_gen_pruned(const char * in_);
void AYfatalerror(char * msg, int exit_code=1);
void AYfatalerror(const char *msg, int exit_code=1);


#endif
