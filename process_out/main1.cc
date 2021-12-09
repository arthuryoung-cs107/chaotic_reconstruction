#include <sys/stat.h>
#include <sstream>
#include <string>
#include <fstream>

#include "blas.h"
#include "omp.h"
#include "AYlinalg.hh"

extern "C"
{
  #include "nrutil.h"
  #include "knuth_lcg.h"
  #include "auxiliary_functions.h"
}

void preliminary_test1()
{
  AYmat mat1(3, 5);
  mat1.init_123();
  mat1.print_mat();
  char name1[50]; name_gen(name1, 50, "mat1");
  char name2[50]; name_gen(name2, 50, "mat2");
  mat1.fprintf_mat(name1);

  AYmat * mat2 = aysml_read(name2);
  mat2->print_mat();
}
void preliminary_test2()
{
  int W=5;
  int M=4;
  int N=3;
  int i, j, k ;

  double *** T = AYd3tensor(W, M, N);
  for ( i = 0; i < W*M*N; i++) *(**T + i) = (double) i;

  for ( k = 0; k < W; k++)
  {
    printf("k = %d: \n", k);
    for ( i = 0; i < M; i++)
    {
      printf("i = %d: ", i);
      for ( j = 0; j < N; j++)
      {
        int check = (int) T[k][i][j];
        printf("%d ", check);
      }
      printf("\n");
    }
    printf("\n");
  }
  free_AYd3tensor(T);
}
void preliminary_test3()
{
  int W=5;
  int M=4;
  int N=3;
  AYtens T1(W, M, N); T1.init_123();
  AYtens T2(W, M, N); T2.init_mats123();
  AYmat M1(M, N); M1.init_123();

  printf("T1\n");
  T1.print_tens();

  printf("T2\n");
  T2.print_tens();

  printf("M1\n");
  M1.print_mat();

  char name1[50]; name_gen(name1, 50, "./dat_dir/tens1");
  char name2[50]; name_gen(name2, 50, "./dat_dir/mat1");
  M1.fprintf_mat(name2);
  T1.fprintf_tens(name1);

  char name1_back[50]; name_gen(name1_back, 50, "./dat_dir/tens1_back");
  char name2_back[50]; name_gen(name2_back, 50, "./dat_dir/mat1_back");

  printf("T1 back\n");
  AYtens * T1_back = aysml_read_tens(name1_back);
  T1_back->print_tens();
  printf("M1 back\n");
  AYmat * M1_back = aysml_read(name2_back);
  M1_back->print_mat();
}

// void convert_RYdata()
// {
//   char * check;
//   int len = 1200 + 1;
//   char rydat_dir[50]; name_gen(rydat_dir, 50, "../twin/circ6.odr/");
//   char specs_name[100]; name_gen(specs_name, 100, "dat_dir/cir6_specs");
//   char tens_name[100]; name_gen(tens_name, 100, "dat_dir/cir6_data");
//
//   char file_it[100];
//   double t, cx, cy, wall_sca;
//   double x, y, z, q1, q2, q3, q4;
//
//   AYmat specs_mat(len, 4);
//
//   std::ifstream ry_dat;
//   ry_dat.open(file_it);
//   int lines = 0;
//   for(std::string line; std::getline(ry_dat, line);) ++lines;
//   int dat_lines=lines-1;
//   ry_dat.close();
//
//   AYtens rytensor(len, dat_lines, 7);
//
//   for ( int i = 0; i < 2; i++)
//   {
//     memset(file_it, 0, 99); snprintf(file_it, 100, "%spts.%d", rydat_dir, i);
//     std::ifstream ry_dat;
//     ry_dat.open(file_it);
//     int count = 0;
//     for(std::string line; std::getline(ry_dat, line);)
//     {
//       std::istringstream in(line);
//       if (count == 0)
//       {
//         in >> t >> cx >> cy >> wall_sca;
//         specs_mat.set(i, 0, t);
//         specs_mat.set(i, 1, cx);
//         specs_mat.set(i, 2, cy );
//         specs_mat.set(i, 3, wall_sca);
//       }
//       else
//       {
//         in >> x >> y >> z >> q1 >> q2 >> q3 >> q4;
//         rytensor.M_ptr[i]->set(count-1, 0, x);
//         rytensor.M_ptr[i]->set(count-1, 1, y);
//         rytensor.M_ptr[i]->set(count-1, 2, z);
//         rytensor.M_ptr[i]->set(count-1, 3, q1);
//         rytensor.M_ptr[i]->set(count-1, 4, q2);
//         rytensor.M_ptr[i]->set(count-1, 5, q3);
//         rytensor.M_ptr[i]->set(count-1, 6, q4);
//       }
//       count++;
//     }
//     ry_dat.close();
//   }
//   specs_mat.fprintf_mat(specs_name);
//   rytensor.fprintf_tens(tens_name);
// }

void process_circ6()
{

}

int main()
{
  // preliminary_test1();
  // preliminary_test2();
  preliminary_test3();
  // process_circ6();
  // convert_RYdata();

  return 0;
}
