#include <sys/stat.h>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>

#include "RYdat2AYdat.hh"

extern "C"
{
  #include "AYaux.h"
}

ODR_struct::ODR_struct(char name_[], int Frames_): AYdata(Frames_, 2)
{
    memset(name, 0, 49); snprintf(name, 50, "%s", name_);
    memset(rydat_dir, 0, 74); snprintf(rydat_dir, 75, "../twin/%s", name);
    char file_it[100]; memset(file_it, 0, 99); snprintf(file_it, 100, "%spts.0", rydat_dir);

    std::ifstream ry_dat;
    ry_dat.open(file_it);
    int lines = 0;
    for(std::string line; std::getline(ry_dat, line);) ++lines; // determining the number of line in the file
    P=lines-1;
    ry_dat.close();
    set_dims();

    data = AYd3tensor(Frames, P, len_dat); data_alloc_flag(true);
    specs = AYdmatrix(Frames, len_specs); specs_alloc_flag(true);

    double t, cx, cy, wall_sca;
    double x, y, z, q1, q2, q3, q4;

    for ( int i = 0; i < Frames; i++)
    {
      memset(file_it, 0, 99); snprintf(file_it, 100, "%spts.%d", rydat_dir, i);
      std::ifstream ry_dat;
      ry_dat.open(file_it);
      int count = 0;
      for(std::string line; std::getline(ry_dat, line);)
      {
        std::istringstream in(line);
        if (count == 0)
        {
          in >> t >> cx >> cy >> wall_sca;
          specs[i][0] = t;
          specs[i][1] = cx;
          specs[i][2] = cy;
          specs[i][3] = wall_sca;
        }
        else
        {
          in >> x >> y >> z >> q1 >> q2 >> q3 >> q4;
          data[i][count-1][0] = x;
          data[i][count-1][1] = y;
          data[i][count-1][2] = z;
          data[i][count-1][3] = q1;
          data[i][count-1][4] = q2;
          data[i][count-1][5] = q3;
          data[i][count-1][6] = q4;
        }
        count++;
      }
      ry_dat.close();
    }
}

ODR_struct::ODR_struct(const char *full_dir_): AYdata()
{
  size_t l=strlen(full_dir_)+1;

  // Allocate space for the directory filename and copy it. Allocate additional
  // space for assembling the output filenames.
  full_dir=new char[2*l+32];
  memcpy(odir,odir_,sizeof(char)*l);
  obuf=odir+l;

  // Make the output directory if it doesn't already exist
  mkdir(odir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

}

ODR_struct::~ODR_struct()
{
  free_AYd3tensor(data); free_AYdmatrix(specs);
}

void ODR_struct::set_dims()
{
  dims[0][0] = 1; dims[0][1] = 1; dims[0][2] = len_specs;
  dims[1][0] = 1; dims[1][1] = P; dims[1][2] = len_dat;
}

void ODR_struct::AYdata_rysml_gen(char name_[], int split_)
{
  char smlfile[300]; memset(smlfile, 0, 299); snprintf(smlfile, 300, "%s.rysml", name_);
  FILE * aysml_file = fopen(smlfile, "w");
  fprintf(aysml_file, "%d %d %d \n", split_, Frames, depth);
  for (int i = 0; i < depth; i++)
  {
    fprintf(aysml_file, "%d %d %d\n", dims[i][0], dims[i][1], dims[i][2]);
  }
  fclose(aysml_file);
}


void ODR_struct::fprintf_split(char name_[], bool verbose_)
{
  char aysml_name[300]; memset(aysml_name, 0, 299); snprintf(aysml_name, 300, "%s%s", directory, name_);
  char specfile[300];

  for (int i = 0; i < Frames; i++)
  {
    memset(specfile, 0, 299); snprintf(specfile, 300, "%s%s.%d.aydat", directory, name_, i);
    FILE * data_file = fopen(specfile, "wb");
    fwrite(specs[i], sizeof(double), dims[0][1]*dims[0][2], data_file);
    fwrite(data[i][0], sizeof(double), dims[1][1]*dims[1][2], data_file);
    fclose(data_file);
  }
  AYdata_rysml_gen(aysml_name, 1);
  if (verbose_) printf("AYdata: wrote file(s) %s.aydat/aysml\n", name);
}

void ODR_struct::prepare_datdir(char name_[])
{memset(directory, 0, 99); snprintf(directory, 100, "%s%s", name_, name); mkdir(directory, S_IRWXU);
printf("made directory  %s\n", directory);}
