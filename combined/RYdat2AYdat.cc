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

ODR_struct::ODR_struct(char *rydat_loc_, char *rydat_dir_, char *file_name_, int Frames_): AYdata(Frames_, 2), rydat_loc(string_gen_pruned(rydat_loc_)), rydat_dir(string_gen_pruned(rydat_dir_)), file_name(string_gen_pruned(file_name_))
{
    size_t in_buf_len = (size_t)(strlen(rydat_loc)+strlen(rydat_dir)+strlen(file_name) + 50);
    in_buf = new char[in_buf_len]; sprintf(in_buf, "%s%s%s", rydat_loc, rydat_dir, file_name);
    ibuf_end = strlen(in_buf);

    char *file_it = in_buf + ibuf_end; //address at the end of the input string
    sprintf(file_it, ".%d", 0);

    std::ifstream ry_dat;
    ry_dat.open(in_buf);
    int lines = 0;
    for(std::string line; std::getline(ry_dat, line);) ++lines; // determining the number of line in the file
    P=lines-1;
    ry_dat.close();
    set_dims();

    data = AYd3tensor(Frames, P, len_dat);
    specs = AYdmatrix(Frames, len_specs);

    double t, cx, cy, wall_sca;
    double x, y, z, q1, q2, q3, q4;

    for ( int i = 0; i < Frames; i++)
    {
      sprintf(file_it, ".%d", i);
      std::ifstream ry_dat;
      ry_dat.open(in_buf);
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

ODR_struct::ODR_struct(const char * proc_loc_, const char * rydat_dir_, const char * file_name_): AYdata(), rydat_dir(string_gen_pruned(rydat_dir_)), file_name(string_gen_pruned(file_name_))
{
  size_t proc_len = (size_t)(strlen(proc_loc_)+strlen(rydat_dir)+1);
  proc_dir = new char[proc_len]; sprintf(proc_dir, "%s%s", proc_loc_, rydat_dir); mkdir(proc_dir, S_IRWXU);

  size_t out_buf_len = (size_t)(strlen(proc_dir)+strlen(file_name)+50);
  out_buf = new char[out_buf_len]; sprintf(out_buf, "%s%s", proc_dir, file_name);
  obuf_end = strlen(out_buf);

  Frames = 0; depth = 2; dims = AYimatrix(depth, 3);
}

ODR_struct::ODR_struct(const char * filin_dir_): AYdata(), filin_dir(string_gen_pruned(filin_dir_)), reading_flag(true)
{
  size_t filin_dir_len = (size_t)(strlen(filin_dir) + 50);
  in_buf = new char[filin_dir_len]; strcpy(in_buf, filin_dir);
  ibuf_end = strlen(in_buf);

  char * file_it = in_buf + ibuf_end;
  sprintf(file_it, ".fisml");
  int val1, val2, val3;
  std::ifstream sml_stream; std::string line; sml_stream.open(in_buf);
  std::getline(sml_stream, line); std::istringstream dimline0(line);
  dimline0 >> val1 >> depth >> val3;
  std::getline(sml_stream, line); std::istringstream dimline1(line);
  dimline1 >> val1 >> val2 >> Frames;
  std::getline(sml_stream, line); std::istringstream dimline2(line);
  dimline2 >> val1 >> P >> val3;

  sml_stream.close();
}

ODR_struct::~ODR_struct()
{
  if (data!=NULL) free_AYd3tensor(data);
  if (specs!=NULL) free_AYdmatrix(specs);
  if (in_buf!=NULL) delete in_buf;
  if (out_buf!=NULL) delete out_buf;
  if (rydat_loc!=NULL) delete rydat_loc;
  if (rydat_dir!=NULL) delete rydat_dir;
  if (file_name!=NULL) delete file_name;
  if (proc_dir!=NULL) delete proc_dir;
}

void ODR_struct::set_dims()
{
  dims[0][0] = 1; dims[0][1] = 1; dims[0][2] = len_specs;
  dims[1][0] = 1; dims[1][1] = P; dims[1][2] = len_dat;
}

void ODR_struct::fprintf_split(bool verbose_)
{
  char * file_it = out_buf + obuf_end;

  for (int i = 0; i < Frames; i++)
  {
    sprintf(file_it, ".%d.aydat", i);
    FILE * data_file = fopen(out_buf, "wb");
    fwrite(specs[i], sizeof(double), dims[0][1]*dims[0][2], data_file);
    fwrite(data[i][0], sizeof(double), dims[1][1]*dims[1][2], data_file);
    fclose(data_file);
  }

  sprintf(file_it, ".rysml");
  FILE * aysml_file = fopen(out_buf, "w");
  fprintf(aysml_file, "%d %d %d\n", 1, Frames, depth);
  for (int i = 0; i < depth; i++) fprintf(aysml_file, "%d %d %d\n", dims[i][0], dims[i][1], dims[i][2]);
  fclose(aysml_file);

  if (verbose_) printf("ODR_struct: wrote file(s) %s%s.aydat/rysml\n", proc_dir, file_name);
}

void ODR_struct::prepare_datdir(char *proc_loc_)
{
  size_t proc_len = (size_t)(strlen(proc_loc_)+strlen(rydat_dir)+1);
  proc_dir = new char[proc_len]; sprintf(proc_dir, "%s%s", proc_loc_, rydat_dir); mkdir(proc_dir, S_IRWXU);
  size_t obuf_len = (size_t)(strlen(proc_dir)+strlen(file_name)+50);
  out_buf = new char[obuf_len]; sprintf(out_buf, "%s%s", proc_dir, file_name);
  obuf_end = strlen(out_buf);

  printf("made directory  %s\n", proc_dir);
}

void ODR_struct::write_split(int k_, double * specs_vec_, particle * q_)
{
  char * file_it = out_buf + obuf_end;

  if (writing_flag)
  {
    sprintf(file_it, ".%d.aydat", k_);
    FILE * data_out = fopen(out_buf, "wb");
    fwrite(specs_vec_, sizeof(double), (size_t)len_specs, data_out);

    for (int i = 0; i < P; i++)
    {
      double data_vector[] = {q_[i].x, q_[i].y, q_[i].z, q_[i].q0, q_[i].q1, q_[i].q2, q_[i].q3};
      fwrite(data_vector, sizeof(double), (size_t)len_dat, data_out);
    }
    fclose(data_out);
    Frames++;
  }
}

void ODR_struct::write_split(int k_, double * specs_vec_, particle * q_, double ctheta_)
{
  char * file_it = out_buf + obuf_end;

  if (writing_flag)
  {
    sprintf(file_it, ".%d.aydat", k_);
    FILE * data_out = fopen(out_buf, "wb");
    fwrite(specs_vec_, sizeof(double), (size_t)len_specs, data_out);

    double cx_loc = specs_vec_[1]; double cy_loc = specs_vec_[2];

    for (int i = 0; i < P; i++)
    {
      double data_vector[] = {q_[i].x, q_[i].y, q_[i].z, q_[i].q0, q_[i].q1, q_[i].q2, q_[i].q3};
      fwrite(data_vector, sizeof(double), (size_t)len_dat, data_out);
      xs.push_back(cx_im+cl_im*(q_[i].x-cx_loc));
      xs.push_back(cy_im+cl_im*(q_[i].y-cy_loc));
    }
    fclose(data_out);
    ts.push_back(specs_vec_[0]*t_phys);
    d_ang.push_back(ctheta_);

    Frames++;
  }
}

void ODR_struct::end_writing(bool verbose_)
{
  char * file_it = out_buf + obuf_end;

  sprintf(file_it, ".rysml");
  FILE * aysml_file = fopen(out_buf, "w");
  fprintf(aysml_file, "%d %d %d\n", 1, Frames, depth);
  for (int i = 0; i < depth; i++) fprintf(aysml_file, "%d %d %d\n", dims[i][0], dims[i][1], dims[i][2]);
  fclose(aysml_file);

  if (make_filter_inputs_flag)
  {
    double vidspecs_vec[] = {t_phys, cx_im, cy_im, cl_im};
    sprintf(file_it, ".filin");
    FILE * filin_dat = fopen(out_buf, "wb");
    fwrite(&ts[0],sizeof(double),Frames,filin_dat);
    fwrite(&xs[0],sizeof(double),2*P*Frames,filin_dat);
    fwrite(&d_ang[0],sizeof(double),Frames,filin_dat);
    fwrite(vidspecs_vec,sizeof(double),4,filin_dat);
    fclose(filin_dat);

    sprintf(file_it, ".fisml");
    FILE * filin_sml = fopen(out_buf, "wb");
    fprintf(filin_sml, "%d %d %d\n", 0, 1, 4);
    fprintf(filin_sml, "%d %d %d\n", 1, 1, Frames);
    fprintf(filin_sml, "%d %d %d\n", Frames, P, 2);
    fprintf(filin_sml, "%d %d %d\n", 1, 1, Frames);
    fprintf(filin_sml, "%d %d %d\n", 1, 1, 4);
    fclose(filin_sml);
  }

  writing_flag = false;
}

void ODR_struct::set_vidspecs(double t_phys_, double cx_im_, double cy_im_, double cl_im_)
{
  t_phys = t_phys_;
  cx_im = cx_im_;
  cy_im = cy_im_;
  cl_im = cl_im_;
}

void ODR_struct::print_time_rotation()
{for(int i=0;i<Frames;i++) printf("%f %f\n",ts[i],d_ang[i]);}

void ODR_struct::load_filter(double *ts_, double *xs_, double *d_ang_, int offset_)
{
  if (reading_flag)
  {
    int nsnaps = Frames - offset_;
    char * file_name = in_buf + ibuf_end;
    sprintf(file_name, ".filin");
    FILE * inputs = fopen(in_buf, "r"); // or "rb"?

    fseek_safe(inputs, sizeof(double)*offset_, SEEK_SET);
    fread_safe(ts_, sizeof(double), nsnaps, inputs);

    fseek_safe(inputs, sizeof(double)*2*offset_*P, SEEK_CUR);
    fread_safe(xs_, sizeof(double), 2*nsnaps*P, inputs);

    fseek_safe(inputs, sizeof(double)*offset_, SEEK_CUR);
    fread_safe(d_ang_, sizeof(double), nsnaps, inputs);
    fclose(inputs);
  }
  else printf("ODR_struct: not staged for reading binary filter inputs\n");
}

void ODR_struct::read_filin(int offset_)
{
  if (reading_flag)
  {
    int nsnaps = Frames - offset_;
    char * file_name = in_buf + ibuf_end;

    AYvec t_vec(nsnaps), d_vec(nsnaps);
    AYtens x_tens(nsnaps, 2, P);

    double *ts_ = t_vec.A_ptr;
    double *d_ang_ = d_vec.A_ptr;
    double *xs_ = x_tens.T_AT[0][0];

    sprintf(file_name, ".filin");
    FILE * inputs = fopen(in_buf, "r"); // or "rb"?

    fseek_safe(inputs, sizeof(double)*offset_, SEEK_SET);
    fread_safe(ts_, sizeof(double), nsnaps, inputs);

    fseek_safe(inputs, sizeof(double)*2*offset_*P, SEEK_CUR);
    fread_safe(xs_, sizeof(double), 2*nsnaps*P, inputs);

    fseek_safe(inputs, sizeof(double)*offset_, SEEK_CUR);
    fread_safe(d_ang_, sizeof(double), nsnaps, inputs);
    fclose(inputs);

    x_tens.mat[0].print_mat();

  }
  else printf("ODR_struct: not staged for reading binary filter inputs\n");
}
