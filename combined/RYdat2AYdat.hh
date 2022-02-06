#ifndef RYDAT2AYDAT_HH
#define RYDAT2AYDAT_HH

#include "AYlinalg.hh"
#include "particle.hh"

class ODR_struct : public AYdata
{
    public:
        double ***data, ** specs;

        bool  writing_flag = false,
          data_alloc_flag = false,
          specs_alloc_flag = false,
          ibuf_alloc_flag = false,
          obuf_alloc_flag = false,
          rydat_loc_alloc_flag = false,
          rydat_dir_alloc_flag = false,
          file_name_alloc_flag = false,
          proc_dir_alloc_flag = false;

        char *rydat_loc, *rydat_dir, *proc_dir, *file_name, *in_buf, *out_buf;

        int P;
        int len_specs = 4;
        int len_dat = 7;
        int counter=0;

        ODR_struct(char *rydat_loc_, char *rydat_dir_, char *file_name_, int Frames_);
        ODR_struct(const char * proc_loc_, const char * rydat_dir_, const char * file_name_);
        ~ODR_struct();
        void set_dims();
        void fprintf_split(bool verbose_=false);
        void prepare_datdir(char name_[]);
        void write_split(int k_, double * specs_, particle * q_);

};

#endif
