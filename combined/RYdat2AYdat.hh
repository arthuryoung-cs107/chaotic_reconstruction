#ifndef RYDAT2AYDAT_HH
#define RYDAT2AYDAT_HH

#include "AYlinalg.hh"

class ODR_struct : public AYdata
{
    public:
        double *** data;
        double ** specs;
        int counter=0;

        bool writing = false;
        bool data_alloc_flag = false;
        bool specs_alloc_flag = false; 

        char name[50];
        char rydat_dir[75];
        char directory[100];
        int P;
        int len_specs = 4;
        int len_dat = 7;

        ODR_struct(char name_[], int len_=1201);
        ODR_struct(const char *odir_);
        ~ODR_struct();
        void set_dims();
        void AYdata_rysml_gen(char name_[], int split_=1);
        void fprintf_split(char name_[], bool verbose_=false);
        void prepare_datdir(char name_[]);
};

#endif
