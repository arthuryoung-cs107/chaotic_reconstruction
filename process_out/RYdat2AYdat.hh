#ifndef RYDAT2AYDAT_HH
#define RYDAT2AYDAT_HH

#include "AYlinalg.hh"

class ODR_struct : public AYdata
{
    public:
        AYtens * data;
        AYmat * specs;
        int counter=0;

        bool writing = false;

        char name[50];
        char rydat_dir[75];

        int P;

        ODR_struct(char name_[], int len_=1201);
        ~ODR_struct();
        void set_dims();
};

#endif
