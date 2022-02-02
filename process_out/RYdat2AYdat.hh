#ifndef RYDAT2AYDAT_HH
#define RYDAT2AYDAT_HH

#include "AYlinalg.hh"

class ODR_struct
{
    public:
        AYtens * data;
        AYmat * specs;

        int len;
        char name[50];
        char rydat_dir[75];

        ODR_struct(char name_[], int len_=1201);
        ~ODR_struct();

};

#endif
