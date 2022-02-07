#ifndef RYDAT2AYDAT_HH
#define RYDAT2AYDAT_HH

#include "AYlinalg.hh"
#include "particle.hh"
#include <vector>

class ODR_struct : public AYdata
{
    public:
        double ***data=NULL, ** specs=NULL;

        bool  writing_flag = false, make_filter_inputs_flag;

        char *rydat_loc=NULL, *rydat_dir=NULL, *proc_dir=NULL, *file_name=NULL, *in_buf=NULL, *out_buf=NULL, *filin_dir=NULL;

        size_t ibuf_end, obuf_end;

        int P;
        int len_specs = 4;
        int len_dat = 7;

        double t_phys;
        /** The dish x center in the images. */
        double cx_im;
        /** The dish y center in the images. */
        double cy_im;
        /** The bead diameter in the images. */
        double cl_im;

        ODR_struct(char *rydat_loc_, char *rydat_dir_, char *file_name_, int Frames_);
        ODR_struct(const char * proc_loc_, const char * rydat_dir_, const char * file_name_);
        ODR_struct(const char * filin_dir_);
        ~ODR_struct();
        void set_dims();
        void prepare_datdir(char name_[]);
        void fprintf_split(bool verbose_=false);
        void write_split(int k_, double * specs_, particle * q_);
        void write_split(int k_, double * specs_, particle * q_, double ctheta_);
        void end_writing(bool verbose_=false);
        void set_vidspecs(double t_phys_, double cx_im_=402.0, double cy_im_ = 380.0, double cl_im_= 37.6);
        void print_time_rotation();
      private:
        /** The time points of the snapshots. */
        std::vector<double> ts;
        /** The bead positions at the snapshots. */
        std::vector<double> xs;
        /** The array of pointers to swirling simulations. */
        std::vector<double> d_ang;

};

#endif
