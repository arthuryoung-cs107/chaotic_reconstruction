#ifndef RYDAT2AYDAT_HH
#define RYDAT2AYDAT_HH

#include "AYlinalg.hh"
#include "particle.hh"
#include "swirl_param.hh"
#include <vector>

int set_special_params(int id_, double *vec_);
int set_special_params(const char *id_, double *vec_);

class ODR_struct : public AYdata
{
    public:
        double ***data=NULL, ** specs=NULL;

        bool  writing_flag = false, reading_flag = false, make_filter_inputs_flag, write_split_flag=true;

        char *rydat_loc=NULL, *rydat_dir=NULL, *proc_dir=NULL, *file_name=NULL, *filin_dir=NULL,
        *in_buf=NULL, *out_buf=NULL;

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
        ODR_struct(): AYdata() {}
        void init_process(char *rydat_loc_, char *rydat_dir_, char *file_name_, int Frames_);
        void init_swirl(char * proc_loc_, char * rydat_dir_, char * file_name_, bool write_split_flag_=true);
          void init_swirl(const char * proc_loc_, const char * rydat_dir_, const char * file_name_, bool write_split_flag_=true)
            {init_swirl((char*)proc_loc_,(char*)rydat_dir_,(char*)file_name_, write_split_flag_);}
        void init_filter(char * filin_dir_);
          void init_filter(const char * filin_dir_)
            {init_filter((char*) filin_dir_);}
        void init_race(char * proc_loc_, char * rydat_dir_, char * file_name_);
          void init_race(const char *proc_loc_, const char *rydat_dir_, const char *file_name_)
            {init_race((char*)proc_loc_,(char*)rydat_dir_,(char*)file_name_);}
        ~ODR_struct();
        void set_dims();
        void prepare_datdir(char name_[]);
        void fprintf_split(bool verbose_=false);
        void write(int k_, double * specs_, particle * q_);
        void write(int k_, double * specs_, particle * q_, double ctheta_);
        void end_writing(bool verbose_=false);
        inline void set_vidspecs(double t_phys_, double cx_im_, double cy_im_, double cl_im_) {t_phys=t_phys_, cx_im=cx_im_, cy_im=cy_im_,cl_im=cl_im_;}
        void print_time_rotation();
        void load_filter(double *ts_, double *xs_, double *d_ang_, int offset_=0);
        void read_filin(int offset_=0);
        void stage_filout() {}
        ODR_struct * spawn_swrlbest(char * name_);
        void write_sparam(swirl_param * sparam_, char * name_);
          void write_sparam(swirl_param * sparam_, const char * name_)
            {write_sparam(sparam_,(char*)name_);}
      private:
        /** The time points of the snapshots. */
        std::vector<double> ts;
        /** The bead positions at the snapshots. */
        std::vector<double> xs;
        /** The array of pointers to swirling simulations. */
        std::vector<double> d_ang;

        FILE * file_ptr=NULL;
};

#endif
