#include "swirl.hh"
#include "common.hh"

#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

/** Sets up a directory for saving simulation output.
 * \param[in] odir_ the name of the directory. */
void swirl::setup_output_dir(const char *odir_) {
    size_t l=strlen(odir_)+1;

    // Allocate space for the directory filename and copy it. Allocate additional
    // space for assembling the output filenames.
    odir=new char[2*l+32];
    memcpy(odir,odir_,sizeof(char)*l);
    obuf=odir+l;

    // Make the output directory if it doesn't already exist
    mkdir(odir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
}

/** Saves a snapshot of the particle positions and rotations. */
void swirl::output(int k) {

    // Open the output file
    if(odir==NULL) return;
    sprintf(obuf,"%s/pts.%d",odir,k);
    FILE *fp=safe_fopen(obuf,"w");
    output(fp);
    fclose(fp);

    if(dstore!=NULL) dstore->snapshot(time,q,n,cx,cy,ctheta);
}

void swirl::output(FILE *fp) {

    // Output a header line containing the time, dish position, and wall
    // scaling factor
    fprintf(fp,"%g %g %g %g\n",time,cx,cy,wall_sca);

    // Output the particle positions and rotations
    for(int i=0;i<n;i++) {
        q[i].normalize_q();
        fprintf(fp,"%g %g %g %g %g %g %g\n",q[i].x,q[i].y,q[i].z,q[i].q0,q[i].q1,q[i].q2,q[i].q3);
    }
}

/** Imports the particle positions from a text file, applying the conversion
 * from pixels to simulation units. */
void swirl::import(const char* filename) {
    FILE *fp=fopen(filename,"r");
    double icl=1./cl_im,x,y;
    int id;
    for(int i=0;i<n;i++) {
        if(fscanf(fp,"%d %lf %lf\n",&id,&x,&y)!=3)
            fatal_error("Can't read particle information\n",1);
        q[i].set_pos((x-cx_im)*icl+cx,(y-cy_im)*icl+cy,rad);
    }
    fclose(fp);
}

/** Saves the simulation parameters as single-precision floating point numbers
 * in binary format.
 * \param[in] fp the file hand to write to. */
void swirl::output_state(FILE *fp) {
    float buf[sp_num_vparams+1];
    double *d=&(this->Kn);

    // Assemble the data in the buffer, converting to single precision, and
    // then write it out
    *buf=static_cast<float>(wei);
    for(int k=0;k<sp_num_vparams;k++) buf[k+1]=d[k];
    fwrite(buf,sizeof(float),sp_num_vparams+1,fp);
}

/** Adds the parameters and their squares to an array, for computing the mean
 * and standard deviations of the parameters across all simulations and
 * outputting them to file.
 * \param[in] tdig a pointer to the array to write to. */
void swirl::accumulators(double *tdig) {
    double *d=&(this->Kn);
    for(int k=0;k<sp_num_vparams;k++) {
        tdig[k]+=wei*d[k];
        tdig[sp_num_vparams+k]+=wei*d[k]*d[k];
    }
}
