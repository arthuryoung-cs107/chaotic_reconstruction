#include "filter.hh"
#include "common.hh"

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

/** Sets up a directory for saving simulation output.
 * \param[in] fflags_ the types of file output to store.
 * \param[in] odir_ the name of the directory.
 * \param[in] state_freq_ the frequency at which to store state information (if
 *                        used). */
void filter::setup_output_info(unsigned int fflags_,const char *odir_,int state_freq_) {
    fflags=fflags_;state_freq=state_freq_;
    size_t l=strlen(odir_)+1;

    // Allocate space for the directory filename and copy it. Allocate additional
    // space for assembling the output filenames.
    odir=new char[2*l+64];
    memcpy(odir,odir_,sizeof(char)*l);
    obuf=odir+l;

    // Make the output directory if it doesn't already exist
    mkdir(odir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Open the digest file if it is in use
    if(fflags&128) {
        sprintf(obuf,"%s/digest",odir);
        fdigest=safe_fopen(obuf,"wb");
    }
}

/** Writes a selection of files about the swirling simulations and the
 * associated experimental data. */
void filter::write_files() {

    // Check for the case when no output flags have been set
    if(fflags==0) return;

    // Output the most probable particle, plus the associated experimental data
    if(fflags&3) {
        swirl *swa=sw[most_likely()];
        if(fflags&1) output_particle(swa,"pa"); // mode particle sample
        if(fflags&2) output_data(swa,"ea");
    }

    // Output the fourth-power weighted average of the particles, plus the
    // associated experimental data
    if(fflags&12) {
        swirl swb(sw,npar,4);
        if(fflags&4) output_particle(&swb,"pb"); // something in between, check the weighting in output_particle
        if(fflags&8) output_data(&swb,"eb");
    }

    // Output the weighted average of the particles, plus the
    // associated experimental data
    if(fflags&48) {
        swirl swc(sw,npar,1);
        if(fflags&16) output_particle(&swc,"pc"); // mean particle sample. most visually compelling
        if(fflags&32) output_data(&swc,"ec");
    }

    // Output the state information
    if((fflags&64)&&frame%state_freq==0) output_states();

    // Save the digest information
    write_digest();
}

/** Saves a snapshot of the bead positions and rotations in a simulation.
 * \param[in] swp a pointer to the simulation to save.
 * \param[in] prefix the filename prefix. */
void filter::output_particle(swirl *swp,const char* prefix) {
    sprintf(obuf,"%s/%s.%d",odir,prefix,frame);
    FILE *fp=safe_fopen(obuf,"w");
    swp->output(fp);
    fclose(fp);
}

/** Saves a snapshot of experimental data, converted into simulation units.
 * \param[in] swp a pointer to the simulation to obtain the length conversion
 *                factors from. */
void filter::output_data(swirl *swp,const char* prefix) {

    // Open the output file
    if(odir==NULL) return;
    sprintf(obuf,"%s/%s.%d",odir,prefix,frame);
    FILE *fp=safe_fopen(obuf,"w");

    // Output header information
    fprintf(fp,"%g %g %g %g\n",swp->time,swp->cx,swp->cy,swp->wall_sca);

    // Output the experimental bead positions
    float *f=xs+2*n*frame;
    double sx=swp->cx_im,sy=swp->cy_im,icl=1./swp->cl_im;
    for(int i=0;i<n;i++,f+=2)
        fprintf(fp,"%g %g %g\n",(*f-sx)*icl+swp->cx,(f[1]-sy)*icl+swp->cy,swp->rad);
    fclose(fp);
}

/** Saves all of the parameters for all of the particles to a file in binary
 * format. */
void filter::output_states() {

    // Open the output file
    if(odir==NULL) return;
    sprintf(obuf,"%s/st.%d",odir,frame);
    FILE *fp=safe_fopen(obuf,"wb");

    // Output the states of the all the simulations, and close the output file
    for(int k=0;k<npar;k++) sw[k]->output_state(fp);
    fclose(fp);
}

/** Writes the means and standard deviations of the parameters to the file. */
void filter::write_digest() {

    // Allocate and clear the global digest array
    double gdig[2*sp_num_vparams];
    for(int i=0;i<2*sp_num_vparams;i++) gdig[i]=0.;

#pragma omp parallel
    {

        // Clear local accumulators for this thread
        double tdig[2*sp_num_vparams];
        for(int i=0;i<2*sp_num_vparams;i++) tdig[i]=0.;

        // Add the contributions from each simulation
#pragma omp for
        for(int k=0;k<npar;k++) sw[k]->accumulators(tdig);

        for(int i=0;i<2*sp_num_vparams;i++) {
#pragma omp atomic
            gdig[i]+=tdig[i];
        }
    }

    // Assemble the output array, converting to single precision
    float odig[2*sp_num_vparams+3],*op=odig+1;
    *odig=ts[frame];
    for(int i=0;i<sp_num_vparams;i++) {
        op[i]=gdig[i];
        op[sp_num_vparams+i]=sqrt(gdig[sp_num_vparams+i]-gdig[i]*gdig[i]);
    }
    odig[2*sp_num_vparams+1]=min_l2;
    odig[2*sp_num_vparams+2]=min_linf;

    // Write the results to the digest file
    fwrite(odig,sizeof(float),2*sp_num_vparams+3,fdigest);
    fflush(fdigest);
}
