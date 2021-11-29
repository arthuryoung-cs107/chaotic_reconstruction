#include "filter.hh"
#include "common.hh"

#include <cstdlib>
#include <cfloat>

#include <gsl/gsl_randist.h>

#ifdef _OPENMP
#include "omp.h"
#endif

/** Initializes the filtering problem, setting class constants and reading in
 * the experimental data from a file.
 * \param[in] fparam the filtering parameters.
 * \param[in] sp_min_ the minimum swirling simulation parameters.
 * \param[in] sp_max_ the maximum swirling simulation parameters.
 * \param[in] sp_rnd_ the random displacements to the apply to the swirling
 *                    parameters.
 * \param[in] wl_ a list of walls in the simulation.
 * \param[in] t_phys_ the physical time elapsed in one simulation time unit.
 * \param[in] filename the name of the file from which to read in the
 *                     experimental data.
 * \param[in] offset the number of frames in the experimental data to offset
 *                   by. */
filter::filter(fil_param &fparam,swirl_param &sp_min_,swirl_param &sp_max_,swirl_param &sp_rnd_,wall_list &wl_,double t_phys_,const char* filename,int offset)
    : fil_param(fparam), sp_min(sp_min_), sp_max(sp_max_), sp_rnd(sp_rnd_),
    t_phys(t_phys_), min_l2(0.), min_linf(0.), nfail(0), sw(NULL), fflags(0),
#ifdef _OPENMP
    nt(omp_get_max_threads()),
#else
    nt(1),
#endif
    ttab(new int[nt+1]), rloc(new int[nt+1]), rng(new gsl_rng*[nt]), wl(wl_),
    pg(new proximity_grid*[nt]), odir(NULL), fdigest(NULL) {

    // Read in the header information and check it makes sense
    FILE *fp=safe_fopen(filename,"rb");
    safe_fread(&(this->n),sizeof(int),2,fp,"header information");
    if(n<=0) fatal_error("Number of particles must be positive",1);
    nsnap-=offset;
    if(nsnap<=0) fatal_error("Number of snapshots must be positive",1);

    // Set up memory for the snapshot times and read the data from the file
    ts=new float[nsnap];
    safe_fseek(fp,sizeof(float)*offset,SEEK_CUR);
    safe_fread(ts,sizeof(float),nsnap,fp,"snapshot times");

    // Read in the bead (x,y) positions
    xs=new float[2*n*nsnap];
    safe_fseek(fp,sizeof(float)*2*n*offset,SEEK_CUR);
    safe_fread(xs,sizeof(float),2*n*nsnap,fp,"bead (x,y) positions");

    // Read in the angle data
    d_ang=new float[nsnap];
    safe_fseek(fp,sizeof(float)*offset,SEEK_CUR);
    safe_fread(d_ang,sizeof(float),nsnap,fp,"angle data");
    fclose(fp);

    // Set up the random number generators
#pragma omp parallel
    {
        int t=thread_num();
        rng[t]=gsl_rng_alloc(gsl_rng_taus2);
        pg[t]=new proximity_grid();
        gsl_rng_set(rng[t],t+1);
    }

    // Set the zeroth entries of the thread table and resampling location table
    // to zero, since they are always held at this value
    *ttab=*rloc=0;
}

/** The class destructor frees the dynamically allocated memory. */
filter::~filter() {
    if(odir!=NULL) {
        delete [] odir;
        if(fdigest!=NULL) fclose(fdigest);
    }
    if(sw!=NULL) {
        for(int i=0;i<npar;i++) delete sw[i];
        delete [] sw;
        delete [] rsbuf;
    }
    for(int i=0;i<nt;i++) {
        gsl_rng_free(rng[i]);
        delete pg[i];
    }
    delete [] pg;
    delete [] rng;
    delete [] rloc;
    delete [] ttab;
    delete [] d_ang;
    delete [] xs;
    delete [] ts;
}

/** Runs the particle filter for the swirling systems.
 * \param[in] frames the total number of experimental data snapshots to filter
 *                   to. */
void filter::run(int frames) {
    if(frame==0) write_files();

    // Step through the experimental data snapshots
    for(int lframe=frame+frames;frame<lframe;) {

        // Advance the swirling simulations forward
        step_frame();

        // Perform simulation outptut
        write_files();
    }
}

/** Resamples the particle distribution by deleting particles with low weight,
 * and duplicating particles with high weight. */
void filter::resample() {

    // Loop over the existing particles and keep track of their cumulative
    // weight in cw. Check for when the cumulative weight hits a target weight
    // stored in tw.
    int l=0,j=0,p=0,tc=1,*rsbp=rsbuf;
    double cw=0,hw=1./npar,tw=hw*gsl_rng_uniform(*rng);
    while(j<npar) {

        // Add up weights of the existing particles until hitting the target
        // weight
        while(l<npar) {
            cw+=sw[l]->wei;
            if(cw<tw) {
                sw[l]->wei=-1;p++;
                if(++l==ttab[tc]) rloc[tc++]=p;
            } else break;
        }

        // Copy the pointer to particle that matched target weight over to the
        // new array
        tw=hw*(++j+gsl_rng_uniform(*rng));

        // Check to see if this particle also hits subsequent target weights, and
        // if so duplicate the particle in the new array
        while(cw>=tw&&j<npar) {
            *(rsbp++)=l;
            tw=hw*(++j+gsl_rng_uniform(*rng));
        }
        if(++l==ttab[tc]) rloc[tc++]=p;
    }
    while(tc<=nt) rloc[tc++]=p;

    // In parallel, get the threads to replace any deleted simulations with
    // copies from the resampling buffer
#pragma omp parallel
    {
        int t=thread_num(),v=rloc[t];
        for(int l=ttab[t];l<ttab[t+1];l++) {
            if(sw[l]->wei<0)
                sw[l]->copy(*sw[rsbuf[v++]]);
            sw[l]->wei=1./npar;
        }
        if(v!=rloc[t+1]) fatal_error("resampling error",1);
    }
}

/** Initializes the filtering problem, by creating many swirling simulations
 * (a.k.a. "particles") to match the experimental initial conditions.
 * \param[in] npar_ the number of swirling simulations to create. */
void filter::init(int npar_) {

    // Allocate memory for the resampling buffer and the swirling simulations
    rsbuf=new int[npar=npar_];
    sw=new swirl*[npar];
    for(int i=1;i<=nt;i++) ttab[i]=long(i)*long(npar)/long(nt);

    // Initialize the swirling simulations
    frame=0;
    double time=ts[frame]/t_phys;
    float *f=xs+2*n*frame;
#pragma omp parallel
    {
        gsl_rng *r=rng[thread_num()];
#pragma omp for
        for(int i=0;i<npar;i++) {
            swirl_param sp(sp_min,sp_max,r);
            sw[i]=new swirl(sp,NULL,wl,n);
            sw[i]->wei=1./npar;
            sw[i]->init_positions(time,d_ang[frame],f,r);
        }
    }
}

/** Steps all of the particles forward until the next frame marker. */
void filter::step_frame() {
    double dur=(ts[frame+1]-ts[frame])/t_phys,sum_l2=DBL_MAX,min_linf_=DBL_MAX,
           ctheta=d_ang[frame],comega=d_ang[frame+1]-d_ang[frame],sp=0,spp=0;
    if(comega>M_PI) comega-=2*M_PI;else if(comega<-M_PI) comega+=2*M_PI;
    comega/=dur;
    float *f=xs+2*n*++frame;

#pragma omp parallel reduction(min:sum_l2) reduction(min:min_linf_)
    {
        int t=thread_num();
        double lnorm[3];
        gsl_rng *r=rng[t];
        proximity_grid *pgp=pg[t];
#pragma omp for
        for(int i=0;i<npar;i++) {

            // Link the proximity grid of this thread to the swirling
            // simulation and advance it forward by the time to the next frame
            sw[i]->pg=pgp;
            sw[i]->advance(dur,ctheta,comega,dt_sim);

            // Apply random perturbations to the simulation
            sw[i]->jiggle(sp_min,sp_max,sp_rnd,r,v_pert,omega_pert);

            // Calculate the change in weight of the particle, based on how
            // well it matches the bead data
            sw[i]->update_weight(f,dur,t_wheels,lnorm);
            if(lnorm[1]<sum_l2) sum_l2=lnorm[1];
            if(lnorm[2]<min_linf_) min_linf_=lnorm[2];
        }
    }

    // Check for failure if the minimum L2 norm is beyond a complete bead
    // diameter
    min_l2=sqrt(sum_l2/n);min_linf=min_linf_;
    if(min_l2>sp_max.cl_im) {
        puts("# Inaccurate bead positions in all particles");
        if(nfail++==filter_nfail_thresh) {
            fputs("# Too many consecutive inaccurate frames\n",stderr);
            exit(1);
        }
    } else nfail=0;

    // Apply the changes to the weights, shifting the exponential to ensure
    // that at least one applied factor is non-zero
#pragma omp parallel for reduction(+:sp) reduction(+:spp)
    for(int i=0;i<npar;i++) {
        sw[i]->wei*=exp(gau_coeff*(sum_l2-sw[i]->logfac));

        // Compute first and second moments of the weights
        sp+=sw[i]->wei;
        spp+=sw[i]->wei*sw[i]->wei;
    }

    // Normalize the particle weights to sum to one
    double isp=1./sp,neff=sp*sp/spp;
    for(int i=0;i<npar;i++) sw[i]->wei*=isp;

    // If the effective number of particles is less than a fraction of the
    // total particles, then resample the particles
    printf("# Frame %d, n_eff=%g   {%g} {%g}",frame,neff,min_l2,min_linf);
    if(neff<rs_thresh*npar) {puts(" [resample]");resample();}
    else putchar('\n');
}

/** Finds the current most likely swirling simulation.
 * \return The index of the particle with the largest weight. */
int filter::most_likely() {
    int j=0;double wmax=sw[0]->wei;
    for(int i=1;i<npar;i++)
        if(sw[i]->wei>wmax) {wmax=sw[i]->wei;j=i;}
    return j;
}
