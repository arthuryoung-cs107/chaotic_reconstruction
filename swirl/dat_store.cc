#include "dat_store.hh"
#include "common.hh"

#include <cstdio>

void dat_store::snapshot(double time,particle *q,int n_,double cx,double cy,double ctheta) {
    if(nsnap++==0) n=n_;
    else if(n!=n_) fatal_error("Particle number mismatch",1);
    ts.push_back(time*t_phys);
    for(int i=0;i<n;i++) {
        xs.push_back(cx_im+cl_im*(q[i].x-cx));
        xs.push_back(cy_im+cl_im*(q[i].y-cy));
    }
    d_ang.push_back(ctheta);
}

void dat_store::write(const char* filename) {
    FILE *fp=safe_fopen(filename,"wb");
    fwrite(&n,sizeof(int),2,fp);
    fwrite(&ts[0],sizeof(float),nsnap,fp);
    fwrite(&xs[0],sizeof(float),2*n*nsnap,fp);
    fwrite(&d_ang[0],sizeof(float),nsnap,fp);
}

void dat_store::print_theta_info() {
    for(int i=0;i<nsnap;i++) printf("%g %g\n",ts[i],d_ang[i]);
}
