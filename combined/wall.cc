#include "wall.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/** Doubles the memory allocation for shear records. */
void wall_list::add_memory() {

    // Check that the maximum limit on shear records hasn't been reached
    if(mem>=wall_list_max_memory) {
        fputs("Maximum wall list memory limit reached\n",stderr);
        exit(1);
    }

    // Allocate a new array with double the size
    mem<<=1;
    wall **na=new wall*[mem];

    // Copy the information from the old array into the new one
    for(int i=0;i<n;i++) na[i]=a[i];

    // Delete the old array and update the pointer
    delete [] a;
    a=na;
}

bool wall_par_planes::sep(double rx,double ry,double rz,double &dx,double &dy,double &dz,double rad,double wall_sca) {
    double s=rx*nx+ry*ny+rz*nz;
    double w=s-(s>0?d:-d);
    if(fabs(w)>rad*wall_sca) return false;
    dx=nx*w;
    dy=ny*w;
    dz=nz*w;
    return true;
}

bool wall_cylinder::sep(double rx,double ry,double rz,double &dx,double &dy,double &dz,double rad,double wall_sca) {
    double rr=rx*rx+ry*ry,drad=crad*wall_sca-rad;
    if(rr<drad*drad) return false;
    rr=sqrt(rr);
    dx=rx*(1-drad/rr);
    dy=ry*(1-drad/rr);
    dz=0;
    return true;
}

bool wall_floor::sep(double rx,double ry,double rz,double &dx,double &dy,double &dz,double rad,double wall_sca) {
    double ez=rz-fz;
    if(fabs(ez)>rad) return false;
    dx=dy=0;
    dz=ez;
    return true;
}
