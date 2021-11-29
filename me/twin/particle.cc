#include "particle.hh"

#include <cstdio>
#include <cstdlib>

/** Doubles the memory allocation for shear records. */
void particle::add_shear_memory() {

    // Check that the maximum limit on shear records hasn't been reached
    if(mem>=particle_max_shear_records) {
        fputs("Maximum shear records reached\n",stderr);
        exit(1);
    }

    // Allocate the new arrays with double the size
    mem<<=1;
    int *nsj=new int[mem];
    double *nsr=new double[3*mem];

    // Copy the information from the old arrays into the new ones
    for(int i=0;i<co;i++) {
        nsj[i]=sj[i];
        nsr[3*i]=sr[3*i];
        nsr[3*i+1]=sr[3*i+1];
        nsr[3*i+2]=sr[3*i+2];
    }

    // Delete the old arrays and update the pointers
    delete [] sr;
    delete [] sj;
    sj=nsj;
    sr=nsr;
}

void particle::copy(particle &p) {
    x=p.x;y=p.y;z=p.z;
    vx=p.vx;vy=p.vy;vz=p.vz;
    ax=p.ax;ay=p.ay;az=p.az;
    q0=p.q0;q1=p.q1;q2=p.q2;q3=p.q3;
    omegax=p.omegax;
    omegay=p.omegay;
    omegaz=p.omegaz;
    tx=p.tx;ty=p.ty;tz=p.tz;
    mask=p.mask;
    co=p.co;
    while(co>mem) add_shear_memory();
    for(int i=0;i<co;i++) {
        sj[i]=p.sj[i];
        sr[3*i]=p.sr[3*i];
        sr[3*i+1]=p.sr[3*i+1];
        sr[3*i+2]=p.sr[3*i+2];
    }
}

void particle::verlet(double dt) {
    double hdt=0.5*dt;
    x+=dt*(vx+hdt*ax);
    y+=dt*(vy+hdt*ay);
    z+=dt*(vz+hdt*az);
    vx+=hdt*(ax+nax);
    vy+=hdt*(ay+nay);
    vz+=hdt*(az+naz);
    ax=nax;
    ay=nay;
    az=naz;
    double dq0=-q1*omegax-q2*omegay-q3*omegaz,
           dq1=q0*omegax+q3*omegay-q2*omegaz,
           dq2=-q3*omegax+q0*omegay+q1*omegaz,
           dq3=q2*omegax-q1*omegay+q0*omegaz;
    q0+=hdt*dq0;
    q1+=hdt*dq1;
    q2+=hdt*dq2;
    q3+=hdt*dq3;
    omegax+=hdt*(tx+ntx);
    omegay+=hdt*(ty+nty);
    omegaz+=hdt*(tz+ntz);
    tx=ntx;
    ty=nty;
    tz=ntz;
}

void particle::set_gravity() {
    nax=0;nay=0;naz=-1;
    ntx=0;nty=0;ntz=0;
    mask=0;
}

int particle::find_shear(int j) {
    int i=0;
    for(;i<co;i++) if(sj[i]==j) break;
    if(i==co) {
        if(co==mem) add_shear_memory();
        co++;
        sj[i]=j;
        sr[3*i]=sr[3*i+1]=sr[3*i+2]=0;
    }
    mask|=1<<i;
    return i;
}

void particle::purge_shear() {
    int i=0;
    while(i<co) {
        if((mask&(1<<i))==0) {
            sj[i]=sj[--co];
            sr[3*i]=sr[3*co];
            sr[3*i+1]=sr[3*co+1];
            sr[3*i+2]=sr[3*co+2];
            mask|=(mask&(1<<co))>>(co-i);
        } else i++;
    }
}
