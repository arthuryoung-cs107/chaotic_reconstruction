#include "swirl.hh"

#include <cstdio>
#include <cmath>

#include <gsl/gsl_rng.h>

/** Initializes the particle swirling simulation.
 * \param[in] sp the physical parameters.
 * \param[in] pg_ a pointer to the proximity grid to use.
 * \param[in] wl_ the list of walls.
 * \param[in] n_ the total number of particles. */
swirl::swirl(swirl_param &sp,proximity_grid *pg_,wall_list &wl_,int n_) : swirl_param(sp), n(n_),
    q(new particle[n]), time(0.), ctheta(0.), comega(0.), pg(pg_), dstore(NULL),
    wl(wl_), odir(NULL) {
    set_dish_coords();
}

/** Initializes the particles swirling simulation as a weighted average of an
 * array of others.
 * \param[in] sw an array of pointers to swirling simulations.
 * \param[in] npar the size of the sw array.
 * \param[in] gamma the gamma weighting to apply to the simulation weights. */
swirl::swirl(swirl **sw,int npar,double gamma) : swirl_param(*(*sw),0.),
    n((*sw)->n), q(new particle[n]), time(0.), ctheta(0.), comega(0.),
    pg(NULL), wl((*sw)->wl), odir((*sw)->odir) {

    // Clear all bead information
    for(int i=0;i<n;i++) q[i].zero();

    // Loop over all of the swirling simulations
    double cw=0.,gw,cth=0,sth=0;
    for(int k=0;k<npar;k++) {
        swirl *swp=sw[k];

        // Compute gamma-transformed weight
        cw+=gw=pow(swp->wei,gamma);

        // Add the contribution from this swirling simulation
        add_param(*swp,gw);
        for(int i=0;i<n;i++) q[i].add(swp->q[i],gw);
        time+=gw*swp->time;
        cth+=gw*cos(swp->ctheta);
        sth+=gw*sin(swp->ctheta);
        comega+=gw*swp->comega;
    }

    // Normalize the data
    cw=1./cw;
    scale_param(cw);
    for(int i=0;i<n;i++) q[i].scale(cw);
    time*=cw;ctheta=atan2(sth,cth);comega*=cw;
    set_dish_coords();
}

/** The class destructor frees the dynamically allocated memory. */
swirl::~swirl() {
    if(odir!=NULL) delete [] odir;
    delete [] q;
}

/** Copies the swirling simulation from another, assuming the same number of
 * particles.
 * \param[in] sw a reference to another swirling simulation. */
void swirl::copy(swirl &sw) {
    copy_param(sw);
    ctheta=sw.ctheta;
    comega=sw.comega;
    set_dish_coords();
    for(int i=0;i<n;i++) q[i].copy(sw.q[i]);
}

/** Simulates the swirling particle system, periodically saving output.
 * \param[in] dur the duration to simulate over.
 * \param[in] dt the timestep to use.
 * \param[in] frames the number of frames to output. */
void swirl::solve(double dur,double dt,int frames) {

    // Physical constants
    double g_phys=9.804,                 // Gravity (m/s^2)
           d_phys=0.00635,               // Diameter (m)
           t_phys=sqrt(d_phys/g_phys),   // Time unit (s)
           comega_max=36*t_phys;

    double fint=dur/frames;
    int l=static_cast<int>(fint/dt)+1;
    double adt=fint/l;
    printf("# %d steps per frame\n",l);

    puts("# Output frame 0");
    output(0);

    for(int j=1;j<=frames;j++) {

        for(int k=0;k<l;k++) {
            step_forward(adt);
            comega+=0.015*adt;if(comega>comega_max) comega=comega_max;
        }

        printf("# Output frame %d (%g)\n",j,comega);
        output(j);
    }
}

void swirl::advance(double dur,double ctheta_,double comega_,double dt) {
    set_swirl(ctheta_,comega_);
    int l=static_cast<int>(dur/dt)+1,k;
    double adt=dur/l;
    for(k=0;k<l;k++) step_forward(adt);
}

/** Steps the particle positions forward. */
void swirl::step_forward(double dt) {

    // Advance the time and update the dish center position
    time+=dt;
    update_swirl(dt);

    // Set up the proximity grid, periodically scanning the particle positions
    pg->setup(q,n);

    // Add the particles to the proximity grid
    pg->clear();
    for(int i=0;i<n;i++) {
        pg->put(i,q[i].x,q[i].y,q[i].z);
        q[i].set_gravity();
    }
    pg->resolve_overflows();

    // Search for all inter-particle contacts
    int i,j,k,ijk,l,*co=pg->co;
    for(ijk=k=0;k<pg->o;k++) for(j=0;j<pg->n;j++) for(i=0;i<pg->m;i++,ijk++)
        for(l=0;l<co[ijk];l++) {
            point_info &p=pg->p[ijk][l];
            int id=p.id,li,ui,lj,uj,lk,uk;

            // Find all blocks the could have a sphere in contact with the current one
            pg->subregion(p.x,p.y,p.z,2*rad,li,ui,lj,uj,lk,uk);
            for(int ck=lk;ck<=uk;ck++) for(int cj=lj;cj<=uj;cj++)
                for(int ci=li;ci<=ui;ci++) {
                int ijk=ci+pg->m*(cj+pg->n*ck);
                point_info *pip=pg->p[ijk],*pie=pip+pg->co[ijk];
                for(;pip<pie;pip++) if(pip->id<id) {
                    double delx=p.x-pip->x,
                           dely=p.y-pip->y,
                           delz=p.z-pip->z,
                           rsq=delx*delx+dely*dely+delz*delz;
                    if(rsq<4*rad*rad) contact(id,pip->id,delx,dely,delz,rsq,dt);
                }
            }
    }

    // Handle all contacts between the particles and the walls
    double dx,dy,dz;
    for(i=0;i<n;i++) for(l=0;l<wl.n;l++)
        if(wl.a[l]->sep(q[i].x-cx,q[i].y-cy,q[i].z-cz,dx,dy,dz,rad,wall_sca))
            contact(i,-1-l,dx,dy,dz,dx*dx+dy*dy+dz*dz,dt);

    // Integrate the particle position, velocity, and spin. Purge any shear
    // records for any contacts that have been lost.
    for(int i=0;i<n;i++) {
        q[i].purge_shear();
        q[i].verlet(dt);
    }
}

/** Computes the contact force between a particle and another particle or wall.
 * \param[in] id1 the ID of the particle to consider.
 * \param[in] id2 the ID of the second particle to consider; if negative, then
 *                the routine does a wall contact.
 * \param[in] (delx,dely,delz) a minimum separation vector between.
 * \param[in] dt the simulation timestep, used to update the shear information.
 */
void swirl::contact(int id1,int id2,double delx,double dely,double delz,double rsq,double dt) {
    bool w=id2<0;
    particle &q1=q[id1],*qp2=w?NULL:q+id2;
    double r=sqrt(rsq),rinv=1./r,rsqinv=1./rsq,

    // Relative translational velocity
           vr1=q1.vx-(w?cvx:qp2->vx),
           vr2=q1.vy-(w?cvy:qp2->vy),
           vr3=q1.vz-(w?cvz:qp2->vz),

    // Normal component
           vnnr=vr1*delx+vr2*dely+vr3*delz,
           vn1=delx*vnnr*rsqinv,
           vn2=dely*vnnr*rsqinv,
           vn3=delz*vnnr*rsqinv,

    // Tangential component
           vt1=vr1-vn1,vt2=vr2-vn2,vt3=vr3-vn3,

    // Relative rotational velocity
           wr1=(rad*q1.omegax+(w?0:rad*qp2->omegax))*rinv,
           wr2=(rad*q1.omegay+(w?0:rad*qp2->omegay))*rinv,
           wr3=(rad*q1.omegaz+(w?0:rad*qp2->omegaz))*rinv,

    // Effective mass of pair of particles
           meff=w?mass:0.5*mass,
           gamman=w?(id2==-1?gamman_base:gamman_wall):gamman_bead,
           gammat=0.5*gamman,

    // Normal forces = Hookean contact + normal velocity damping
           damp=meff*gamman*vnnr*rsqinv,
           ccel=Kn*((w?rad:2*rad)-r)*rinv-damp,

    // Relative velocities
           vtr1=vt1-(delz*wr2-dely*wr3),
           vtr2=vt2-(delx*wr3-delz*wr1),
           vtr3=vt3-(dely*wr1-delx*wr2);

    // Locate the shear information
    int j=q1.find_shear(id2);
    double *sr=q1.sr+3*j;

    // Shear history effects
    *sr+=vtr1*dt;
    sr[1]+=vtr2*dt;
    sr[2]+=vtr3*dt;
    double shrmag=sqrt(sr[0]*sr[0]+sr[1]*sr[1]+sr[2]*sr[2]),

    // Rotate shear displacements
           rsht=(sr[0]*delx+sr[1]*dely+sr[2]*delz)*rsqinv;
    *sr-=rsht*delx;
    sr[1]-=rsht*dely;
    sr[2]-=rsht*delz;

    // Get relevant friction coefficient
    double xmu=w?(id2==-1?mu_base:mu_wall):mu_bead,
           Kt=2/7.*Kn,

    // Tangential forces = Shear + Tangential velocity damping
           fs1=-Kt*sr[0]-meff*gammat*vtr1,
           fs2=-Kt*sr[1]-meff*gammat*vtr2,
           fs3=-Kt*sr[2]-meff*gammat*vtr3,

    // Rescale the frictional displacements and forces if a Coulomb criterion
    // is reached
           fs=sqrt(fs1*fs1+fs2*fs2+fs3*fs3),
           fn=xmu*fabs(ccel*r);
    if(fs>fn) {
        if(shrmag!=0.) {
            double sf=meff*gammat/Kt;
            sr[0]=(fn/fs)*(sr[0]+sf*vtr1)-sf*vtr1,
            sr[1]=(fn/fs)*(sr[1]+sf*vtr2)-sf*vtr2,
            sr[2]=(fn/fs)*(sr[2]+sf*vtr3)-sf*vtr3;
            fs1*=fn/fs;
            fs2*=fn/fs;
            fs3*=fn/fs;
        } else fs1=fs2=fs3=0;
    }

    // Apply forces
    double fx=delx*ccel+fs1,
           fy=dely*ccel+fs2,
           fz=delz*ccel+fs3;
    q1.nax+=fx/mass;
    q1.nay+=fy/mass;
    q1.naz+=fz/mass;

    // Apply torques
    double tor1=rinv*(dely*fs3-delz*fs2),
           tor2=rinv*(delz*fs1-delx*fs3),
           tor3=rinv*(delx*fs2-dely*fs1);
    q1.ntx-=rad*tor1/momi;
    q1.nty-=rad*tor2/momi;
    q1.ntz-=rad*tor3/momi;

    // If this is a particle-particle contact, then apply equal and opposite
    // forces to the other particle
    if(!w) {
        qp2->nax-=fx/mass;
        qp2->nay-=fy/mass;
        qp2->naz-=fz/mass;
        qp2->ntx-=rad*tor1/momi;
        qp2->nty-=rad*tor2/momi;
        qp2->ntz-=rad*tor3/momi;
    }
}

/** Initializes the positions of the beads to match an experimental snapshot
 * with a small displacement.
 * \param[in] time_ the simulation time.
 * \param[in] ctheta_ the dish angle.
 * \param[in] f a pointer to the experimental data, consisting of (x,y)
 *              positions measured in pixels.
 * \param[in] r a pointer to the GSL random number generator to use. */
void swirl::init_positions(double time_,double ctheta_,float *f,gsl_rng *r) {

    // Initialize time and the dish parameters
    time=time_;
    set_swirl(ctheta_,0);

    // Set the initial positions of the beads to match the experimental data
    // with some random perturbations
    double icl=1./cl_im;
    for(int i=0;i<n;i++,f+=2)
        q[i].set_pos((*f+(2*gsl_rng_uniform(r)-1)-cx_im)*icl+cx,
                     (f[1]+(2*gsl_rng_uniform(r)-1)-cy_im)*icl+cy,rad);
}

/** Updates the weight of the simulation depending on how well it matches the
 * experimental data.
 * \param[in] f a pointer to the experimental data, consisting of (x,y)
 *              positions measured in pixels.
 * \param[in] dur the duration between the current calculation and the previous one.
 * \param[in] t_wheels the training wheels parameter, setting the amount of
 *                     drift toward the true positions.
 * \param[in] lnorm an array of three doubles for computing distance
 *                  diagnostics in the L_1, L_2, and L_infinity norms (without
 *                  normalization). */
void swirl::update_weight(float *f,double dur,double t_wheels,double *lnorm) {
    double fac=t_wheels/cl_im,idur=1./dur,rmax=cl_im*cl_im;
    *lnorm=lnorm[1]=lnorm[2]=0;logfac=0;
    for(int i=0;i<n;i++,f+=2) {

        // Calculate the difference between the simulation position and the
        // experimental data
        double xt=(q[i].x-cx)*cl_im+cx_im-*f,
               yt=(q[i].y-cy)*cl_im+cy_im-f[1],
               rsq=xt*xt+yt*yt,r=sqrt(rsq);

        // Compute the diagnostic measurements
        *lnorm+=r;
        lnorm[1]+=rsq;
        if(r>lnorm[2]) lnorm[2]=r;
        logfac=rsq>rmax?rmax:rsq;

        // Adjust the position of bead to drift toward the experimental
        // position
        q[i].tweak_pos(-fac*xt,-fac*yt,idur);
    }
    logfac=lnorm[1];
}

/** Applies random adjustments to the simulation parameters and particle
 * velocities and spins.
 * \param[in] sp_min the minimum values for the swirling parameters.
 * \param[in] sp_max the maximum values for the swirling parameters.
 * \param[in] sp_rnd the random perturbations to apply to the swirling
 *                   parameters.
 * \param[in] r a pointer to the GSL random number generator to use.
 * \param[in] v_pert the size of the velocity perturbation to apply in each
 *                   direction.
 * \param[in] w_pert the size of the angular velocity perturbation to apply in
 *                   each direction. */
void swirl::jiggle(swirl_param &sp_min,swirl_param &sp_max,swirl_param &sp_rnd,gsl_rng *r,double v_pert,double omega_pert) {
    jiggle_param(sp_min,sp_max,sp_rnd,r);

    // Apply random displacements to the bead velocities and angular velocities
    for(int i=0;i<n;i++) {
        q[i].vx+=v_pert*(2*gsl_rng_uniform(r)-1);
        q[i].vy+=v_pert*(2*gsl_rng_uniform(r)-1);
        q[i].vz+=v_pert*(2*gsl_rng_uniform(r)-1);
        q[i].omegax+=omega_pert*(2*gsl_rng_uniform(r)-1);
        q[i].omegay+=omega_pert*(2*gsl_rng_uniform(r)-1);
        q[i].omegaz+=omega_pert*(2*gsl_rng_uniform(r)-1);
    }
}

/** Sets the position and velocity of dish from its current angle and angular
 * velocity. The routine also remaps the dish angle to the range from 0 to
 * 2*pi.*/
void swirl::set_dish_coords() {

    // Remap the dish angle from 0 to 2*pi
    if(ctheta<0) ctheta+=2*M_PI;
    else if(ctheta>=2*M_PI) ctheta-=2*M_PI;

    // Set the dish position
    double cc=cos(ctheta),ss=sin(ctheta);
    cx=dsamp*cc;
    cy=dsamp*ss;
    cz=0;

    // Set the dish velocity
    cvx=-dsamp*comega*ss;
    cvy=dsamp*comega*cc;
    cvz=0;
}
