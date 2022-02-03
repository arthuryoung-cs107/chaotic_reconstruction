#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <cmath>

const int particle_init_shear_records=4;
const int particle_max_shear_records=32;

class particle {
    public:
        /** The particle position. */
        double x,y,z;
        /** The particle velocity. */
        double vx,vy,vz;
        /** The current particle acceleration. */
        double ax,ay,az;
        /** The new particle acceleration. */
        double nax,nay,naz;
        /** The quaternion describing the rotation of the particle. */
        double q0,q1,q2,q3;
        /** The particle angular velocity. */
        double omegax,omegay,omegaz;
        /** The current particle torque. */
        double tx,ty,tz;
        /** The new particle torque. */
        double ntx,nty,ntz;
        particle() : co(0), mem(particle_init_shear_records),
            sj(new int[mem]), sr(new double[3*mem]) {}
        ~particle() {
            delete [] sr;
            delete [] sj;
        }
        void verlet(double dt);
        void set_gravity();
        void copy(particle &p);
        int find_shear(int j);
        void purge_shear();
        inline void zero() {
            x=y=z=q0=q1=q2=q3=0;
        }
        inline void add(particle &p,double gw) {
            x+=gw*p.x;
            y+=gw*p.y;
            z+=gw*p.z;
            q0+=gw*p.q0;
            q1+=gw*p.q1;
            q2+=gw*p.q2;
            q3+=gw*p.q3;
        }
        inline void scale(double cw) {
            x*=cw;y*=cw;z*=cw;
            normalize_q();
        }
        inline void set_pos(double x_,double y_,double z_) {
            x=x_;
            y=y_;
            z=z_;clear_q();
        }
        inline void tweak_pos(double delx,double dely,double idur) {
            x+=delx;y+=dely;
            vx+=delx*idur;vy+=dely*idur;
        }
        inline void clear_q() {
            q0=1;
            q1=q2=q3=0;
        }
        inline void normalize_q() {
            double qmod=1./sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
            q0*=qmod;q1*=qmod;q2*=qmod;q3*=qmod;
        }
        void add_shear_memory();
        /** The shear record mask. */
        unsigned int mask;
        /** The total number of shear records. */
        int co;
        /** The current memory allocation of shear records. */
        int mem;
        /** The shear record indices. */
        int *sj;
        /** The shear records. */
        double *sr;
};

#endif
