#ifndef WALL_HH
#define WALL_HH

const int wall_list_init_memory=8;
const int wall_list_max_memory=128;

struct wall {
    virtual bool sep(double rx,double ry,double rz,double &dx,double &dy,double &dz,double rad,double wall_sca) = 0;
};

class wall_list {
    public:
        /** The total number of walls. */
        int n;
        /** The memory allocation for walls. */
        int mem;
        /** An array of pointers to walls. */
        wall **a;
        wall_list() : n(0), mem(wall_list_init_memory), a(new wall*[mem]) {}
        ~wall_list() {delete [] a;}
        inline void add_wall(wall *w) {
            if(n==mem) add_memory();
            a[n++]=w;
        }
    private:
        void add_memory();
};

struct wall_par_planes : public wall {
    /* The x component of the normal vector. */
    const double nx;
    /* The y component of the normal vector. */
    const double ny;
    /* The z component of the normal vector. */
    const double nz;
    /* The separation of the planes along the normal vector. */
    const double d;
    wall_par_planes(double nx_,double ny_,double nz_,double d_)
        : nx(nx_), ny(ny_), nz(nz_), d(d_) {}
    virtual bool sep(double rx,double ry,double rz,double &dx,double &dy,double &dz,double rad,double wall_sca);
};

struct wall_cylinder : public wall {
    /** The radius of the cylinder. */
    const double crad;
    wall_cylinder(double crad_) : crad(crad_) {}
    virtual bool sep(double rx,double ry,double rz,double &dx,double &dy,double &dz,double rad,double wall_sca);
};

struct wall_floor : public wall {
    /** The z displacement of the floor. */
    const double fz;
    wall_floor(double fz_) : fz(fz_) {}
    virtual bool sep(double rx,double ry,double rz,double &dx,double &dy,double &dz,double rad,double wall_sca);
};

#endif
