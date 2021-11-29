#ifndef P_GRID_HH
#define P_GRID_HH

#include "pg_structs.hh"

class particle;

/** The default memory to allocate for each grid block. */
const int pgrid_default_memory=4;

/** The maximum memory to allocate in a grid block (or overflow buffer) before
 * creating an error. */
const int pgrid_max_memory=16777216;

/** The average number of points to aim for in each grid block, for optimal
 * performance. */
const double pgrid_optimal_points=2.;

/** The number of calculations to perform before trying to rebuild the grid
 * data structure. */
const int pgrid_build_interval=8192;

/** The initial size of the point overflow buffer. */
const int pgrid_init_overflow_buffer=64;

class proximity_grid {
    public:
        /** The number of grid blocks in the x direction. */
        int m;
        /** The number of grid blocks in the y direction. */
        int n;
        /** The number of grid blocks in the z direction. */
        int o;
        /** The total number of grid blocks. */
        int mno;
        /** A counter for the number of calculations that have been performed.
         */
        int count;
        /** The lower x coordinate of the grid structure. */
        double ax;
        /** The lower y coordinate of the grid structure. */
        double ay;
        /** The lower z coordinate of the grid structure. */
        double az;
        /** The size of a block in the x direction. */
        double dx;
        /** The size of a block in the y direction. */
        double dy;
        /** The size of a block in the z direction. */
        double dz;
        /** The reciprocal size of a block in the x direction. */
        double xsp;
        /** The reciprocal size of a block in the y direction. */
        double ysp;
        /** The reciprocal size of a block in the z direction. */
        double zsp;
        /** An array of the number of points in each block. */
        int* co;
        /** An array of the amount of memory allocated to each block. */
        int* mem;
        /** An array of particle IDs and positions in each block. */
        point_info** p;
        proximity_grid();
        ~proximity_grid();
        void setup(particle *q,int c);
        void resolve_overflows();
        void allocate(particle *q,int c,int imem=pgrid_default_memory,double opt_points=pgrid_optimal_points);
        void put(int id,double x,double y,double z);
        void free_grid();
        inline void subregion(double x,double y,double z,double r,int &li,int &ui,int &lj,int &uj,int &lk,int &uk) {
            li=(x-ax-r)*xsp,lj=(y-ay-r)*ysp,lk=(z-az-r)*zsp;
            ui=(x-ax+r)*xsp,uj=(y-ay+r)*ysp,uk=(z-az+r)*zsp;
            map_inside(li,lj,lk);
            map_inside(ui,uj,uk);
        }
        inline void clear() {
            co_ov=0;
#pragma omp parallel for
            for(int *cop=co;cop<co+mno;cop++) *cop=0;
        }
    private:
        inline void map_inside(int &i,int &j,int &k) {
            if(i<0) i=0;else if(i>m-1) i=m-1;
            if(j<0) j=0;else if(j>n-1) j=n-1;
            if(k<0) k=0;else if(k>o-1) k=o-1;
        }
        void add_memory(int ijk);
        void add_overflow_memory();
        /** The total number of points in the overflow buffer. */
        int co_ov;
        /** The size of the overflow buffer. */
        int mem_ov;
        /** The overflow buffer. */
        oflow_info *ov;
};

#endif
