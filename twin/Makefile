# Include the configuration file
include ../../config.mk

# Create iflag to include files in the shared folder
iflags=-I../shared `gsl-config --cflags`
lflags=`gsl-config --libs`

# Lists of files to be built
objs=p_grid.o swirl.o particle.o common.o filter_io.o swirl_io.o wall.o filter.o swirl_param.o dat_store.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=swirl_test fi_test

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

swirl_test: swirl_test.cc $(objs)
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

fi_test: fi_test.cc $(objs)
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

data:
	mkdir -p data
	rsync -r olympus.seas.harvard.edu:"/data/files/twin/iteration_0_output.npz.xz /data/files/twin/iteration_0_trajectories.npz.xz" data
	unxz -f data/iteration_0_output.npz.xz
	unxz -f data/iteration_0_trajectories.npz.xz

clean:
	rm -f $(objs) $(execs)

.PHONY: clean depend data
