CXX=g++-11
CC=gcc-11
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

IDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

AYOBJS = AYaux.o AYlinalg.o AYmat.o AYvec.o AYrng.o AYtens.o AYdecomps.o AYdata.o

RYOBJS =  RYdat2AYdat.o

SWOBJS = p_grid.o swirl.o particle.o common.o swirl_io.o wall.o swirl_param.o dat_store.o

FIOBJS = filter_io.o filter.o

PTOBJS = $(RYOBJS) $(AYOBJS)

STOBJS = $(SWOBJS) $(RYOBJS) $(AYOBJS)

FTOBJS = $(FIOBJS) $(SWOBJS) $(RYOBJS) $(AYOBJS)

%.o:%.c
	$(CC) $(CFLAGS) $(IDIR) -c $<

%.o:%.cc
	$(CXX) $(CFLAGS) $(IDIR) -c $<

all: process_test swirl_test filter_test stat_test

process_test: process_test.cc $(PTOBJS)
	$(CXX) $(CFLAGS) $(IDIR) -o $@ $^ $(LINK)

swirl_test: swirl_test.cc $(STOBJS)
	$(CXX) $(CFLAGS) $(IDIR) -o $@ $^ $(LINK)

filter_test: filter_test.cc $(FTOBJS)
	$(CXX) $(CFLAGS) $(IDIR) -o $@ $^ $(LINK)

stat_test: stat_test.cc $(STOBJS)
	$(CXX) $(CFLAGS) $(IDIR) -o $@ $^ $(LINK)

clean:
	rm -f process_test
	rm -f swirl_test
	rm -f filter_test
	rm -f stat_test
	rm -f *.o
clean_dat:
	rm -f *.aysml
	rm -f *.rysml
	rm -f *.aydat
clean_datdir:
	rm -f ./dat_dir/*.aysml
	rm -f ./dat_dir/*.rysml
	rm -f ./dat_dir/*.aydat

clean_all: clean clean_dat
