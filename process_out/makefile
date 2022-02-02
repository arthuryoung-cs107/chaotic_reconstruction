CXX=g++-11
CC=gcc-11
OS=OSX
CFLAGS= -O3 -I. -fopenmp -DGSL -DKNUTH

IDIR= -I/usr/local/include/
LINK= -L/usr/local/lib/ -lm -lgsl -lgslcblas -lfftw3 -lblas

CDEPS = AYaux.h

CPPDEPS = AYlinalg.hh RYdat2AYdat.hh

OBJS1= AYaux.o AYlinalg.o AYmat.o AYvec.o AYrng.o AYtens.o AYdecomps.o RYdat2AYdat.o

%.o:%.c
	$(CC) $(CFLAGS) $(IDIR) -c $<

%.o:%.cc
	$(CXX) $(CFLAGS) $(IDIR) -c $<

all: test1

test1: main1.cc $(OBJS1)
	$(CXX) $(CFLAGS) $(IDIR) -o $@ $^ $(LINK)

clean:
	rm -f test1
	rm -f *.o
clean_dat:
	rm -f *.dat
	rm -f *.aysml
	rm -f *.aydat
clean_datdir:
	rm -f ./dat_dir/*.dat
	rm -f ./dat_dir/*.aysml
	rm -f ./dat_dir/*.aydat

clean_all: clean clean_dat
