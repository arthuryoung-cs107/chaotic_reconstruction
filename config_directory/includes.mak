CXX:=g++-11
CC:=gcc-11
OS:=OSX
CFLAGS:= -march=native -O3 -fopenmp -DGSL -DKNUTH

IBASIC:=-I/usr/local/include/
LINK:=-L/usr/local/lib/
LIBS:=-lm -lgsl -lgslcblas -lfftw3 -lblas

# directories
AY_SRC:= AYlinalg/linalg/
AY_DIR:= AYlinalg/linalg_objs/

SW_SRC:= swirl/
FI_SRC:= filter/
RA_SRC:= race/
WA_SRC:= walk/
RE_SRC:= relay/
MH_SRC:= Metropolis_Hastings/

SW_DIR:=$(addsuffix objs/, $(SW_SRC))
FI_DIR:=$(addsuffix objs/, $(FI_SRC))
RA_DIR:=$(addsuffix objs/, $(RA_SRC))
WA_DIR:=$(addsuffix objs/, $(WA_SRC))
RE_DIR:=$(addsuffix objs/, $(RE_SRC))
MH_DIR:=$(addsuffix objs/, $(MH_SRC))

TEST_SRC:=tests/

WRITING_SRC:=writing_tests/

IDIR:=$(IBASIC) -I$(AY_SRC) -I$(SW_SRC) -I$(FI_SRC) -I$(RA_SRC) -I$(WA_SRC) -I$(RE_SRC) -I$(MH_SRC)
