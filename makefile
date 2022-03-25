include config_directory/makefile.includes

AY_DIR:=AYlinalg/objs/
AY_SRC:=AYlinalg/src/
SRC_DIR:=src/
OBJ_DIR:=objs/
TEST_DIR:=tests/

AYOBJS:=$(addprefix $(AY_DIR), $(addsuffix .o, $(AYLINALG)))
PTOBJS:=$(AYOBJS) $(addprefix $(OBJ_DIR), $(addsuffix .o, $(RYOBJS)))
STOBJS:=$(PTOBJS) $(addprefix $(OBJ_DIR), $(addsuffix .o, $(SWOBJS)))
FTOBJS:=$(STOBJS) $(addprefix $(OBJ_DIR), $(addsuffix .o, $(FIOBJS)))
RTOBJS:=$(STOBJS) $(addprefix $(OBJ_DIR), $(addsuffix .o, $(RAOBJS)))

$(AY_DIR)%.o: $(AY_SRC)%.c
	$(CC) $(IDIR) $(CFLAGS) -c $< -o $@

$(AY_DIR)%.o: $(AY_SRC)%.cc
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)%.o: $(SRC_DIR)%.cc
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

all: process_test swirl_test filter_test stat_test race_test

process_test: $(TEST_DIR)process_test.cc $(PTOBJS) | $(OBJ_DIR)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

swirl_test: $(TEST_DIR)process_test.cc $(STOBJS) | $(OBJ_DIR)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

filter_test: $(TEST_DIR)filter_test.cc $(FTOBJS) | $(OBJ_DIR)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

stat_test: $(TEST_DIR)stat_test.cc $(STOBJS) | $(OBJ_DIR)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

race_test: $(TEST_DIR)race_test.cc $(RTOBJS) | $(OBJ_DIR)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

clean:
	rm -f *_test
	rm -f ./objs/*.o

clean_AY:
	rm -f $(AY_DIR)*.o

clean_all: clean clean_AY
