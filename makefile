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

.PHONY: all clean

all: process_test swirl_test filter_test stat_test race_test

$(AY_DIR)%.o: $(AY_SRC)%.c | $(AY_DIR)
	$(CC) $(IDIR) $(CFLAGS) -c $< -o $@

$(AY_DIR)%.o: $(AY_SRC)%.cc | $(AY_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)%.o: $(SRC_DIR)%.cc | $(OBJ_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

process_test: $(TEST_DIR)process_test.cc $(PTOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

swirl_test: $(TEST_DIR)swirl_test.cc $(STOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

filter_test: $(TEST_DIR)filter_test.cc $(FTOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

stat_test: $(TEST_DIR)stat_test.cc $(STOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

race_test: $(TEST_DIR)race_test.cc $(RTOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

$(AY_DIR) $(OBJ_DIR):
	mkdir -p $@

clean_:
	rm -f *_test
	rm -f ./objs/*.o

clean_AY:
	rm -f $(AY_DIR)*.o

clean: clean_ clean_AY
