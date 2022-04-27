include config_directory/makefile.includes

# AYlinalg object files
AYOBJS:= $(addprefix $(AY_DIR), $(addsuffix .o, $(AYLINALG)))

SWOBJS:= $(AYOBJS) $(addprefix $(SW_DIR), $(addsuffix .o, $(SWIRL)))

FIOBJS:= $(SWOBJS) $(addprefix $(FI_DIR), $(addsuffix .o, $(FILTER)))

RAOBJS:= $(SWOBJS) $(addprefix $(RA_DIR), $(addsuffix .o, $(RACE)))

WAOBJS:= $(SWOBJS) $(addprefix $(WA_DIR), $(addsuffix .o, $(WALK)))

REOBJS:= $(SWOBJS) $(addprefix $(RE_DIR), $(addsuffix .o, $(RELAY)))

# rules for each directory
# AYlinalg rules
$(AY_DIR)%.o: $(AY_SRC)%.c | $(AY_DIR)
	$(CC) $(IDIR) $(CFLAGS) -c $< -o $@

$(AY_DIR)%.o: $(AY_SRC)%.cc | $(AY_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

$(SW_DIR)%.o: $(SW_SRC)%.cc | $(SW_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

$(FI_DIR)%.o: $(FI_SRC)%.cc | $(FI_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

$(RA_DIR)%.o: $(RA_SRC)%.cc | $(RA_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

$(WA_DIR)%.o: $(WA_SRC)%.cc | $(WA_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

$(RE_DIR)%.o: $(RE_SRC)%.cc | $(RE_DIR)
	$(CXX) $(IDIR) $(CFLAGS) -c $< -o $@

all: process_test swirl_test filter_test stat_test race_test walk_test relay_test

process_test: $(TEST_SRC)process_test.cc $(SWOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

swirl_test: $(TEST_SRC)swirl_test.cc $(SWOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

filter_test: $(TEST_SRC)filter_test.cc $(FIOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

stat_test: $(TEST_SRC)stat_test.cc $(SWOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

race_test: $(TEST_SRC)race_test.cc $(RAOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

walk_test: $(TEST_SRC)walk_test.cc $(WAOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

relay_test: $(TEST_SRC)relay_test.cc $(REOBJS)
	$(CXX) $(IDIR) $(CFLAGS) $(LINK) $^ $(LIBS) -o $@

$(AY_DIR) $(SW_DIR) $(FI_DIR) $(RA_DIR) $(WA_DIR) $(RE_DIR):
	mkdir -p $@

clean_:
	rm -f *_test

clean_AY:
	rm -f $(AY_DIR)*.o

clean_swirl:
	rm -f $(SW_DIR)*.o

clean_filter:
	rm -f $(FI_DIR)*.o

clean_race:
	rm -f $(RA_DIR)*.o

clean_walk:
	rm -f $(WA_DIR)*.o

clean_relay:
	rm -f $(RE_DIR)*.o

clean_process: clean_swirl

clean: clean_ clean_swirl clean_filter clean_race clean_relay

clean_all: clean_ clean_AY clean_swirl clean_filter clean_race clean_relay
