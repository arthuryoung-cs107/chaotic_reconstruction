clear
close all
run AYfigprops.m

nbeads=5;
par_id=0;
MH_id=0;

mhgen = MH_genetic_data(nbeads, par_id, MH_id);
mhgen.event_block0 = mhgen.read_event_block_i(0);
mhgen.Class0 = mhgen.read_Class_i(0);
mhgen.gen0 = mhgen.read_gen_i(0);
