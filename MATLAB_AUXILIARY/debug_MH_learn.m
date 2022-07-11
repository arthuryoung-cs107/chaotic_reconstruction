clear
close all
run AYfigprops.m

nbeads=5;
par_id=0;
MH_id=0;

mhgen = MH_genetic_data(nbeads, par_id, MH_id);
mhgen.event_block0 = MH_genetic_data.read_event_block_i(mhgen,0);
mhgen.Class0 = MH_genetic_data.read_Class_i(mhgen,0);
mhgen.gen0 = MH_genetic_data.read_gen_i(mhgen,0);
