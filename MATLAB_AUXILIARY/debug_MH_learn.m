clear
close all
run AYfigprops.m

pos_full = [1 1 1728 1000];

full_fig_specs = AYfig.specs_gen('full_fig', pos_full);

nbeads=5;
par_id=0;
MH_id=0;

stat = read_stat(3,'maxmin',0);
swtrue = stat.sw0;
utrue_full=swtrue.params;
utrue=utrue_full(3:end);

mhgen = MH_genetic_data(nbeads, par_id, MH_id);
mhgen.event_block0 = MH_genetic_data.read_event_block_i(mhgen,0);

gen0 = MH_genetic_data.read_gen_i(mhgen,0);
gen1 = MH_genetic_data.read_gen_i(mhgen,1);
gen2 = MH_genetic_data.read_gen_i(mhgen,2);
gen3 = MH_genetic_data.read_gen_i(mhgen,3);
gen4 = MH_genetic_data.read_gen_i(mhgen,4);
gen5 = MH_genetic_data.read_gen_i(mhgen,5);
gen6 = MH_genetic_data.read_gen_i(mhgen,6);

Class0 = MH_genetic_data.read_Class_i(mhgen,0);
Class1 = MH_genetic_data.read_Class_i(mhgen,1);
Class2 = MH_genetic_data.read_Class_i(mhgen,2);
Class3 = MH_genetic_data.read_Class_i(mhgen,3);
Class4 = MH_genetic_data.read_Class_i(mhgen,4);
Class5 = MH_genetic_data.read_Class_i(mhgen,5);

Class10 = MH_genetic_data.read_Class_i(mhgen,10);
Class20 = MH_genetic_data.read_Class_i(mhgen,20);
Class30 = MH_genetic_data.read_Class_i(mhgen,30);
Class39 = MH_genetic_data.read_Class_i(mhgen,39);

mhgen.gen0 = gen0;
mhgen.Class0 = Class0;

% pool_list=[gen0.pool gen1.pool gen2.pool gen3.pool gen4.pool gen5.pool];
% leaders_list=[Class0.leaders Class1.leaders Class2.leaders Class3.leaders Class4.leaders Class5.leaders];
leaders_list=[Class0.leaders Class5.leaders Class10.leaders Class20.leaders Class30.leaders Class39.leaders];

full_fig=AYfig(full_fig_specs,false);
full_fig.init_tiles([2,3]);
% debug_MH_plots.plot_rec_uerr(full_fig, pool_list, utrue, blue5)
debug_MH_plots.plot_rec_uerr(full_fig, leaders_list, utrue, green4)
