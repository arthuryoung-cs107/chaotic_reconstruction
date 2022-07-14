clear
close all
run AYfigprops.m

pos_full = [1 1 1728 1000];

full_fig_specs = AYfig.specs_gen('full_fig', pos_full);

nbeads=3;
par_id=0;
MH_id=0;
e_last = 1;
c_last = 32;
g_last=112;

% nbeads=3;
% par_id=0;
% MH_id=0;
% e_last = 1;
% c_last = 32;
% g_last=112;


stat = read_stat(3,'maxmin',0);
swtrue = stat.sw0;
utrue_full=swtrue.params;
utrue=utrue_full(3:end);

mhgen = MH_genetic_data(nbeads, par_id, MH_id);
mhgen.event_block0 = MH_genetic_data.read_event_block_i(mhgen,0);
mhgen.gen0 = MH_genetic_data.read_gen_i(mhgen,0);
mhgen.Class0 = MH_genetic_data.read_Class_i(mhgen,0);

full_fig=AYfig(full_fig_specs,false);
full_fig.init_tiles([3,4]);

for i = 1:c_last
    i
    Class_i = MH_genetic_data.read_Class_i(mhgen,i);
    MH_genetic_data.print_Class(Class_i);
    debug_MH_plots.plot_Class_diagnostics(full_fig,mhgen,Class_i,utrue,green4)

    pause
end


% debug_MH_plots.plot_reclist_uerr(full_fig, leaders_list, utrue, green4)
% debug_MH_plots.plot_nClass_uerr(full_fig, mhgen,utrue,c_last,green4)
