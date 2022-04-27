clear
close all
run AYfigprops.m

fig_all = AYfig(AYfig.specs_gen('stat_test_plots', [3 34 1726 967]));
fig_all.init_tiles([2, 3]);

nbeads = 3;
test = 'maxmin'; %% parameter perturbation
ran_id = 1;
nlead = 500;


tic
stat = read_stat(nbeads,test,ran_id);
toc

error_plots.plot_param_error(fig_all.ax_tile(1), green4, stat);
