clear
close all
run AYfigprops.m

fig_all = AYfig(AYfig.specs_gen('stat_test_plots', [3 34 1726 967]));
fig_all.init_tiles([2, 2]);

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

nbeads = 3;
test = 'maxmin'; %% parameter perturbation
ran_id = 0;

%% this effort is to produce a learned quantity that is FULLY PHYSICALLY INTERPRETABLE

stat = read_stat(nbeads,test,ran_id);
stgp = stat.spawn_swirlgroup();

[stat.par_err, stat.pos_err, stat.pos_err_acc] = stgp.compute_error();
[stat.I_best, stat.I_truest, stat.par_cov] = stgp.compute_statistics(stat.par_err, stat.pos_err, stat.pos_err_acc);

stgp.plot_frame_error(fig_all.ax_tile(1), blue5, stat.pos_err, stat.I_best(1), stat.I_truest(1));
stgp.plot_param_error(fig_all.ax_tile(2), green4, stat.par_err, stat.I_best(1), stat.I_truest(1));
stgp.plot_param_covariance(fig_all.ax_tile(3), orange1, stat.par_cov, )
% stat.plot_param_pos_error(figs(4), red5);

%* NOTE: in pursuit of a fully interpretable model, we need to make a few mathematical connections between what we are plotting between tests. One thing that is worth seeing is how the singular value profile varies with the behaviour of the covariance approximation. If we can make any connection between the values of these graphs, we should be able to gain access to a ton of tools in both linear algebra and statistics/probability/Markov chains.

% stat.make_movieijk(movie1, 1, stat.I_best(1)+1, stat.I_truest(1)+1);
