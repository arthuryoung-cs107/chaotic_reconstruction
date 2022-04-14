clear
close all
run AYfigprops.m

fig_all = AYfig(AYfig.specs_gen('stat_test_plots', [3 34 1726 967]));
fig_all.init_tiles([2, 2]);

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

nbeads = 3;
test = 'maxmin'; %% parameter perturbation
ran_id = 1;

%% this effort is to produce a learned quantity that is FULLY PHYSICALLY INTERPRETABLE

stat = read_stat(nbeads,test,ran_id);

stat.plot_frame_error(fig_all.ax_tile(1), blue5, red5, stat.sw0.t_vec, stat.pos_err, stat.I_best(1), stat.I_truest(1));
stat.plot_param_error(fig_all.ax_tile(2), green4, stat.par_err, stat.sw0.params, stat.I_best(1), stat.I_truest(1));
stat.plot_err_vs_accerr(fig_all.ax_tile(3), orange1, stat.pos_err, stat.I_best(1), stat.I_truest(1));
stat.plot_derr_vs_err(fig_all.ax_tile(4), orange1, stat.sw0.dish(2,1), stat.pos_err, stat.I_best(1), stat.I_truest(1));

% sw0 = stgp.sw0;
% stbest = stgp.spawn_swirl_i(stat.I_best(1));
% sttruest = stgp.spawn_swirl_i(stat.I_truest(1));
%
% movie_data = cell([3,2]);
% [movie_data{:,1}] = deal(sw0.pos(:,1:2,:),stbest.pos(:,1:2,:),sttruest.pos(:,1:2,:));
% [movie_data{:,2}] = deal(ones(sw0.beads,3).*green4,ones(stbest.beads,3).*blue5,ones(sttruest.beads,3).*red5);
%
% movie_specs = struct('Frames', stgp.sw0.Frames, 'dish', sw0.dish);
%
% swirl_group.make_movie_comp(movie1, movie_data, movie_specs);

%* NOTE: in pursuit of a fully interpretable model, we need to make a few mathematical connections between what we are plotting between tests. One thing that is worth seeing is how the singular value profile varies with the behaviour of the covariance approximation. If we can make any connection between the values of these graphs, we should be able to gain access to a ton of tools in both linear algebra and statistics/probability/Markov chains.
