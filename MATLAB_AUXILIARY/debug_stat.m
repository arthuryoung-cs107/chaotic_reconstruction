clear
close all
run AYfigprops.m

pos_top_row = [1 551 1728 460];
pos_full = [1 1 1728 1000];
pos_bottom_row = [0 1 1728 460];

stat_fig = AYfig(AYfig.specs_gen('stat_test_plots', pos_top_row));
stat_fig.init_tiles([1, 4]);

nbeads = 3;
test = 'maxmin'; %% parameter perturbation
ran_id = 1;

tic
stat = read_stat(nbeads,test,ran_id);
toc
swtrue = stat.sw0;
[t_vec, pos_res] = deal(swtrue.t_vec, stat.pos_res);
[INTpos_res, Dpos_res] = deal(swirl_group.compute_INT(t_vec, pos_res),swirl_group.compute_D(t_vec, pos_res));
[contact_f, bead_contact_f, wall_contact_f] = swtrue.find_contact_frames;
frame_stats = swirl_group.compute_residual_time_stats(INTpos_res,pos_res,Dpos_res,contact_f);

error_plots.plot_INTres_vs_time(stat_fig.ax_tile(1), red5, stat);
error_plots.plot_res_vs_frames(stat_fig.ax_tile(2), blue5, stat);
error_plots.plot_Dres_vs_time(stat_fig.ax_tile(3), green4, stat);
error_plots.plot_residual_stats(stat_fig.ax_tile(1), stat_fig.ax_tile(2), stat_fig.ax_tile(3), t_vec, contact_f, frame_stats);

error_plots.plot_param_error(stat_fig.ax_tile(4), purple1, stat);

%% video stuff

pause

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

sw0 = swtrue;
stbest = stat.spawn_swirl_i(stat.I_best(1));
sttruest = stat.spawn_swirl_i(stat.I_truest(1));
stleader = stat.spawn_swirl_i(stat.I_leader(1));

movie_data = cell([4,2]);
[movie_data{:,1}] = deal(sw0.pos(:,1:2,:),stbest.pos(:,1:2,:),sttruest.pos(:,1:2,:),stleader.pos(:,1:2,:));
[movie_data{:,2}] = deal(ones(sw0.beads,3).*green4,ones(stbest.beads,3).*orange1,ones(sttruest.beads,3).*red5, ones(stleader.beads,3).*blue5);

% movie_specs = struct('Frame_vec', 1:(stat.frscores(stat.I_leader(1))), 'dish', sw0.dish);
movie_specs = struct('Frame_vec', bead_contact_f, 'dish', sw0.dish);

swirl_group.make_movie_comp(movie1, movie_data, movie_specs, 'watch');
pause
% movie1.play_movie(10,30);
movie1.frame_by_frame(1:length(movie_specs.Frame_vec), 'wait');

% * NOTE: in pursuit of a fully interpretable model, we need to make a few mathematical connections between what we are plotting between tests. One thing that is worth seeing is how the singular value profile varies with the behaviour of the covariance approximation. If we can make any connection between the values of these graphs, we should be able to gain access to a ton of tools in both linear algebra and statistics/probability/Markov chains.
