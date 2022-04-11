clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];

stat_diagnostics_fig = AYfig(AYfig.specs_gen('stat_diagnostics',pos_top_row));
stat_diagnostics_fig.init_tiles([1, 3]);

race_diagnostics_fig = AYfig(AYfig.specs_gen('race_diagnostics',pos_bottom_row));
race_diagnostics_fig.init_tiles([1, 3]);

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

nbeads = 3;
par_id = 0;
stat_test_type = 'maxmin';
ran_id = 0;

sw = read_swirl(nbeads, par_id);
stat = read_stat(nbeads, stat_test_type, ran_id);
race = read_race(nbeads, par_id);

stgp = stat.spawn_swirlgroup();
stat = stat_data.init_statistics(stat, stgp);
stat = stat_data.init_stat_generation(stat, stgp);

stat.plot_gen_scores(stat_diagnostics_fig.ax_tile(1));

race.F = 1200;
race.lambda = (1.0)*log(((1e16)-1)/(100-1));

race.plot_diagnosticsi(race_diagnostics_fig, 1);

sw0 = stgp.sw0;
swbst = race.spawn_swbest();
stbest = stgp.spawn_swirl_i(stat.I_best(1));
sttruest = stgp.spawn_swirl_i(stat.I_truest(1));

movie_data = cell([4,2]);
[movie_data{:,1}] = deal(sw0.pos(:,1:2,:),swbst.pos(:,1:2,:), stbest.pos(:,1:2,:),sttruest.pos(:,1:2,:));
[movie_data{:,2}] = deal(ones(sw0.beads,3).*green4,ones(swbst.beads,3).*blue5, ones(stbest.beads,3).*orange1,ones(sttruest.beads,3).*red5);

movie_specs = struct('Frames', stgp.sw0.Frames, 'dish', sw0.dish);

swirl_group.make_movie_comp(movie1, movie_data, movie_specs);
