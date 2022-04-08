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

nbeads = 3;
par_id = 0;
stat_test_type = 'maxmin';
ran_id = 1;

sw = read_swirl(nbeads, par_id);
stat = read_stat(nbeads, stat_test_type, ran_id);
race = read_race(nbeads, par_id);

stgp = stat.spawn_swirlgroup();
stat = stat_data.init_statistics(stat, stgp);
stat = stat_data.init_stat_generation(stat, stgp);

stat.plot_gen_scores(stat_diagnostics_fig.ax_tile(1));

swbst3 = race.spawn_swbest();

race.F = 1200;
race.lambda = (1.0)*log(((1e16)-1)/(100-1));

race.plot_diagnosticsi(race_diagnostics_fig, 1);
