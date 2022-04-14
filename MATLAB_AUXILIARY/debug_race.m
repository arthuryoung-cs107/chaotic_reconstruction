clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];

race_diagnostics_fig = AYfig(AYfig.specs_gen('race_diagnostics',pos_top_row));
race_diagnostics_fig.init_tiles([1, 3]);

% movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

nbeads = 3;
par_id = 0;
stat_test_type = 'maxmin';
ran_id = 0;

sw = read_swirl(nbeads, par_id);
race = read_race(nbeads, par_id);

race.plot_diagnosticsi(race_diagnostics_fig, 1);

swbst = race.spawn_swbest();
