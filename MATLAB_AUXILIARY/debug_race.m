clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);
pos_bottom_row = [1 1 1440 345];

diagnostics_fig = AYfig(AYfig.specs_gen('generation diagnostics',pos_bottom_row));
diagnostics_fig.init_tiles([1, 3]);

nbeads = 3;
par_id = 0;
stat_test_type = 'gauss';
ran_id = 0;

sw3 = read_swirl(nbeads, par_id);
stat3 = read_stat(nbeads, stat_test_type, ran_id);
race3 = read_race(nbeads, par_id);

stat3.
swbst3 = race3.spawn_swbest();

% race.F = 1200;
% race.lambda = (1.0)*log(((1e16)-1)/(100-1));

% for i=1:race.gen_count
%   race.plot_diagnosticsi(diagnostics_fig, i);
%   pause
% end
