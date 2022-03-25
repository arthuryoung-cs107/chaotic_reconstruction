clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);
pos_bottom_row = [1 1 1440 345];

diagnostics_fig = AYfig(AYfig.specs_gen('generation diagnostics',pos_bottom_row));
diagnostics_fig.init_tiles([1, 3]);

nbeads = 3;
par_id = 0;
ran_id = 0;

sw = read_swirl(nbeads, par_id);
sb = sw.spawn_best();
race = sw.read_race_data();
race.F = 1200;
race.lambda = (1.0)*log(((1e16)-1)/(100-1));

for i=1:race.gen_count
  race.plot_diagnosticsi(diagnostics_fig, i);
  pause
end
