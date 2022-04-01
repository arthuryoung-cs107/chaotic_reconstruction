clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);
pos_bottom_row = [1 1 1440 345];

diagnostics_fig = AYfig(AYfig.specs_gen('generation diagnostics',pos_bottom_row));
diagnostics_fig.init_tiles([1, 3]);

nbeads = 3;
par_id = 0;
stat_test_type = 'maxmin'
ran_id = 0;

swirl3 = read_swirl(nbeads, par_id);
race3 = swirl3.read_race();
stat3 = swirl3.read_stat(par_id, stat_test_type, ran_id);

race3 = read_race(nbeads, par_id);

odr3 = read_ODR(nbeads, par_id);
stat3 = odr.spawn_stat_data();
odr
.spawn_best();
race = ODR.read_race();

sw = ODR.spawn_read_race();
sw
race.F = 1200;
race.lambda = (1.0)*log(((1e16)-1)/(100-1));


stat_data.get_swirl_group;

for i=1:race.gen_count
  race.plot_diagnosticsi(diagnostics_fig, i);
  pause
end
