clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);

figs(1) = AYfig(AYfig.specs_gen('frscores', fig_pos(6, :) ));
figs(2) = AYfig(AYfig.specs_gen('sigmas', fig_pos(4, :) ));

nbeads = 3;
par_id = 0;
ran_id = 0;

sw = read_swirl(nbeads, par_id);
sb = sw.spawn_best();
race = sw.read_race_data();


for i=1:race.gen_count
  race.frscore_histogrami(figs(1), i)
  race.sigmas_ploti(figs(2), i)
  pause
end
