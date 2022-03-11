clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

nbeads = 3;
par_id = 0;
ran_id = 0;

sw = read_swirl(nbeads, par_id);
sb = sw.spawn_best();
race = sw.read_race_data();
stat = read_stat(nbeads, 'maxmin', ran_id);
stat_best = stat.sw(stat.I_best(1)+1);
stat_true = stat.sw(stat.I_truest(1)+1);

sw.make_movie_comp3(movie1, sb, stat_best, stat_true);
