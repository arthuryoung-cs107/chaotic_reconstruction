clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);
movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

nbeads = 3;
par_id = 0;
ran_id = 0;

sw30_0 = read_swirl(30, 0);
sw30_3 = read_swirl(30, 3);

sw30_0.make_movie_comp(movie1, sw30_3)
