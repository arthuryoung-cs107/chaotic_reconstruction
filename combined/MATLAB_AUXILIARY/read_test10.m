clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));
movie2 = AYfig(AYfig.specs_gen('playback', [fig_pos(1, 1:2), 500,500] ));

race_name = 'race_3beads.odr/';
race10_name = 'race_10beads.odr/';

exp_name = race10_name;
dat_dir_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'pts';

odr = ODR_data(dat_dir_name, exp_name, dat_name);
% odr.load_filin();

odr.make_movie(movie1);
movie1.play_movie(10, 60);


% swbest = ODR_data(dat_dir_name, racebest_name, dat_name);
% odr.make_movie_comp(movie1, swbest);
