clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));
movie2 = AYfig(AYfig.specs_gen('playback', [fig_pos(1, 1:2), 500,500] ));

exp_name1 = 'race_30beads.odr/';
exp_name2 = 'race2_30beads.odr/';

dat_dir_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'pts';

odr1 = ODR_data(dat_dir_name, exp_name1, dat_name);
odr2 = ODR_data(dat_dir_name, exp_name2, dat_name);
% odr.load_filin();

odr1.make_movie_comp(movie1, odr2);
movie1.play_movie(10, 60);


% swbest = ODR_data(dat_dir_name, racebest_name, dat_name);
% odr.make_movie_comp(movie1, swbest);
