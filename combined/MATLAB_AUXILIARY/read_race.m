clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));
movie2 = AYfig(AYfig.specs_gen('playback', [fig_pos(1, 1:2), 500,500] ));

nbeads = 3;
par_id = 0;
proc_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'pts';

exp_name = ['swirl' num2str(nbeads) '.odr/'];
file_name = [dat_name '.' num2str(par_id)];

sw = ODR_data(proc_name, exp_name, file_name);
sw.load_filin();
swb = sw.spawn_best();

sw.make_movie_comp(movie1, swb);
% movie1.play_movie(10, 30);
