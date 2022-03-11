clear
close all
run AYfigprops.m

fig_pos = AYfig.fig_pos_gen(2, 3);
movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

proc_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'pts';

nbeads = 20;
par_id = 0;

exp_name = ['swirl' num2str(nbeads) '.odr/'];
file_name = [dat_name '.' num2str(par_id)];
sw20 = ODR_data(proc_name, exp_name, file_name);
sw20.make_movie(movie1);
