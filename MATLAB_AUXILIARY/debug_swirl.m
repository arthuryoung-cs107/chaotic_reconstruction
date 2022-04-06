clear
close all
run AYfigprops.m

fig_pos = AYfig.fig_pos_gen(2, 3);
movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

proc_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'pts';

nbeads = 3;
par_id = 0;

sw = read_swirl(nbeads, par_id);
