clear
close all
run AYfigprops.m

fig_pos = AYfig.fig_pos_gen(2, 3);
movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500]));

proc_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'pts';

nbeads = 3;

sw0 = read_swirl(nbeads, 0);
sw5 = read_swirl(nbeads, 5);

movie_data = cell([2,2]);
[movie_data{:,1}] = deal(sw0.pos(:,1:2,:), sw5.pos(:,1:2,:));
[movie_data{:,2}] = deal(ones(nbeads,3).*green4, ones(nbeads,3).*blue5);

movie_specs = struct('Frame_vec', 1:sw0.Frames, 'dish', sw0.dish);

swirl_group.make_movie_comp(movie1, movie_data, movie_specs, 'watch');
pause
movie1.play_movie(10,30);
