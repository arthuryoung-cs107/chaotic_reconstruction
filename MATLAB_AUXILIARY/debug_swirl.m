clear
close all
run AYfigprops.m

fig_pos = AYfig.fig_pos_gen(2, 3);
movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500]));

proc_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'pts';

nbeads = 3;
par_id = 0;

sw = read_swirl(nbeads, par_id);

[contact_f, bead_contact_f, wall_contact_f] = sw.find_contact_frames;

movie_data = cell([1,2]);
movie_data{:,1} = sw.pos(:,1:2,:);
movie_data{:,2} = ones(sw.beads,3).*green4;

movie_specs = struct('Frame_vec', contact_f, 'dish', sw.dish);

swirl_group.make_movie_comp(movie1, movie_data, movie_specs, 'watch');
pause
% movie1.play_movie(10,30);
movie1.frame_by_frame(1:length(movie_specs.Frame_vec), 'wait');
