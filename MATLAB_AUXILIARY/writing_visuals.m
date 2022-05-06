clear
close all
run AYfigprops.m

write_figs = false;
write_all_figs = true;
make_movie = true;
make_plots = true;
figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';

fig_pos = AYfig.fig_pos_gen(2, 3);
pos_full = [1 1 1728 1000];
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];
swirl_path_specs = AYfig.specs_gen('swirl_path', [fig_pos(5, 1:2), 500,500]);
movie_specs = AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500]);

movie1 = AYfig(movie_specs);

proc_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'pts';

par_id = 0;

sw1 = read_swirl(1, par_id);
sw30 = read_swirl(30, par_id);

sw_watch = sw1;

[sw1.contact_f, sw1.bead_contact_f, sw1.wall_contact_f] = sw1.find_contact_frames;

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  get data  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

if (make_plots)
    swirl_path_fig = AYfig(swirl_path_specs);

    writing_plots.plot_1bead_swirl(sw1,swirl_path_fig, green4, red5);
elseif (make_movie)
    movie1 = AYfig(movie_specs);

    movie_data = cell([1,2]);
    movie_data{:,1} = sw_watch.pos(:,1:2,:);
    movie_data{:,2} = ones(sw_watch.beads,3).*green4;

    % movie_specs = struct('Frame_vec', contact_f, 'dish', sw.dish);
    movie_specs = struct('Frame_vec', 1:sw_watch.Frames, 'dish', sw_watch.dish);

    swirl_group.make_movie_comp(movie1, movie_data, movie_specs, 'watch');
    pause
    movie1.play_movie(10,30);
    % movie1.frame_by_frame(1:length(movie_specs.Frame_vec), 'wait');
end

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  end plots  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

if (write_figs)
  % if (write_all_figs)
  %   figs_to_write = 1:length(figs);
  % end
  % AYfig.save_figs(figs, figs_to_write, save_type, save_dir);
  AYfig.save_fig(diagnostics_fig.fig, save_type, save_dir);
end
