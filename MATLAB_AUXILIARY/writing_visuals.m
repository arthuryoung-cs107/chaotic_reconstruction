clear
close all
run AYfigprops.m

write_figs = false;
write_all_figs = true;

make_beadalpha_plots=false;
make_alpha_plots=false;
make_collision_plots=true;
make_stat_plots=false;
make_convergence_plots=false;
make_path_plots = false;
make_movie = false;

figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';

fig_pos = AYfig.fig_pos_gen(2, 3);
page_pos1x2 = [1 701 920 300];
page_pos2x2 = [1 1 920 600];
page_pos3x2 = [1 1 920 1000];
halfpage = [1 1 730 600];
pos_full = [1 1 1728 1000];
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];

beadalpha_specs = AYfig.specs_gen('beadalpha_plots', page_pos1x2);
alpha_specs = AYfig.specs_gen('alpha_plots', page_pos2x2);
collision_specs = AYfig.specs_gen('collision_plots', [1, 1, 700, 700]);
stat_specs = AYfig.specs_gen('baseline_res_plots', page_pos1x2);
convergence_specs = AYfig.specs_gen('convergence_plots', pos_top_row);
swirl_path_specs = AYfig.specs_gen('swirl_path', [fig_pos(5, 1:2), 500,500]);
movie_specs = AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500]);

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

if (make_beadalpha_plots)
    beadalpha_fig = AYfig(beadalpha_specs, false);
    beadalpha_fig.init_tiles([1, 3]);

    stat3 = read_stat(3,'maxmin',0);
    test3 = test_group(3, 0, 2, 2);

    writing_plots.plot_alpha_bead(beadalpha_fig, purple1, purple2, purple3, test3, stat3.sw0.t_vec)

    % AYfig.save_fig(beadalpha_fig.fig, save_type, save_dir);

elseif (make_alpha_plots)
    alpha_fig = AYfig(alpha_specs, false);
    alpha_fig.init_tiles([2, 2]);

    stat3 = read_stat(3,'maxmin',1);
    test3 = test_group(3, 0, 2, 2);

    [t_vec, cl] = deal(stat3.sw0.t_vec, stat3.sw0.params(13));
    pos_res = cl*cl*stat3.pos_res;

    INTpos_res = swirl_group.compute_INT(t_vec, pos_res);
    testINTpos_res = swirl_group.compute_INT(t_vec, test3.pos_res);
    alpha_INTpos_res = swirl_group.compute_alpha_logD(t_vec, INTpos_res);
    testalpha_INTpos_res = swirl_group.compute_alpha_logD(t_vec, testINTpos_res);

    writing_plots.plot_alpha_truenoise(alpha_fig, orange1, purple1, red5, blue1, t_vec,alpha_INTpos_res,testalpha_INTpos_res, INTpos_res, testINTpos_res)

    % AYfig.save_fig(alpha_fig.fig, save_type, save_dir);

    % swirl_group.write_test(stat3.params_mat,3,0,2,2);
elseif (make_collision_plots)
    collision_fig = AYfig(collision_specs, false);
    writing_plots.set_collision_fig(collision_fig);

    stat3 = read_stat(3,'maxmin',1);
    [stat3.sw0.contact_f, stat3.sw0.bead_contact_f, stat3.sw0.wall_contact_f] = stat3.sw0.find_contact_frames;

    writing_plots.plot_collision_frames(collision_fig, stat3, blue5, red5, green4, orange1);
    AYfig.save_fig(collision_fig.fig, save_type, save_dir);
elseif (make_stat_plots)
    stat_fig = AYfig(stat_specs);
    stat_fig.init_tiles([1, 2]);
    stat3 = read_stat(3,'maxmin',1);
    writing_plots.plot_stat_residuals(stat_fig, blue5, red5, stat3);
    % AYfig.save_fig(stat_fig.fig, save_type, save_dir);
elseif (make_convergence_plots)
    convergence_fig = AYfig(convergence_specs);
    convergence_fig.init_tiles([1, 3]);

    stat = read_stat(3,'maxmin',0);
    swtrue = stat.sw0;
    true_params = swtrue.params(3:end);

    relay3 = read_relay(3, 0, 5);
    gen_last3=324;
    [relay3.gen_ind, relay3.gen] = relay3.read_generations(1:gen_last3, gen_last3);
    writing_plots.plot_convergence(convergence_fig.ax_tile(1), green4, relay3, 1:gen_last3, true_params)
    writing_plots.plot_K_vs_gamma(convergence_fig.ax_tile(2), green4, relay3, 1:gen_last3, true_params)

    relay20 = read_relay(20, 0, 2);
    gen_last20=350;
    [relay20.gen_ind, relay20.gen] = relay20.read_generations(1:gen_last20, gen_last20);
    writing_plots.plot_convergence(convergence_fig.ax_tile(3), blue3, relay20, 1:gen_last20, true_params)
    writing_plots.plot_K_vs_gamma(convergence_fig.ax_tile(2), blue3, relay20, 1:gen_last20, true_params)

elseif (make_path_plots)
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
  % AYfig.save_fig(diagnostics_fig.fig, save_type, save_dir);
end
