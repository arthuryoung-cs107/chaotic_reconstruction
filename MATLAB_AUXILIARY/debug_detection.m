clear
close all
run AYfigprops.m

write_figs = false;
write_all_figs = true;
figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';

make_movie = false;

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  get data  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

nbeads = 3;
par_id = 0;
relay_id = 1;

relay = read_relay(nbeads, par_id, relay_id);
stat = read_stat(nbeads,'maxmin',1);
swtrue = stat.sw0;

det = relay.read_relay_test(swtrue,0,0);

Frame_end = det.Frame_end;
Frame_vec = 1:Frame_end;
Frame_limvec = 1:170;
swtrue_pos = swtrue.pos(:,1:2,Frame_vec);
stat_parmat = stat.params_mat(3:end, :);
swtrue_params = swtrue.params(3:end);

det_bcell = det.frame_posmat2beadcell(det.sim_pos_mat);
sim_beadres = det.frame_simbeadres(swtrue_pos,det_bcell);

beadres = det.frame_resmat2beadcell(det.pos_res_mat);
beadalpha = det.frame_resmat2beadcell(det.alpha_mat);
beadINTres = det.frame_resmat2beadcell(det.INTmat);
beadres_matcomp = det.beadres_matcomp();


%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  begin plots  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

fig_pos = AYfig.fig_pos_gen(2, 3);
pos_full = [1 1 1728 1000];
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];

test_fig = AYfig(AYfig.specs_gen('test_diagnostics',pos_full));
test_fig.init_tiles([3, 3]);


relay_plots.plot_cell_vs_frames(test_fig.ax_tile(1:3), blue5, Frame_vec, beadres, 'Frames', 'position residual', 'position_residual_vs_Frames_bead')
relay_plots.plot_cell_vs_frames(test_fig.ax_tile(4:6), red5, Frame_vec, beadres_matcomp, 'Frames', 'matcomp position residual', 'posres_matcomp_vs_Frames_bead')
relay_plots.plot_cell_vs_frames(test_fig.ax_tile(7:9), orange1, Frame_vec, beadalpha, 'Frames', 'alpha', 'alpha_vs_Frames')


% relay_plots.plot

% relay_plots.plot_param_error(test_fig.ax_tile(7), green4,stat_parmat,swtrue_params)
% relay_plots.plot_param_error(test_fig.ax_tile(8), green4,det.params,swtrue_params)

if (make_movie)
    det_pcell = det.frame_posmat2poolcell(det.sim_pos_mat);
    movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));


    [refsimpos, simpos] = det.datpos2simpos(swtrue.dish,det.frame_posmat2poolcell(det.pos_mat), 1);


    movie_data = cell([3,2]);
    [movie_data{:,1}] = deal(swtrue_pos, simpos, refsimpos);
    [movie_data{:,2}] = deal(ones(swtrue.beads,3).*green4,ones(swtrue.beads,3).*orange1, ones(swtrue.beads,3).*purple1);

    movie_specs = struct('Frame_vec', 1:det.Frame_end, 'dish', swtrue.dish);

    swirl_group.make_movie_comp(movie1, movie_data, movie_specs, 'watch');
    pause
    movie1.play_movie(10,30);

end

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  end plots  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------


if (write_figs)
  % if (write_all_figs)
  %   figs_to_write = 1:length(figs);
  % end
  % AYfig.save_figs(figs, figs_to_write, save_type, save_dir);
  AYfig.save_fig(test_fig.fig, save_type, save_dir);
end
