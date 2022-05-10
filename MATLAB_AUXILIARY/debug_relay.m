clear
close all
run AYfigprops.m

write_figs = true;
write_all_figs = true;
figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';

fig_pos = AYfig.fig_pos_gen(2, 3);
pos_full = [1 1 1728 1000];
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];

% event_fig = AYfig(AYfig.specs_gen('event_data',pos_full));
% event_fig.init_tiles([3, 3]);

diagnostics_fig = AYfig(AYfig.specs_gen('relay_diagnostics',pos_full));
diagnostics_fig.init_tiles([4, 3]);

convergence_fig = AYfig(AYfig.specs_gen('convergence',fig_pos(1,:)));

figs = [diagnostics_fig, convergence_fig];
%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  get data  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

nbeads = 20;
par_id = 0;
relay_id = 3;

relay = read_relay(nbeads, par_id, relay_id);
stat = read_stat(3,'maxmin',1);
swtrue = stat.sw0;
[contact_f, bead_contact_f, wall_contact_f] = swtrue.find_contact_frames;

% gen_last = relay.specs.gen_max;
gen_last = 154;
true_params = swtrue.params(3:end);
last_indices = (gen_last-length(diagnostics_fig.ax_tile)+1):gen_last;
first_indices = 1:length(diagnostics_fig.ax_tile);
uniform_indices = round(linspace(1, double(gen_last),length(diagnostics_fig.ax_tile)));
[uniform_indices(1), uniform_indices(length(diagnostics_fig.ax_tile))] = deal(1, gen_last);

gen_ind = uniform_indices;

[relay.gen_ind, relay.gen] = relay.read_generations(1:gen_last, gen_last);

% relay.write_relay_test(stat.params_mat(3:end, :),0,1);

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  begin plots  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

% relay_plots.plot_3bead_event_stats(event_fig.ax_tile, [blue5; red5; green4], relay);

relay_plots.plot_convergence(convergence_fig.ax, green4, relay, 1:gen_last, true_params)
% writing_plots.plot_param_error(convergence_fig.ax, green4, stat.params_mat(3:end,:), true_params)

% relay_plots.plot_gen_param_error(diagnostics_fig.ax_tile, green4, relay, true_params, gen_ind);
% relay_plots.plot_gen_netresiduals(diagnostics_fig.ax_tile, red5, relay, gen_ind);
% relay_plots.plot_gen_probabilities(diagnostics_fig.ax_tile, purple1, relay, gen_ind);
% relay_plots.plot_gen_zeta(diagnostics_fig.ax_tile, purple1, relay, gen_ind);
% relay_plots.plot_gen_weights(diagnostics_fig.ax_tile, pink1, relay, gen_ind);
% relay_plots.plot_gen_duplication_count(diagnostics_fig.ax_tile, purple2, relay, gen_ind);
% relay_plots.plot_gen_weight_vs_dup(diagnostics_fig.ax_tile, blue1, relay, gen_ind);

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  end plots  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

if (write_figs)
  if (write_all_figs)
    figs_to_write = 1:length(figs);
  end
  for i=1:length(figs)
      AYfig.save_fig(figs(i).fig,save_type, save_dir);
  end
  % AYfig.save_figs(figs., figs_to_write, save_type, save_dir);
  % AYfig.save_fig(diagnostics_fig.fig, save_type, save_dir);
end
