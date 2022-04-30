clear
close all
run AYfigprops.m

write_figs = false;
write_all_figs = true;
figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  get data  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

nbeads = 3;
par_id = 0;
relay_id = 1;

relay = read_relay(nbeads, par_id, relay_id);
stat = read_stat(nbeads,'maxmin',0);
swtrue = stat.sw0;

% relay.test = relay.read_relay_test(0,0);

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  begin plots  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

fig_pos = AYfig.fig_pos_gen(2, 3);
pos_full = [1 1 1728 1000];
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];

% test_fig = AYfig(AYfig.specs_gen('test_diagnostics',pos_full));
% test_fig.init_tiles([2, 3]);
%
% relay_plots.plot_posres_vs_time(test_fig.ax_tile(2), blue5, relay)
%
% pause

diagnostics_fig = AYfig(AYfig.specs_gen('relay_diagnostics',pos_full));
diagnostics_fig.init_tiles([3, 3]);

gen_last = relay.specs.gen_last;
true_params = swtrue.params(3:end);
last_indices = (gen_last-length(diagnostics_fig.ax_tile)+1):gen_last;
first_indices = 1:length(diagnostics_fig.ax_tile);
plot_indices =first_indices;

% relay_plots.plot_3bead_event_stats(diagnostics_fig.ax_tile, blue5, red5, green4, relay);
relay_plots.plot_gen_param_error(diagnostics_fig.ax_tile, green4, relay, true_params, plot_indices);
% relay_plots.plot_gen_weights(diagnostics_fig.ax_tile, pink1, relay, plot_indices);
% relay_plots.plot_gen_duplication_count(diagnostics_fig.ax_tile, purple1, relay, plot_indices);
% relay_plots.plot_gen_weight_vs_dup(diagnostics_fig.ax_tile, blue1, relay, plot_indices);



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
