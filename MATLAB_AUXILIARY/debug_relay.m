clear
close all
run AYfigprops.m

write_figs = true;
write_all_figs = true;
figs_to_write = 0;
save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
save_type = 'pdf';

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  begin plots  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

fig_pos = AYfig.fig_pos_gen(2, 3);
pos_full = [1 1 1728 1000];
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];

relay_diagnostics_fig = AYfig(AYfig.specs_gen('relay_diagnostics',pos_full));
relay_diagnostics_fig.init_tiles([2, 3]);

nbeads = 3;
par_id = 0;
relay_id = 0;

relay = read_relay(nbeads, par_id, relay_id);

relay_plots.plot_3bead_event_stats(relay_diagnostics_fig.ax_tile, blue5, red5, relay);

%%%%%%%% ----------------------------------------------------------------------------------
%%%%%%%% ------------------------------  end plots  ---------------------------------
%%%%%%%% ---------------------------------------------------------------------------

if (write_figs)
  % if (write_all_figs)
  %   figs_to_write = 1:length(figs);
  % end
  % AYfig.save_figs(figs, figs_to_write, save_type, save_dir);
  AYfig.save_fig(relay_diagnostics_fig.fig, save_type, save_dir); 
end
