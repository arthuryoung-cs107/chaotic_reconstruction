clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);
pos_full = [1 1 1728 1000];

relay_diagnostics_fig = AYfig(AYfig.specs_gen('race_diagnostics',pos_full));
relay_diagnostics_fig.init_tiles([2, 3]);

nbeads = 3;
par_id = 0;
relay_id = 0;

relay = read_relay(nbeads, par_id, relay_id);
