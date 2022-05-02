clear
close all
run AYfigprops.m

pos_top_row = [1 551 1728 460];
pos_full = [1 1 1728 1000];
pos_bottom_row = [0 1 1728 460];

res_fig = AYfig(AYfig.specs_gen('stat_test_plots', pos_top_row));
res_fig.init_tiles([1, 3]);

alpha_fig = AYfig(AYfig.specs_gen('alpha_fig', pos_bottom_row));
alpha_fig.init_tiles([1, 4]);

nbeads = 3;
test = 'maxmin'; %% parameter perturbation
ran_id = 1;

tic
stat = read_stat(nbeads,test,ran_id);
toc
swtrue = stat.sw0;
[t_vec, pos_res] = deal(swtrue.t_vec, stat.pos_res);
[INTpos_res, Dpos_res] = deal(swirl_group.compute_INT(t_vec, pos_res),swirl_group.compute_D(t_vec, pos_res));
alpha_INTpos_res = swirl_group.compute_alpha_logD(t_vec, INTpos_res);
[contact_f, bead_contact_f, wall_contact_f] = swtrue.find_contact_frames;
frame_stats = swirl_group.compute_residual_time_stats(INTpos_res,pos_res,Dpos_res,contact_f);

error_plots.plot_INTres_vs_time(res_fig.ax_tile(1), red5, stat);
error_plots.plot_res_vs_frames(res_fig.ax_tile(2), blue5, stat);
error_plots.plot_Dres_vs_time(res_fig.ax_tile(3), green4, stat);
error_plots.plot_residual_stats(res_fig.ax_tile(1), res_fig.ax_tile(2), res_fig.ax_tile(3), t_vec, contact_f, frame_stats);

error_plots.plot_alphaINTres_vs_time(alpha_fig.ax_tile(1), orange1, stat);
error_plots.plot_alphaINTres_vs_time_bead(alpha_fig.ax_tile(2:end), orange1, stat);

% pause

% flim =[1, contact_f(1)+1];
% limind=flim(1):(contact_f(1)-3);
flim =[1, 400];
limind=flim(1):flim(2);

tlim = reshape(t_vec(flim), 1,[]);

alpha_lims = error_plots.lims_minmax(alpha_INTpos_res(:,limind));
ylim_mat = [error_plots.lims_minmax(INTpos_res(:,limind)); error_plots.lims_minmax(Dpos_res(:,limind)); alpha_lims; alpha_lims; alpha_lims; alpha_lims];

error_plots.adjust_lims([res_fig.ax_tile(1) res_fig.ax_tile(3) alpha_fig.ax_tile(1) alpha_fig.ax_tile(2) alpha_fig.ax_tile(3) alpha_fig.ax_tile(4)], ones(6,2).*tlim, ylim_mat);
error_plots.adjust_lims(res_fig.ax_tile(2), flim, error_plots.lims_minmax(pos_res(:,limind)));

% error_plots.plot_Daccres_vs_accres(fig_all.ax_tile(4), orange1, stat);
% error_plots.plot_INTres_vs_res(fig_all.ax_tile(5), purple1, stat);
% error_plots.plot_accres_vs_res(fig_all.ax_tile(6), pink1, stat);
%
% error_plots.plot_INTres_vs_time(fig_all.ax_tile(4), red5, stat);
% error_plots.plot_Dres_vs_time(fig_all.ax_tile(6), green4, stat);
