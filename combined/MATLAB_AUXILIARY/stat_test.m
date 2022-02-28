clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);

figs(1) = AYfig.figure(AYfig.specs_gen('position error vs frame', fig_pos(1, : )));
xlabel('Frame', 'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('position error', 'Interpreter', 'Latex', 'Fontsize', 14)
set(gca, 'YScale', 'log')
hold on

figs(2) = AYfig.figure(AYfig.specs_gen('parameter error vs index', fig_pos(2, : )));
xlabel('parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('param error', 'Interpreter', 'Latex', 'Fontsize', 14)
hold on

figs(3) = AYfig.figure(AYfig.specs_gen('position error vs parameter error', fig_pos(3, : )));
xlabel('param error', 'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('position error', 'Interpreter', 'Latex', 'Fontsize', 14)
% set(gca, 'YScale', 'log')
hold on

figs(4) = AYfig.figure(AYfig.specs_gen('position error vs parameter error', fig_pos(4, : )));
xlabel('param error', 'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('position error', 'Interpreter', 'Latex', 'Fontsize', 14)
set(gca, 'YScale', 'log')
hold on

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

stat_name_maxmin = 'stat3_maxmin.odr/';
stat_name_gauss = 'stat3_gauss.odr/';

stat_name = stat_name_gauss;
dat_dir_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'rand';

stat3 = stat_data(dat_dir_name, stat_name, dat_name, 100);
stat3.plot_frame_error(figs(1), blue5);
stat3.plot_param_error(figs(2), green4);
stat3.plot_param_index_error(figs(3), orange1)
stat3.plot_param_pos_error(figs(4), red5);

stat3.make_POVray_inputsi(1, pov_dir);

% stat3.make_moviei(movie1, 1);

% movie1.play_movie(10, 60)
