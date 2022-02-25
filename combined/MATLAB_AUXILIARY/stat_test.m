clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 2);

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
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold on

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(4, 1:2), 350, 400] ));

dat_dir_name = '../dat_dir/';
stat_name = 'stat3.odr/';
dat_name = 'rand';

stat3 = stat_data(dat_dir_name, stat_name, dat_name, 50);
stat3.plot_frame_error(figs(1), blue5);
stat3.plot_param_error(figs(2), green4);
stat3.plot_param_frame_error(figs(3), red5);
stat3.make_moviei(movie1, 1);
