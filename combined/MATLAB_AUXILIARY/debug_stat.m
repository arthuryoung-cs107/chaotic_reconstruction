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

figs(3) = AYfig.figure(AYfig.specs_gen('position error vs parameter index', fig_pos(3, : )));
xlabel('parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('position error', 'Interpreter', 'Latex', 'Fontsize', 14)
% set(gca, 'YScale', 'log')
hold on

figs(4) = AYfig.figure(AYfig.specs_gen('position error vs parameter error', fig_pos(4, : )));
xlabel('param error', 'Interpreter', 'Latex', 'Fontsize', 14)
ylabel('position error', 'Interpreter', 'Latex', 'Fontsize', 14)
set(gca, 'YScale', 'log')
hold on

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

nbeads = 3;
par_id = 0;
test = 'maxmin';
ran_id = 0;
dat_dir_name = '../dat_dir/';
pov_dir = '../POV_AUXILIARY/';
dat_name = 'rand';

exp_name = ['stat', num2str(nbeads), '_' test num2str(ran_id) '.odr/'];

stat = stat_data(dat_dir_name, exp_name, dat_name, 0);
stat.plot_frame_error(figs(1), blue5);
stat.plot_param_error(figs(2), green4);
stat.plot_param_index_error(figs(3), orange1)
stat.plot_param_pos_error(figs(4), red5);

stat.make_movieijk(movie1, 1, stat.I_best(1)+1, stat.I_truest(1)+1);
% movie1.play_movie(10, 45)

% stat3.make_POVray_inputsi(1, pov_dir);
