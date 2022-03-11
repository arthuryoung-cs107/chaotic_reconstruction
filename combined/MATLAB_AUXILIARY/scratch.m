clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);

movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

movie1.init_movie(1);
figdims = movie1.get_dims()
lims = [-1, 1, -1, 1];

hold(movie1.ax, 'on');
plot(movie1.ax, 0, 0, '.k', 'MarkerSize', figdims(4)*1);
plot(movie1.ax, 0, 0, 'ok', 'MarkerSize', figdims(4)*1);

figdims = movie1.get_dims()
