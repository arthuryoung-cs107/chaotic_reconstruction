clear
close all
run AYfigprops.m
fig_pos = AYfig.fig_pos_gen(2, 3);
pos_full = [1 1 1728 1000];
pos_top_row = [1 551 1728 460];
pos_bottom_row = [0 1 1728 460];

lineage_fig = AYfig(AYfig.specs_gen('race_diagnostics',pos_full));
race_diagnostics_fig = AYfig(AYfig.specs_gen('race_diagnostics',pos_bottom_row));
lineage_fig.init_tiles([2, 3]);
race_diagnostics_fig.init_tiles([1, 3]);


nbeads = 3;
par_id = 0;
stat_test_type = 'maxmin';
ran_id = 1;
race_id = 1;

stat = read_stat(nbeads, stat_test_type, ran_id);
race = read_race(nbeads, par_id, race_id);

race.plot_diagnosticsi(race_diagnostics_fig, race.gen_count);

race.plot_final_lineage(lineage_fig, stat);

% for i = 1:race.gen_count
%     race.plot_diagnosticsi(race_diagnostics_fig, i);
%     pause
% end
%
% swbst = race.spawn_swbest();
% movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));
% movie_data = cell([2,2]);
% [movie_data{:,1}] = deal(sw.pos(:,1:2,:),swbst.pos(:,1:2,:));
% [movie_data{:,2}] = deal(ones(sw.beads,3).*green4, ones(swbst.beads,3).*blue5);

% movie_specs = struct('Frames', sw.Frames, 'dish', sw.dish);
% swirl_group.make_movie_comp(movie1, movie_data, movie_specs, 'watch');
% pause
