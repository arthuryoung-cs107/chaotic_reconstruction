function movie1 = watch_swirl(nbeads, par_id)
  run AYfigprops.m
  fig_pos = AYfig.fig_pos_gen(2, 3);
  movie1 = AYfig(AYfig.specs_gen('playback', [fig_pos(5, 1:2), 500,500] ));

  sw = read_swirl(nbeads, par_id);
  sw.load_filin();
  sw.make_movie(movie1, 'watch');
end
