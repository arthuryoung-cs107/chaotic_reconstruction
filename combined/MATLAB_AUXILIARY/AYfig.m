classdef AYfig < handle
  properties
    fig;
    ax;
    props_in;


    %% movie stuff
    movie_gen;
  end
  methods
    function obj = AYfig(props_in_)
      obj.fig = figure;
      obj.ax = gca;
      obj.props_in = props_in_;
      for i=1:size(props_in_, 1)
        obj.fig.set(props_in_{i, 1}, props_in_{i, 2});
      end
    end
    function init_movie(obj, Frames_)
      str(Frames_) = struct('cdata', [], 'colormap', []);
      obj.movie_gen = str;
      obj.ax.NextPlot = 'replaceChildren';
      % obj.fig.Visible = 'off';
    end
    function play_movie(obj)
      obj.fig.Visible = 'on';
      movie(obj.ax, obj.movie_gen);
    end
  end
  methods(Static)
    function fig_out = figure(props_in_)
      fig_out = figure;
      for i=1:size(props_in_, 1)
        fig_out.set(props_in_{i, 1}, props_in_{i, 2});
      end
    end
    function struct_out = specs_gen(name_in_, pos_in_)
      struct_out = {'Name', name_in_; 'Renderer', 'painters'; 'Position', pos_in_;};
    end
    function fig_pos_out = fig_pos_gen(rows_, cols_)
      rows = rows_;
      if rows > 10;
        rows = 10;
      end
      cols = cols_;
      if cols > 10;
        cols = 10;
      end
      fig_pos_out = zeros(rows*cols, 4);
      H = 850;
      W = 1500;
      del = 10;
      hbar = 80;
      h = floor(H/rows-hbar);
      w = floor((W)/cols);
      for i=1:rows*cols
        fig_pos_out(i, 3) = w;
        fig_pos_out(i, 4) = h;
      end
      k = 1;
      p = rows*cols;
      for i=1:rows
        for j=1:cols
          fig_pos_out(k, 1) = w*(j-1);
          fig_pos_out(p, 2) = (h + hbar)*(i-1);
          k = k + 1;
          p = p - 1;
        end
      end
    end

  end
end
