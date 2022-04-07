classdef AYfig < handle
  properties
    fig;
    ax;
    props_in;


    %% movie stuff
    movie_gen;

    %% tiled layout stuff
    tile;
    tile_dims;
    tile_num;
    ax_tile;
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
      axdims = get(obj.ax, 'Position');
      axdims(3:4) = min(axdims(3:4));
      set(obj.ax, 'Position', axdims);
      % obj.fig.Visible = 'off';
    end
    function init_tiles(obj, tile_dims_)
      clf(obj.fig);
      obj.tile_dims = tile_dims_;
      obj.tile_num = obj.tile_dims(1)*obj.tile_dims(2);
      obj.tile = tiledlayout(obj.fig, obj.tile_dims(1), obj.tile_dims(2));
      obj.tile.TileSpacing = 'compact';
      obj.tile.Padding = 'compact';

      obj.ax_tile = gobjects(obj.tile_num, 1);
      for i=1:obj.tile_num
        obj.ax_tile(i) = nexttile(obj.tile);
      end
    end
    function play_movie(obj, nplay_, fps_)
        if (nargin==1)
            nplay=5;
            fps=90;
        elseif (nargin==3)
            nplay=nplay_;
            fps=fps_;
        end
        obj.fig.Visible = 'on';
        movie(obj.ax, obj.movie_gen, nplay, fps);
    end
    function dims_out = get_dims(obj)
      curunits = get(obj.ax, 'Units');
      set(obj.ax, 'Units', 'Points');
      dims_out = get(obj.ax, 'Position');
      set(obj.ax, 'Units', curunits);
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
    function save_fig(fig_in, save_type, save_dir_)
      if (nargin==2)
        save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
      else
        save_dir = save_dir_;
      end
      if (strcmp(save_type, 'pdf'))
        exportgraphics(fig_in, [save_dir fig_in.Name '.pdf'], 'ContentType', 'vector')
      else
        saveas(fig_it, [save_dir fig_it.Name], save_type);
      end
    end
    function save_figs(figs, figs_to_write, save_type, save_dir_)
      if (nargin==2)
        save_dir = [getenv('HOME') '/Desktop/MATLAB_OUTPUT/'];
      else
        save_dir = save_dir_;
      end
      if (strcmp(save_type, 'pdf'))
        for i=reshape(figs_to_write, 1, [])
          fig_it = figs(i);
          exportgraphics(fig_it, [save_dir fig_it.Name '.pdf'], 'ContentType', 'vector');
        end
      else
        for i=reshape(figs_to_write, 1, [])
          fig_it = figs(i);
          saveas(fig_it, [save_dir fig_it.Name], save_type);
        end
      end
    end
  end
end
