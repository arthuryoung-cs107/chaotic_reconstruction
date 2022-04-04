classdef swirl_group < swirl
  properties
      group_cell;
      gp_len;
      frame_vec;
      params_mat;
  end
  methods(Static)
  end
  methods
    function obj = swirl_group(sw_, gp_len_, frame_vec_)
      obj=obj@swirl(sw_, sw_.specs); % in memory, will correspond to a order
      if (nargin>1)
        obj.gp_len = gp_len_;
        obj.group_cell = cell([gp_len_, 2]);
        obj.params_mat = nan(gp_len_, obj.params_len); 
        if (nargin==2)
          [Frames, len_dat, len_specs, beads] = deal(obj.Frames, obj.len_dat, obj.len_specs, obj.beads);
          obj.frame_vec = Frames*ones(gp_len_, 1);
          obj.group_cell(:, 1) = deal(nan(beads*len_dat, Frames));
          obj.group_cell(:, 2) = deal(nan(len_specs, Frames));
        elseif (nargin==3)
          obj.frame_vec = frame_vec_;
          dim = size(sw_);
        end
      end
    end
    function make_moviei(obj, AYfig_in, i )
      obj.sw(i).make_movie(AYfig_in);
    end
    function make_movieij(obj, AYfig_in, i, j)
      obj.sw(i).make_movie_comp(AYfig_in, obj.sw(j));
    end
    function make_movieijk(obj, AYfig_in, i, j, k)
      obj.sw(i).make_movie_comp2(AYfig_in, obj.sw(j), obj.sw(k));
    end
  end
end
