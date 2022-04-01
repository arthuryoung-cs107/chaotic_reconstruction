classdef swirl_group < swirl
  properties
      group_cell;
      gp_len;
      frame_vec;
  end
  methods(Static)
  end
  methods
    function obj = swirl_group(sw_, gp_len_, frame_num_)
      obj=obj@swirl(sw_); % in memory, will correspond to a order
      if (nargin>1)
        obj.gp_len = gp_len_;
        obj.group_cell = cell([gp_len_, 3]);
        if (nargin==3)
          obj.frame_vec = frame_num.*ones(gp_len_, 1);
          dim = size(sw_);
          obj.group_cell(:, 1) = deal(nan(dim(1)*dim(2), frame_num_));
          if (sw_.initialized)
          obj.group_cell(:, 2) = deal(dim(3)),);
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
