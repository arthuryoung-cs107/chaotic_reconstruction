classdef swirl_group < handle
  properties
    sw;
    sw_len;
  end
  methods
    function obj = swirl_group(sw_len_)
      obj.sw_len = sw_len_;
      obj.sw = swirl.empty(obj.sw_len, 0);
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
