classdef swirl_group < swirl
%{A class for the management of multiple swirl instances, each being related
to each other by ONE swirl that each of them compare to.%}
  properties
      %{a cell object that contains each swirl's position data and dish data%}
      gp_cell;
      %{the number of swirl instances collected by this swirl group%}
      gp_len;
      %{the frame duration of each swirl in this swirl group%}
      frame_vec;
      %{the parameters associated with each swirl in this swirl group%}
      params_mat;
  end
  methods(Static)
    
  end
  methods
    function obj = swirl_group(sw_, gp_len_)
      %{ the swirl that each member of the group is compared to is the identity element %}
      obj=obj@swirl(sw_, sw_.specs);
      if (nargin>1) %% if we are also given the length of the swirl group
        obj.gp_len = gp_len_;
        obj.gp_cell = cell([gp_len_, 2]);
        obj.params_mat = nan(gp_len_, obj.params_len);
        obj.frame_vec = nan(gp_len_, 1);
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
