classdef record < handle
  properties
    global_index;
    frscore;
    gen;
    parent_gen;
    parent_count;

    l2score;

    params;
  end

  methods
    function obj = record(int_vec_in, double_vec_in, param_vec_in);
      obj.global_index = int_vec_in(1);
      obj.frscore = int_vec_in(2);
      obj.gen = int_vec_in(3);
      obj.parent_gen = int_vec_in(4);
      obj.parent_count = int_vec_in(5);

      obj.l2score = double_vec_in(1);

      obj.params = param_vec_in;
    end
  end

end
