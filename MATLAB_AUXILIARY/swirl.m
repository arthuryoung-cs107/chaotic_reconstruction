classdef swirl < handle
  properties
  % friend class
    odr; % ODR_data.m;

    data;

  end
  methods
    function obj = swirl(sub)

    end
    function init_ODR(obj, odr_)
      obj.odr = odr_;



    end
    function [out, valid_sub_swirl] = spawn_sub_swirl(obj, oth)
      %% on output, we want to give the sub swirl of the input that matches THIS swirl
      out = 0;
      valid_sub_swirl = false;
      [f,not_lost,fcap,tol] = deal(0, true, min(obj.Frames, oth.Frames), 1);
      while ((f<=fcap)&&(not_lost))
        f=f+1;
        diff = obj.data(:, 1:2, f)-oth.data(:, 1:2, f);
        not_lost = max(sqrt(sum((diff.*diff)')))<tol;
      end
      f_final = f;
      if (f_final > 1)
        valid_sub_swirl = true;
        out = swirl();
        out.Frames = f_final;
        out.len_specs = oth.len_specs;
        out.data = oth.data(:, :, 1:f_final);
        out.P = oth.P;
      end
    end
  end
end
