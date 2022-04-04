classdef swirl < double
  properties
  % big brother
    % ODR_data.m;

    %% this class holds the double precision data of a simulation. It IS ODR_data.data. Note that it is NOT a handle

    %% data shared with swirl and
    specs=0; % (t cx cy wall_sca) for each frame
    Frames=0; % number of frames in swirl
    len_specs=0; % length of specs vector in each frame
    len_dat=0; % row length of data matrix; (position, quaternion)
    beads=0; % number of beads

    params;
    params_len=14;

    initialized = false;
  end
  methods(Static)
    function cell_out = tens2cell(tens_);
      cell_out = 0;
    end
  end
  methods
    function data = swirl(data_, specs_, Frames_, len_specs_, len_dat_, beads_)
      if (nargin==0)
        data = data@double(0);
      elseif(nargin<6)
        data=data@double(data_);
        data.beads = size(data_, 1);
        data.len_dat = size(data_, 2);
        data.Frames = size(data_, 3);
        if (nargin==2)
          data.len_specs = size(specs_, 2);
          data.specs=specs_;
        end
        data.initialized=true;
      end
      elseif (nargin==6)
        data = data@double(data_(:,:,1:Frames_));
        data.specs = specs_(:,:,1:Frames_);
        [data.Frames, data.len_specs, data.len_dat, data.beads] = deal(Frames_, len_specs_, len_dat_, beads_);
        data.initialized=true;
      end
    end
    function [out, valid_sub_swirl] = spawn_sub_swirl(obj, oth)
      %% on output, we want to give the sub swirl of the input that matches THIS swirl
      out = 0;
      valid_sub_swirl = false;
      [f,not_lost,fcap,tol] = deal(0, true, min(obj.Frames, oth.Frames), 1);
      diff_full = obj(:, 1:2, :)-oth(:, 1:2, :);
      while ((f<fcap)&&(not_lost))
        f=f+1;
        diff = diff_full(:, :, f);
        not_lost = max(sqrt(sum((diff.*diff)')))<tol;
      end
      f_final = f;
      if (f_final > 1)
        valid_sub_swirl = true;
        out = swirl(oth, oth.specs, f_final, oth.len_specs, oth.len_dat, oth.beads);
      end
    end
  end
end
