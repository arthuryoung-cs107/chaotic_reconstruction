classdef swirl < double
  properties
  % big brother
    % ODR_data.m;

    %% this class holds the double precision data of a simulation. It IS ODR_data data. Note that it is NOT a handle

    %% data shared with swirl ODR data
    Frames=0; % number of frames in swirl
    beads=0; % number of beads
    len_pos=0; % row length of position matrix; (position, quaternion)
    len_dish=0; % length of dish data in each frame
    len_params=14;
    % pos=0; % Note that a swirl IS this tensor
    dish=0; % (t cx cy wall_sca) for each frame
    params; % the parameters associated with the swirl

    initialized = false;
  end
  methods(Static)
    function mat_pos = tens2mat(tens_)
      [beads, len_dat, Frames] = size(tens);
      mat_pos = reshape(tens_, [beads*len_dat, Frames]);
    end
  end
  methods
    function data = swirl(pos_, dish_, Frames_, len_dish_, len_pos_, beads_)
      if (nargin==0)
        data = data@double(0);
      elseif(nargin<6)
        data=data@double(pos_);
        data.beads = size(pos_, 1);
        data.len_pos = size(pos_, 2);
        data.Frames = size(pos_, 3);
        if (nargin==2)
          data.len_dish = size(dish_, 2);
          data.dish=dish_;
        end
        data.initialized=true;
      end
      elseif (nargin==6)
        data = data@double(data_(:,:,1:Frames_));
        data.dish = dish_(:,:,1:Frames_);
        [data.Frames, data.len_dish, data.len_pos, data.beads] = deal(Frames_, len_dish_, len_pos_, beads_);
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
        out = swirl(oth, oth.dish, f_final, oth.len_dish, oth.len_pos, oth.beads);
      end
    end
  end
end
