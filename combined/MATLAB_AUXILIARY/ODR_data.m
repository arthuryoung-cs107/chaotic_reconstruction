classdef ODR_data < handle
  properties
    dat_dir_name;
    exp_name;
    dat_name;

    Frames;
    len_specs;
    len_dat;
    P;

    rysml;
    filin;
    specs;
    data;
  end

  methods
    function obj = ODR_data(dat_dir_name_, exp_name_, dat_name_)
      obj.dat_dir_name = dat_dir_name_;
      obj.exp_name = exp_name_;
      obj.dat_name = dat_name_;
      obj.rysml = dlmread([obj.dat_dir_name obj.exp_name obj.dat_name '.rysml']);
      obj.Frames = obj.rysml(1, 2);
      obj.len_specs = obj.rysml(2, 3);
      obj.len_dat = obj.rysml(3, 3);
      obj.P = obj.rysml(3, 2);

      obj.specs = nan(obj.Frames, obj.len_specs);
      obj.data = nan(obj.P, obj.len_dat, obj.Frames);
      if (obj.rysml(1,1)==1)
        for f=0:(obj.Frames-1)
          id = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.' num2str(f) '.aydat']);
          obj.specs(f+1, :) = (fread( id,[1, obj.len_specs], 'float64=>float64'));
          obj.data(:, :, f+1) = (fread( id,[obj.len_dat, obj.P], 'float64=>float64'))';
          fclose(id);
        end
      elseif (obj.rysml(1,1)==0)
        id = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.rydat']);
        for f=0:(obj.Frames-1)
          obj.specs(f+1, :) = (fread( id,[1, obj.len_specs], 'float64=>float64'));
          obj.data(:, :, f+1) = (fread( id,[obj.len_dat, obj.P], 'float64=>float64'))';
        end
        fclose(id);
      end
    end
    function load_filin(obj)
      obj.filin = filter_inputs([obj.dat_dir_name obj.exp_name obj.dat_name]);
    end
    function verify_filin(obj)
      t_phys = obj.filin.vidspecs(1);
      cx_im = obj.filin.vidspecs(2);
      cy_im = obj.filin.vidspecs(3);
      cl_im = obj.filin.vidspecs(4);

      % works
      % obj.specs(1:10, 1)'
      % obj.filin.t(1:10)/t_phys

      %%works
      % ind = 100;
      % cx_loc = obj.specs(ind, 2);
      % cy_loc = obj.specs(ind, 3);
      % x = obj.data(:, 1, ind);
      % y = obj.data(:, 2, ind);
      % [cx_im + cl_im*(x-cx_loc), cy_im + cl_im*(y-cy_loc)]
      % obj.filin.pos(:, :, ind)

    end
    function err_vec = comp_pos_err(obj, oth)
      %% computes mean bead position error for each frame
      err_vec = zeros(1, obj.Frames-1); %% ignore first frame, assuming equivalent initial conditions
      for i=1:(obj.Frames-1)
        % res = obj.data(:, 1:3, i+1) - oth.data(:, 1:3, i+1);
        % err_vec(i) = mean(sum((res.*res)'));
        err_vec(i) = norm(obj.data(:, 1:3, i+1) - oth.data(:, 1:3, i+1), 'fro')/norm(obj.data(:, 1:3, i+1), 'fro');
      end
    end
  end

end
