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

    end
    function err_vec = comp_pos_err(obj, oth)
      %% computes mean bead position error for each frame
      err_vec = zeros(1, obj.Frames-1); %% ignore first frame, assuming equivalent initial conditions
      for i=1:(obj.Frames-1)
        diff_mat = obj.data(:, 1:3, i+1) - oth.data(:, 1:3, i+1);
        err_vec(i) = sum(sum((diff_mat.*diff_mat)'));
        % err_vec(i) = norm(obj.data(:, 1:3, i+1) - oth.data(:, 1:3, i+1), 'fro')/norm(obj.data(:, 1:3, i+1), 'fro');
      end
    end
    function make_movie(obj, AYfig_in)
      frames = obj.Frames;
      % frames = 200;

      walld = 5.72;
      wallL = (2/sqrt(3))*walld;
      wallv = [-wallL/2 -walld; -walld 0; -wallL/2 walld ; wallL/2 walld; walld 0; wallL/2 -walld; -wallL/2 -walld];

      AYfig_in.init_movie(frames);
      lims = [-walld, walld, -walld, walld];

      figdims = AYfig_in.get_dims();
      r = 0.5;
      MS = figdims(4)/(4*walld);

      for i=1:frames
        plot(AYfig_in.ax, wallv(:, 1), wallv(:, 2), 'k -')
        hold(AYfig_in.ax, 'on');
        plot(AYfig_in.ax, obj.data(:, 1, i)-obj.specs(i, 2), obj.data(:, 2, i)-obj.specs(i, 3), 'o', 'Color', [0 0 1], 'LineWidth', 1.0, 'MarkerSize', MS);
        hold(AYfig_in.ax, 'off');
        axis(AYfig_in.ax, lims);
        drawnow

        AYfig_in.movie_gen(i) = getframe(AYfig_in.ax);
      end
    end
    function make_movie_comp(obj, AYfig_in, oth)
      frames = obj.Frames;
      % frames = 200;

      walld = 5.72;
      wallL = (2/sqrt(3))*walld;
      wallv = [-wallL/2 -walld; -walld 0; -wallL/2 walld ; wallL/2 walld; walld 0; wallL/2 -walld; -wallL/2 -walld];
      AYfig_in.init_movie(frames);
      lims = [-walld, walld, -walld, walld];

      figdims = AYfig_in.get_dims();
      r = 0.5;
      MS = figdims(4)/(4*walld);

      for i=1:frames
        plot(AYfig_in.ax, wallv(:, 1), wallv(:, 2), 'k -')
        hold(AYfig_in.ax, 'on');
        plot(AYfig_in.ax, obj.data(:, 1, i)-obj.specs(i, 2), obj.data(:, 2, i)-obj.specs(i, 3), 'o', 'Color', [0 0 1], 'LineWidth', 1.0, 'MarkerSize', MS);
        plot(AYfig_in.ax, oth.data(:, 1, i)-oth.specs(i, 2), oth.data(:, 2, i)-oth.specs(i, 3), 'o', 'Color', [1 0 0], 'LineWidth', 1.0, 'MarkerSize', MS);
        hold(AYfig_in.ax, 'off');
        axis(AYfig_in.ax, lims);
        drawnow

        AYfig_in.movie_gen(i) = getframe(AYfig_in.ax);
      end
    end
    function make_POVray_inputs(obj, pov_dir_, dat_name_)
      save_dir = [pov_dir_ obj.exp_name];
      status = mkdir(save_dir);
      for i=1:obj.Frames
        id = fopen([save_dir dat_name_ '.' num2str(i-1)], 'w');
        fprintf(id, '%g %g %g %g\n', obj.specs(i, :));
        fprintf(id, '%g %g %g %g %g %g %g\n', obj.data(:, :, i)');
        fclose(id);
      end
    end
  end

end
