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
    contact_frames;
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
    function err_vec = comp_pos_err(obj, oth)
      %% computes mean bead position error for each frame
      err_vec = zeros(1, obj.Frames-1); %% ignore first frame, assuming equivalent initial conditions
      for i=1:(obj.Frames-1)
        diff_mat = obj.data(:, 1:3, i+1) - oth.data(:, 1:3, i+1);
        err_vec(i) = sum(sum((diff_mat.*diff_mat)'));
      end
    end
    function find_contact_frames(obj)
      obj.contact_frames = zeros(obj.Frames, 1);
      fa = sqrt(0.75);
      d = 5.72;
      walls = [0, 1; fa, 0.5; fa , -0.5]';
      for i=1:obj.Frames
        s_mat = (obj.data(:, 1:2, i)-[obj.specs(i, 2), obj.specs(i, 3)])*walls;
        w_mat = s_mat - d*(double(s_mat>0)) + d*(double(s_mat<0));
        contact_beads = find(abs(w_mat)<0.5); %% rycrofts contact criterion
        obj.contact_frames(i) = length(contact_beads)>0;
      end
    end
    function make_movie(obj, AYfig_in)
      frames = obj.Frames;

      walld = 5.72;
      wallL = (2/sqrt(3))*walld;
      wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld ; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];
      AYfig_in.init_movie(frames);
      wsca = max(abs([min(obj.specs(:, 2))-wallL, max(obj.specs(:, 2))+walld, min(obj.specs(:, 3))-wallL, max(obj.specs(:, 3))+walld]));
      lims = [-wsca, wsca, -wsca, wsca];

      rads = 0.5*ones(obj.P,1);
      figdims = AYfig_in.get_dims();
      MS = figdims(4)/(4*wsca); %% radius in pixels
      SS = pi*MS*MS; %% area in pixels
      for i=1:frames
        plot(AYfig_in.ax, wallv(:, 1)+obj.specs(i, 2), wallv(:, 2)+obj.specs(i, 3), 'k -')
        hold(AYfig_in.ax, 'on');
        objdots = scatter(AYfig_in.ax, obj.data(:, 1, i), obj.data(:, 2, i), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1000 0.4440 0.2440], 'MarkerFaceAlpha', 0.5);
        txtbx = annotation(AYfig_in.fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(i-1), 'LineStyle', 'none', 'FontSize', 16);
        hold(AYfig_in.ax, 'off');
        axis(AYfig_in.ax, lims);
        drawnow
        AYfig_in.movie_gen(i) = getframe(AYfig_in.ax);
        delete(txtbx);
      end
      % alternative way of plotting the circles, but not quite as flexible
      % viscircles(AYfig_in.ax, obj.data(:, 1:2, i), rads, 'Color', [0.1000 0.4440 0.2440]);
      % objdots = plot(AYfig_in.ax, obj.data(:, 1, i), obj.data(:, 2, i), 'o', 'ColorMode', 'manual', 'LineWidth', 1, 'MarkerSize', MS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1000 0.4440 0.2440]);
    end
    function make_movie_comp(obj, AYfig_in, oth)
      frames = obj.Frames;
      walld = 5.72;
      wallL = (2/sqrt(3))*walld;
      wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld ; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];
      AYfig_in.init_movie(frames);
      wsca = max(abs([min(obj.specs(:, 2))-wallL, max(obj.specs(:, 2))+walld, min(obj.specs(:, 3))-wallL, max(obj.specs(:, 3))+walld]));
      lims = [-wsca, wsca, -wsca, wsca];

      rads = 0.5*ones(obj.P,1);
      figdims = AYfig_in.get_dims();
      MS = figdims(4)/(4*wsca); %% radius in pixels
      SS = pi*MS*MS; %% area in pixels
      for i=1:frames
        plot(AYfig_in.ax, wallv(:, 1)+obj.specs(i, 2), wallv(:, 2)+obj.specs(i, 3), 'k -')
        hold(AYfig_in.ax, 'on');
        objdots = scatter(AYfig_in.ax, obj.data(:, 1, i), obj.data(:, 2, i), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1000 0.4440 0.2440], 'MarkerFaceAlpha', 0.5);
        othdots = scatter(AYfig_in.ax, oth.data(:, 1, i), oth.data(:, 2, i), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerFaceAlpha', 0.5);
        txtbx = annotation(AYfig_in.fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(i-1), 'LineStyle', 'none');
        txtbx.FontSize = 16;
        hold(AYfig_in.ax, 'off');
        axis(AYfig_in.ax, lims);

        drawnow

        AYfig_in.movie_gen(i) = getframe(AYfig_in.ax);
        delete(txtbx);
      end
    end
    function make_movie_comp2(obj, AYfig_in, oth, odd)
      frames = obj.Frames;
      walld = 5.72;
      wallL = (2/sqrt(3))*walld;
      wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld ; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];
      AYfig_in.init_movie(frames);
      wsca = max([abs(walld) abs(wallL)]);
      lims = [-wsca, wsca, -wsca, wsca];

      rads = 0.5*ones(obj.P,1);
      figdims = AYfig_in.get_dims();
      MS = figdims(4)/(4*wsca); %% radius in pixels
      SS = pi*MS*MS; %% area in pixels
      for i=1:frames
        plot(AYfig_in.ax, wallv(:, 1), wallv(:, 2), 'k -')
        hold(AYfig_in.ax, 'on');
        objdots = scatter(AYfig_in.ax, obj.data(:,1,i)-obj.specs(i, 2), obj.data(:,2,i)-obj.specs(i, 3), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1000 0.4440 0.2440], 'MarkerFaceAlpha', 0.5);
        othdots = scatter(AYfig_in.ax, oth.data(:,1,i)-oth.specs(i, 2), oth.data(:,2,i)-oth.specs(i, 3), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerFaceAlpha', 0.5);
        odddots = scatter(AYfig_in.ax, odd.data(:,1,i)-odd.specs(i, 2), odd.data(:,2,i)-odd.specs(i, 3), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'MarkerFaceAlpha', 0.5);
        txtbx = annotation(AYfig_in.fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(i-1), 'LineStyle', 'none');
        txtbx.FontSize = 16;
        hold(AYfig_in.ax, 'off');
        axis(AYfig_in.ax, lims);

        drawnow
        AYfig_in.movie_gen(i) = getframe(AYfig_in.ax);
        delete(txtbx);
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
    function data_out = filpos2swpos(obj)
      data_in = obj.filin.pos;
      dims = size(data_in);
      data_out = nan(dims);
      for i=1:dims(3)
        data_out(:, 1, i) = (data_in(:, 1, i)-obj.filin.cx_im)/obj.filin.cl_im + obj.specs(i, 2);
        data_out(:, 2, i) = (data_in(:, 2, i)-obj.filin.cy_im)/obj.filin.cl_im + obj.specs(i, 3);
      end
    end
    function sw_out = spawn_best(obj)
      sw_out = ODR_data([obj.dat_dir_name obj.exp_name], 'swirl_best.odr/', obj.dat_name);
    end
  end

end
