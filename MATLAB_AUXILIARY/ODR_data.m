classdef ODR_data < handle
  properties
  % little brother:
    % swirl.m
    dat_dir_name;
    exp_name;
    dat_name;

    filin;

    %% data shared with swirl
    Frames=0; % number of frames in swirl
    beads=0; % number of beads
    len_pos=0; % row length of position matrix; (position, quaternion)
    len_dish=0; % length of dish data in each frame
    len_params=14; % length of swirl params, all inclusive
    pos=0 ; % Note that a swirl IS this tensor
    dish=0; % (t cx cy wall_sca) for each frame
    params; % the parameters associated with the swirl (swirl_params)

    contact_frames;
  end
  methods(Static)
    function sw = construct_swirl(dat_dir_name_, exp_name_, dat_name_)
      rysml = ODR_data.get_rysml(dat_dir_name_, exp_name_, dat_name_);
      rydat = ODR_data.get_rydat(rysml);
      sw = swirl(rydat.pos, rydat.dish, rysml.Frames, rysml.len_dish, rysml.len_pos, rysml.beads);
    end
    function rysml_out = get_rysml(dat_dir_name_,exp_name_,dat_name_)
      rysml = dlmread([dat_dir_name_ exp_name_ dat_name_ '.rysml']);
      rysml_out = struct('dat_dir', dat_dir_name_, 'exp_name',exp_name_,'dat_name',dat_name_,'split',rysml(1,1),'Frames',rysml(1,2),'beads',rysml(3,2),'len_dish',rysml(2, 3),'len_pos',rysml(3,3));
    end
    function rydat_out = get_rydat(rysml)
        dish = nan(rysml.Frames, rysml.len_dish);
        pos = nan(rysml.beads, rysml.len_pos, rysml.Frames);
        if (rysml.split==1) %% if the data is stored across many binary files
            for f=0:(rysml.Frames-1)
                id = fopen([rysml.dat_dir rysml.exp_name rysml.dat_name '.' num2str(f) '.aydat']);
                dish(f+1, :) = (fread( id,[1, rysml.len_dish], 'float64=>float64'));
                pos(:, :, f+1) = (fread( id,[rysml.len_pos, rysml.beads], 'float64=>float64'))';
                fclose(id);
            end
        elseif (rysml.split==0) %% if the data is consolidated into one binary file
            id = fopen([rysml.dat_dir rysml.exp_name rysml.dat_name '.rydat']);
            for f=0:(rysml.Frames-1)
                dish(f+1, :) = fread( id,[1, rysml.len_dish], 'float64=>float64');
                pos(:, :, f+1) = (fread( id,[rysml.len_pos, rysml.beads], 'float64=>float64'))';
            end
            fclose(id);
        end
        rydat_out = struct('pos', pos, 'dish', dish);
    end
  end
  methods
    function obj = ODR_data(dat_dir_name_, exp_name_, dat_name_)
      %% first, construct the swirl that represents this data simulation outputs
      rysml = ODR_data.get_rysml(dat_dir_name_, exp_name_, dat_name_);
      rydat = ODR_data.get_rydat(rysml);

      %% ODR member assignments
      obj.dat_dir_name = rysml.dat_dir;
      obj.exp_name = rysml.exp_name;
      obj.dat_name = rysml.dat_name;

      obj.Frames=rysml.Frames;
      obj.beads=rysml.beads;
      obj.len_pos=rysml.len_pos;
      obj.len_dish=rysml.len_dish;
      % obj.len_params=rysml.len_params;
      obj.pos=rydat.pos;
      obj.dish=rydat.dish;
      % obj.params=rydat.params;
    end
    function sw = spawn_swirl(obj)
      sw = swirl(obj.pos, obj.dish, obj.Frames, obj.len_dish, obj.len_pos, obj.beads);
    end
    function load_filin(obj)
      obj.filin = filter_inputs([obj.dat_dir_name obj.exp_name obj.dat_name]);
    end
    function gen_out = read_geni(obj, gen_count_)
      gen_out = generation(obj.dat_dir_name, obj.exp_name, obj.dat_name, gen_count_);
    end
    function gen_out = read_race_data(obj)
      gen_out = race_data(obj.dat_dir_name, obj.exp_name, obj.dat_name);
    end
    function make_POVray_inputs(obj, pov_dir_, dat_name_)
      save_dir = [pov_dir_ obj.exp_name];
      status = mkdir(save_dir);
      for i=1:obj.Frames
        id = fopen([save_dir dat_name_ '.' num2str(i-1)], 'w');
        fprintf(id, '%g %g %g %g\n', obj.dish(i, :));
        fprintf(id, '%g %g %g %g %g %g %g\n', obj.pos(:, :, i)');
        fclose(id);
      end
    end
    function sw_out = spawn_best(obj)
      sw_out = ODR_data([obj.dat_dir_name obj.exp_name], 'swirl_best.odr/', obj.dat_name);
    end
    function make_movie_comp(obj, AYfig_in, oth)
      frames = obj.Frames;
      walld = 5.72;
      wallL = (2/sqrt(3))*walld;
      wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld ; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];
      AYfig_in.init_movie(frames);
      wsca = max(abs([min(obj.specs(:, 2))-wallL, max(obj.specs(:, 2))+walld, min(obj.specs(:, 3))-wallL, max(obj.specs(:, 3))+walld]));
      lims = [-wsca, wsca, -wsca, wsca];

      rads = 0.5*ones(obj.beads,1);
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

      rads = 0.5*ones(obj.beads,1);
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
    function make_movie_comp3(obj, AYfig_in, oth, odd, awk)
      frames = obj.Frames;
      run AYfigprops.m
      walld = 5.72;
      wallL = (2/sqrt(3))*walld;
      wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld ; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];
      AYfig_in.init_movie(frames);
      wsca = max(abs([min(obj.specs(:, 2))-wallL, max(obj.specs(:, 2))+walld, min(obj.specs(:, 3))-wallL, max(obj.specs(:, 3))+walld]));
      % wsca = max([abs(walld) abs(wallL)]);
      lims = [-wsca, wsca, -wsca, wsca];

      rads = 0.5*ones(obj.beads,1);
      figdims = AYfig_in.get_dims();
      MS = figdims(4)/(4*wsca); %% radius in pixels
      SS = pi*MS*MS; %% area in pixels
      for i=1:frames
        plot(AYfig_in.ax, wallv(:, 1)+obj.specs(i, 2), wallv(:, 2)+obj.specs(i, 3), 'k -')
        hold(AYfig_in.ax, 'on');
        objdots = scatter(AYfig_in.ax, obj.data(:,1,i), obj.data(:,2,i), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', green4, 'MarkerFaceAlpha', 0.5);
        othdots = scatter(AYfig_in.ax, oth.data(:,1,i), oth.data(:,2,i), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', red5, 'MarkerFaceAlpha', 0.5);
        odddots = scatter(AYfig_in.ax, odd.data(:,1,i), odd.data(:,2,i), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', blue4, 'MarkerFaceAlpha', 0.5);
        awkdots = scatter(AYfig_in.ax, awk.data(:,1,i), awk.data(:,2,i), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', orange1, 'MarkerFaceAlpha', 0.5);
        txtbx = annotation(AYfig_in.fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(i-1), 'LineStyle', 'none');
        txtbx.FontSize = 16;
        hold(AYfig_in.ax, 'off');
        axis(AYfig_in.ax, lims);

        drawnow
        AYfig_in.movie_gen(i) = getframe(AYfig_in.ax);
        delete(txtbx);
      end
    end
  end
end
