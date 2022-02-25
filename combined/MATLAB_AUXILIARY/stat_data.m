classdef stat_data < handle
  properties
    dat_dir_name;
    exp_name;
    dat_name;

    pars;
    dels;

    par_len;
    noise_len;

    sw;

    par_err;
    pos_err;
  end

  methods
    function obj = stat_data(dat_dir_name_, exp_name_, dat_name_)
      obj.dat_dir_name = dat_dir_name_;
      obj.exp_name = exp_name_;
      obj.dat_name = dat_name_;

      obj.pars = AYdata.aysml_read([obj.dat_dir_name obj.exp_name obj.dat_name '_pars']);
      obj.dels = AYdata.aysml_read([obj.dat_dir_name obj.exp_name obj.dat_name '_dels']);
      obj.par_len = size(obj.pars, 1);
      obj.noise_len = size(obj.pars, 2);

      obj.sw = ODR_data.empty(obj.noise_len, 0);
      for i=0:(obj.noise_len-1)
        obj.sw(i+1) = ODR_data(obj.dat_dir_name, obj.exp_name, [obj.dat_name '.' num2str(i)]);
      end

      obj.par_err = nan(obj.noise_len-1, 1);
      obj.pos_err = nan(obj.noise_len-1, obj.sw(1).Frames-1);
      for i=1:(obj.noise_len-1)
        obj.par_err(i) = norm((obj.pars(:, 1)-obj.pars(:, i+1))./obj.pars(:, 1));
        obj.pos_err(i, :) = obj.sw(1).comp_pos_err(obj.sw(i+1));
      end
    end
    function plot_frame_error(obj, fig_in, base_color)
      mean_err = mean(obj.pos_err);
      Frames = length(mean_err)+1;
      figure(fig_in.Number)
      for i=1:obj.noise_len-1
        plot(1:Frames-1, obj.pos_err(i, :), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
      end
      plot(1:Frames-1, mean_err, ' -', 'Color', base_color, 'LineWidth', 2);
    end
    function plot_param_error(obj, fig_in, base_color)
      figure(fig_in.Number)
      for i=1:obj.noise_len-1
        plot(1:obj.par_len, (obj.pars(:, 1)-obj.pars(:, i+1))./abs(obj.pars(:, 1)), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
      end

    end
    function plot_param_frame_error(obj, fig_in, base_color)
      Frames = size(obj.pos_err,2 );
      % Frames = 50;
      al = 0.1;
      af = (al-1)/(1-(Frames+1));
      ab = al - af;
      [err_it, I_it] = mink(obj.par_err, length(obj.par_err));
      figure(fig_in.Number)
      for i=1:Frames
        plot(err_it, obj.pos_err(I_it, i), ' -', 'Color', [base_color, ab+(i*af)], 'LineWidth', 2);
      end

    end

  end
end
