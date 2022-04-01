classdef stat_data < handle
  properties
  % friend class
    sgp; % swirl_group.m;


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
    par_cov;
    pos_err_sum;
    pos_err_best;
    par_err_truest;
    I_best;
    I_truest;

  end

  methods
    function obj = stat_data(dat_dir_name_, exp_name_, dat_name_, noise_len_)
      obj.dat_dir_name = dat_dir_name_;
      obj.exp_name = exp_name_;
      obj.dat_name = dat_name_;

      obj.pars = AYdata.aysml_read([obj.dat_dir_name obj.exp_name obj.dat_name '_pars']);
      obj.dels = AYdata.aysml_read([obj.dat_dir_name obj.exp_name obj.dat_name '_dels']);
      if (noise_len_~=0)
        obj.pars = obj.pars(:, 1:(noise_len_+1));
      end
      obj.par_len = size(obj.pars, 1);
      obj.noise_len = size(obj.pars, 2);

      obj.sw = ODR_data.empty(obj.noise_len, 0);
      for i=0:(obj.noise_len-1)
        obj.sw(i+1) = ODR_data(obj.dat_dir_name, obj.exp_name, [obj.dat_name '.' num2str(i)]);
      end

      obj.par_err = nan(obj.noise_len-1, obj.par_len);
      obj.pos_err = nan(obj.noise_len-1, obj.sw(1).Frames-1);
      for i=1:(obj.noise_len-1)
        obj.par_err(i, :) = (obj.pars(:, 1)-obj.pars(:, i+1))./abs(obj.pars(:, 1));
        obj.pos_err(i, :) = obj.sw(1).comp_pos_err(obj.sw(i+1));
      end
      [obj.par_err_truest, obj.I_truest] = mink(sum(abs(obj.par_err)'), obj.noise_len-1);
      obj.pos_err_sum = sum(obj.pos_err');
      [obj.pos_err_best, obj.I_best ] = mink(obj.pos_err_sum, obj.noise_len-1);
      obj.par_cov = (abs(obj.par_err)).*(obj.pos_err_sum)';
    end
    function swout = spawn_swrlindex(obj, id_)
      swout = ODR_data(obj.dat_dir_name, obj.exp_name, [obj.dat_name '.' num2str(i)]);
    end
    function swgp_out = spawn_swirl_group(obj, )

    end
    function get_swirl_group(obj, )

    end
    function plot_frame_error(obj, fig_in, base_color)
      mean_err = mean(obj.pos_err);
      Frames = length(mean_err)+1;
      figure(fig_in.Number)
      for i=1:obj.noise_len-1
        plot(1:Frames-1, obj.pos_err(i, :), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
      end
      plot(1:Frames-1, mean_err, ' -', 'Color', base_color, 'LineWidth', 2);
      plot(1:Frames-1, obj.pos_err(obj.I_best(1), :), ' -', 'Color', [0 0 0], 'LineWidth', 1);
      plot(1:Frames-1, obj.pos_err(obj.I_truest(1), :), ' :', 'Color', [0 0 0], 'LineWidth', 1);
    end
    function plot_param_error(obj, fig_in, base_color)
      figure(fig_in.Number)
      for i=1:obj.noise_len-1
        plot(1:obj.par_len, (obj.pars(:, 1)-obj.pars(:, i+1))./abs(obj.pars(:, 1)), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
      end
      plot(1:obj.par_len, (obj.pars(:, 1)-obj.pars(:, obj.I_best(1)+1))./abs(obj.pars(:, 1)), ' -', 'Color', [0 0 0], 'LineWidth', 1);
      plot(1:obj.par_len, (obj.pars(:, 1)-obj.pars(:, obj.I_truest(1)+1))./abs(obj.pars(:, 1)), ' :', 'Color', [0 0 0], 'LineWidth', 1);
    end
    function plot_param_index_error(obj, fig_in, base_color)
      figure(fig_in.Number)
      for i=1:obj.noise_len-1
        plot(1:obj.par_len, obj.par_cov(i, :), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
      end
      plot(1:obj.par_len, mean(obj.par_cov), ' -', 'Color', base_color, 'LineWidth', 2);
      plot(1:obj.par_len, obj.par_cov(obj.I_best(1), :), ' -', 'Color', [0 0 0 ], 'LineWidth', 1);
      plot(1:obj.par_len, obj.par_cov(obj.I_truest(1), :), ' :', 'Color', [0 0 0 ], 'LineWidth', 1);
    end
    function plot_param_pos_error(obj, fig_in, base_color)
      figure(fig_in.Number)
      plot(sum(abs(obj.par_err)'), obj.pos_err_sum, ' o', 'Color', base_color, 'LineWidth', 2);
      plot(sum(abs(obj.par_err(obj.I_best(1), :))), obj.pos_err_sum(obj.I_best(1)), ' p', 'Color', [0 0 0], 'LineWidth', 2);
      plot(sum(abs(obj.par_err(obj.I_truest(1), :))), obj.pos_err_sum(obj.I_truest(1)), ' h', 'Color', [0 0 0], 'LineWidth', 2);
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
    function make_POVray_inputsi(obj, i_, pov_dir_)
      obj.sw(i_).make_POVray_inputs(pov_dir_, obj.dat_name);
    end
  end
end
