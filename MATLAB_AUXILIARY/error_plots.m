classdef error_plots < swirl_plots
    properties

    end
    methods
        function obj = error_plots()

        end
    end
    methods (Static)
        function plot_frame_stats(ax1_, ax2_, frames_, frame_stats_)
            errorbar(ax1_, frames_', frame_stats_(1,:)', frame_stats_(2,:)', ' p -', 'Color', [0 0 0], 'LineWidth', 2);
            errorbar(ax2_, frames_', frame_stats_(3,:)', frame_stats_(4,:)', ' p -', 'Color', [0 0 0], 'LineWidth', 2);
        end
        function plot_res_vs_frames(ax_, base_color, st)
            [x_vec, pos_res, I_best, I_truest, I_leader] = deal(1:st.sw0.Frames, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, x_vec, pos_res, struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'Frames', 'position residual', 'plot_res_vs_frames');
            % set(ax_, 'YScale', 'log');
        end
        function plot_accres_vs_frames(ax_, base_color, st)
            [x_vec, pos_res, I_best, I_truest, I_leader] = deal(1:st.sw0.Frames, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, x_vec, cumsum(pos_res,2), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'Frames', 'accumulated position residual', 'plot_accres_vs_frames');
            % set(ax_, 'YScale', 'log');
        end
        function plot_delres_vs_frames(ax_, base_color, st)
            [x_vec, pos_res, I_best, I_truest, I_leader] = deal(1:st.sw0.Frames, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, x_vec(2:end), pos_res(:,2:end)-pos_res(:,1:(size(pos_res,2)-1)), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'Frames', 'change in position residual', 'plot_delres_vs_frames');
            % set(ax_, 'YScale', 'log');
        end
        function plot_accres_vs_res(ax_, base_color, st)
            [X_mat, pos_res, I_best, I_truest, I_leader] = deal(cumsum(st.pos_res,2), st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_matmat_frame_data(ax_, X_mat, pos_res, struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'accumulated position error', 'position error', 'plot_accres_vs_res');
            set(ax_, 'XScale', 'log');
            set(ax_, 'YScale', 'log');
        end
        function plot_INTres_vs_res(ax_, base_color, st)
            [t_vec, pos_res, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            INT_mat = pos_res(:,2:end).*(t_vec(2:end)-t_vec(1:(length(t_vec)-1)))';
            plot_matmat_frame_data(ax_, cumsum(INT_mat,2), pos_res(:,2:end), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'integral of residual over time', 'residual', 'plot_INTres_vs_res');
            set(ax_, 'XScale', 'log');
            set(ax_, 'YScale', 'log');
        end
        function plot_INTres_vs_time(ax_, base_color, st)
            [t_vec, pos_res, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            INT_mat = pos_res(:,2:end).*(t_vec(2:end)-t_vec(1:(length(t_vec)-1)))';
            plot_generic_frame_data(ax_, t_vec(2:end), cumsum(INT_mat,2), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'time', 'integral of residual over time', 'plot_INTres_vs_time');
            % set(ax_, 'YScale', 'log');
        end
        function plot_Dres_vs_time(ax_, base_color, st)
            [t_vec, pos_res, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            D_mat = (pos_res(:,2:end)-pos_res(:,1:(size(pos_res,2)-1)))./(t_vec(2:end)-t_vec(1:(length(t_vec)-1)))';
            plot_generic_frame_data(ax_, t_vec(2:end),D_mat,struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'time', 'derivative of position residual over time', 'plot_Dres_vs_time');
            % set(ax_, 'YScale', 'log');
        end
        function plot_err_vs_frames(ax_, base_color, st)
            [x_vec, pos_err, I_best, I_truest, I_leader] = deal(1:st.sw0.Frames, st.pos_err, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, x_vec, pos_err, struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'Frames', 'position error', 'plot_err_vs_frames');
            % set(ax_, 'YScale', 'log');
        end
        function plot_accerr_vs_frames(ax_, base_color, st)
            [x_vec, pos_err, I_best, I_truest, I_leader] = deal(1:st.sw0.Frames, st.pos_err, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, x_vec, cumsum(pos_err,2), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'Frames', 'accumulated position error', 'plot_accerr_vs_frames');
            % set(ax_, 'YScale', 'log');
        end
        function plot_delerr_vs_frames(ax_, base_color, st)
            [x_vec, pos_err, I_best, I_truest, I_leader] = deal(1:st.sw0.Frames, st.pos_err, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, x_vec(2:end), pos_err(:,2:end)-pos_err(:,1:(size(pos_err,2)-1)),struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'Frames', 'change in position error', 'plot_delerr_vs_frames');
            % set(ax_, 'YScale', 'log');
        end
        function plot_alpha_vs_INTres(ax_, base_color, st)
            [t_vec, pos_res, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            INT_mat = pos_res(:,2:end).*(t_vec(2:end)-t_vec(1:(length(t_vec)-1)))';
            plot_matmat_frame_data(ax_, cumsum(INT_mat,2), pos_res(:,2:end), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'integral of residual over time', 'residual', 'plot_accres_vs_res');
            set(ax_, 'XScale', 'log');
            set(ax_, 'YScale', 'log');
        end
    end
end

function plot_generic_frame_data(ax_, x, Y, plotspecs_)
    [len_gp, Frames] = size(Y);
    box(ax_,'on');
    hold(ax_, 'on');
    for i=1:len_gp
        plot(ax_,x, Y(i,:), ' -', 'Color', [plotspecs_.base_color, 0.1], 'LineWidth', 1);
    end
    plot(ax_, x, Y(plotspecs_.I1, :), ' -', 'Color', [0 0 0], 'LineWidth', 2);
    plot(ax_, x, Y(plotspecs_.I2, :), ' :', 'Color', [0 0 0], 'LineWidth', 2);
    plot(ax_, x, Y(plotspecs_.I3, :), ' -.', 'Color', [0 0 0], 'LineWidth', 2);
    xlim(ax_, [min(x), max(x)]);
    ylim(ax_, [min(min(Y)), max(max(Y))]);
end
function plot_matmat_frame_data(ax_, X, Y, plotspecs_)
    [len_gp, Frames] = size(Y);
    box(ax_,'on');
    hold(ax_, 'on');
    for i=1:len_gp
        plot(ax_, X(i,:), Y(i,:), ' -', 'Color', [plotspecs_.base_color, 0.1], 'LineWidth', 1);
    end
    plot(ax_, X(plotspecs_.I1,:), Y(plotspecs_.I1, :), ' -', 'Color', [0 0 0], 'LineWidth', 2);
    plot(ax_, X(plotspecs_.I2,:), Y(plotspecs_.I2, :), ' :', 'Color', [0 0 0], 'LineWidth', 2);
    plot(ax_, X(plotspecs_.I3,:), Y(plotspecs_.I3, :), ' -.', 'Color', [0 0 0], 'LineWidth', 2);
    xlim(ax_, [min(min(X)), max(max(X))]);
    ylim(ax_, [min(min(Y)), max(max(Y))]);
end
function apply_plot_labels(ax_, xname, yname, titlename)
    xlabel(ax_, xname, 'Interpreter', 'Latex', 'Fontsize', 14)
    ylabel(ax_, yname, 'Interpreter', 'Latex', 'Fontsize', 14)
    title(ax_, ['\textbf{', strrep(titlename, '_', ' ') , '}'], 'Interpreter', 'Latex', 'Fontsize', 14)
end

function prime_vec = approx_deriv_weighted_central(t_in, x_in)
  n = length(t_in);
  prime_vec = nan(size(x_in));

  k = 5; %% number of points considered. Must be odd
  l = (k-1)/2; %% number of points to left and right
  p = l+1; %% index of central point

  x = reshape(x_in(1:k), [1 k]);
  t = reshape(t_in(1:k), [1 k]);
  for i=1:l
    ind = [1:(i-1), i+1:k];
    t_hat = t(ind);
    x_hat = x(ind);
    w = (abs((0.5*(t_hat + t(i)))-t(i))).^(-1);
    m = ((x(i)-x_hat))./(t(i)-t_hat);
    prime_vec(i) = (sum(w.*m))/(sum(w));
  end
  for i=p:(n-l)
    x = reshape(x_in(i-l:i+l), [1 k]);
    t = reshape(t_in(i-l:i+l), [1 k]);
    ind = [1:p-1, p+1:k];
    t_hat = t(ind);
    x_hat = x(ind);
    w = (abs((0.5*(t_hat + t(p)))-t(p))).^(-1);
    m = ((x(p)-x_hat))./(t(p)-t_hat);
    prime_vec(i) = (sum(w.*m))/(sum(w));
  end
  x = reshape(x_in(n-k+1:n), [1 k]);
  t = reshape(t_in(n-k+1:n), [1 k]);
  for i=p+1:k
    ind = [1:i-1, i+1:k];
    t_hat = t(ind);
    x_hat = x(ind);
    w = (abs((0.5*(t_hat + t(i)))-t(i))).^(-1);
    m = ((x(i)-x_hat))./(t(i)-t_hat);
    prime_vec(n-k+i) = (sum(w.*m))/(sum(w));
  end
end
