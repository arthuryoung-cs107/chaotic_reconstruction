classdef error_plots < swirl_plots
    properties

    end
    methods
        function obj = error_plots()

        end
    end
    methods (Static)
        %% figure tools
        function lims_ = lims_minmax(data_)
            lims_ = [min(min(data_)), max(max(data_))];
        end
        function adjust_lims(axs_, xlim_, ylim_)
            if (length(axs_) > 1)
                for i = 1:length(axs_)
                    axis(axs_(i), [xlim_(i,:), ylim_(i,:)]);
                end
            else
                axis(axs_, [xlim_, ylim_]);
            end
        end

        %% alpha plots
        function plot_alphaINTres_vs_time(ax_, base_color, st)
            [t_vec, alpha_INTres, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, swirl_group.compute_INT(st.sw0.t_vec,st.pos_res), st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, t_vec, swirl_group.compute_alpha_logD(t_vec,alpha_INTres), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'time', '$$\alpha$$', 'plot_alphaINTres_vs_time');
            % set(ax_, 'YScale', 'log');
        end
        function plot_alphaINTres_vs_time_bead(axs_, base_color, st)
            [t_vec, bead_res_cell, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, st.pos_res_bead, st.I_best(1), st.I_truest(1), st.I_leader(1));

            for i = 1:st.sw0.beads
                posresINTi = swirl_group.compute_INT(t_vec,bead_res_cell{i});
                plot_generic_frame_data(axs_(i), t_vec, swirl_group.compute_alpha_logD(t_vec,posresINTi), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
                apply_plot_labels(axs_(i), 'time', '$$\alpha$$', ['plot_alphaINTres_vs_time_bead' num2str(i)]);
            end

            % set(ax_, 'YScale', 'log');
        end

        %% error bar plots
        function plot_residual_stats(ax1_, ax2_, ax3_, t_vec_, frames_, frame_stats_)
            errorbar(ax1_, t_vec_(frames_), frame_stats_(1,:), frame_stats_(2,:), ' p -', 'Color', [0 0 0], 'LineWidth', 2);
            errorbar(ax2_, frames_, frame_stats_(3,:), frame_stats_(4,:), ' p -', 'Color', [0 0 0], 'LineWidth', 2);
            errorbar(ax3_, t_vec_(frames_), frame_stats_(5,:), frame_stats_(6,:), ' p -', 'Color', [0 0 0], 'LineWidth', 2);
        end


        %% derivative data against data for residuals
        function plot_accres_vs_res(ax_, base_color, st)
            [pos_accres, pos_res, I_best, I_truest, I_leader] = deal(cumsum(st.pos_res,2), st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_matmat_frame_data(ax_, pos_accres, pos_res, struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'accumulated residual', 'position residual', 'plot_accres_vs_res');
            set(ax_, 'XScale', 'log');
            set(ax_, 'YScale', 'log');
        end
        function plot_Daccres_vs_accres(ax_, base_color, st)
            [t_vec, pos_accres, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, cumsum(st.pos_res), st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_matmat_frame_data(ax_, pos_accres, swirl_group.compute_D(t_vec, pos_accres), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'accumulated residual', 'derivative of accumulated residual', 'plot_Daccres_vs_accres');
            set(ax_, 'XScale', 'log');
            set(ax_, 'YScale', 'log');
        end
        function plot_INTres_vs_res(ax_, base_color, st)
            [t_vec, pos_res, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_matmat_frame_data(ax_, swirl_group.compute_INT(t_vec, pos_res), pos_res, struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'integral of residual over time', 'residual', 'plot_INTres_vs_res');
            set(ax_, 'XScale', 'log');
            set(ax_, 'YScale', 'log');
        end

        %% derivative/integrals of residual data against time
        function plot_INTres_vs_time(ax_, base_color, st)
            [t_vec, pos_res, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, t_vec, swirl_group.compute_INT(t_vec, pos_res), struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'time', 'integral of residual over time', 'plot_INTres_vs_time');
            % set(ax_, 'YScale', 'log');
        end
        function plot_Dres_vs_time(ax_, base_color, st)
            [t_vec, pos_res, I_best, I_truest, I_leader] = deal(st.sw0.t_vec, st.pos_res, st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, t_vec, swirl_group.compute_D(t_vec, pos_res),struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'time', 'derivative of position residual over time', 'plot_Dres_vs_time');
            % set(ax_, 'YScale', 'log');
        end

        %% residual data against frames
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

        %% error data against frames
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
        function plot_param_error(ax_, base_color, st)
            [x_vec, Y_mat, I_best, I_truest, I_leader] = deal(1:length(st.sw0.params), ((st.par_err)./abs(st.sw0.params))', st.I_best(1), st.I_truest(1), st.I_leader(1));
            plot_generic_frame_data(ax_, x_vec, Y_mat, struct('base_color', base_color, 'I1', I_best, 'I2', I_truest, 'I3', I_leader));
            apply_plot_labels(ax_, 'indices', 'parameter error', 'plot_param_error');
            % set(ax_, 'YScale', 'log');
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
    plot(ax_, x, mean(Y,1), ' -', 'Color', [0 0 0], 'LineWidth', 2);
    plot(ax_, x, Y(plotspecs_.I1, :), ' --', 'Color', [0 0 0], 'LineWidth', 2);
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
