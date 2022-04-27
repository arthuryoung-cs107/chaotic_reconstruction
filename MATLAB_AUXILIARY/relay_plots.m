classdef relay_plots
    properties

    end
    methods
        function obj = relay_plots()

        end
    end
    methods (Static)

        %% error bar plots
        function plot_3bead_event_stats(axs_, color1, color2, color3, re_)

            histogram(axs_(1), re_.event0.bead_events_full(1,:), 'FaceColor', color1)
            histogram(axs_(2), re_.event0.bead_events_full(2,:), 'FaceColor', color1)
            histogram(axs_(3), re_.event0.bead_events_full(3,:), 'FaceColor', color1)
            apply_plot_labels(axs_(1), 'event frame', 'count', 'bead_1_collision_frame');
            apply_plot_labels(axs_(2), 'event frame', 'count', 'bead_2_collision_frame');
            apply_plot_labels(axs_(3), 'event frame', 'count', 'bead_3_collision_frame');

            histogram(axs_(4), re_.event0.bead_res(1,:), 'FaceColor', color2)
            histogram(axs_(5), re_.event0.bead_res(2,:), 'FaceColor', color2)
            histogram(axs_(6), re_.event0.bead_res(3,:), 'FaceColor', color2)
            apply_plot_labels(axs_(4), 'final residual', 'count', 'bead_1_final_residual');
            apply_plot_labels(axs_(5), 'final residual', 'count', 'bead_2_final_residual');
            apply_plot_labels(axs_(6), 'final residual', 'count', 'bead_3_final_residual');

            histogram(axs_(7), re_.event0.bead_alpha(1,:), 'FaceColor', color3)
            histogram(axs_(8), re_.event0.bead_alpha(2,:), 'FaceColor', color3)
            histogram(axs_(9), re_.event0.bead_alpha(3,:), 'FaceColor', color3)
            apply_plot_labels(axs_(7), 'final alpha', 'count', 'bead_1_final_alpha');
            apply_plot_labels(axs_(8), 'final alpha', 'count', 'bead_2_final_alpha');
            apply_plot_labels(axs_(9), 'final alpha', 'count', 'bead_3_final_alpha');
        end
        function plot_param_error( ax_, base_color, par_mat_, par_true_)
            plot_generic_frame_data(ax_, 1:length(par_true_), ((par_true_-par_mat_)./abs(par_true_))', struct('base_color', base_color));
            apply_plot_labels(ax_, 'indices', 'param err', 'plot_param_error');
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
    xlim(ax_, [min(x), max(x)]);
    ylim(ax_, [min(min(Y)), max(max(Y))]);
end
function plot_significant_indices(ax_, x, Y, plotspecs_)
    plot(ax_, x, Y(plotspecs_.I1, :), ' --', 'Color', [0 0 0], 'LineWidth', 2);
    plot(ax_, x, Y(plotspecs_.I2, :), ' :', 'Color', [0 0 0], 'LineWidth', 2);
    plot(ax_, x, Y(plotspecs_.I3, :), ' -.', 'Color', [0 0 0], 'LineWidth', 2);
end
function apply_plot_labels(ax_, xname, yname, titlename)
    xlabel(ax_, xname, 'Interpreter', 'Latex', 'Fontsize', 14)
    ylabel(ax_, yname, 'Interpreter', 'Latex', 'Fontsize', 14)
    title(ax_, ['\textbf{', strrep(titlename, '_', ' ') , '}'], 'Interpreter', 'Latex', 'Fontsize', 14)
end
