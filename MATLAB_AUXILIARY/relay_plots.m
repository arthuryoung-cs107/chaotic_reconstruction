classdef relay_plots
    properties

    end
    methods
        function obj = error_plots()

        end
    end
    methods (Static)

        %% error bar plots
        function plot_3bead_event_stats(axs_, color1, color2, re_)

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
        end
    end
end

function apply_plot_labels(ax_, xname, yname, titlename)
    xlabel(ax_, xname, 'Interpreter', 'Latex', 'Fontsize', 14)
    ylabel(ax_, yname, 'Interpreter', 'Latex', 'Fontsize', 14)
    title(ax_, ['\textbf{', strrep(titlename, '_', ' ') , '}'], 'Interpreter', 'Latex', 'Fontsize', 14)
end
