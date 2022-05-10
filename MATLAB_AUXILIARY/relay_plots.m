classdef relay_plots
    properties

    end
    methods
        function obj = relay_plots()

        end
    end
    methods (Static)

        function plot_convergence(ax_, base_color, re_, ind_, par_true_)
            [par_lineage, par_mat] = re_.get_lineage_mat(re_.gen(ind_(end)).rec(1));

            x = 1:7;
            Y=((par_true_(x)-par_mat(x,:))./abs(par_true_(x)))';

            [len_gp, Frames] = size(Y);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_,x, Y(i,:), ' - o', 'Color', [base_color, 0.1+(i-1)*(0.9)/(len_gp-1)], 'LineWidth', 2);
            end
            plot(ax_, x, Y(len_gp,:), ' - p', 'Color', [0 0 0], 'LineWidth', 2);
            xlim(ax_, [min(x), max(x)]);
            ylim(ax_, [min(min(Y)), max(max(Y))]);

            names = {'$$K$$';'$$\gamma_b$$';'$$\gamma_f$$';'$$\gamma_w$$';'$$\mu_b$$';'$$\mu_f$$';'$$\mu_w$$'};

            ylabel(ax_, '$$\frac{\hat{x}_i - x_i}{|\hat{x}_i|}$$', 'Interpreter', 'Latex', 'Fontsize', 14)
            % xlabel(ax_, 'parameters', 'Interpreter', 'Latex', 'Fontsize', 14)
            set(ax_,'xtick',x,'xticklabel',names, 'TickLabelInterpreter', 'Latex', 'Fontsize', 14)
        end

        function plot_cell_vs_frames(axs_, base_color, x, Ycell, xname, yname, titlename, xind_)
            plot_struct = struct('base_color', base_color);
            if nargin==8
                for i = 1:length(axs_)
                    plot_generic_frame_data(axs_(i), x(xind_), Ycell{i}(:, xind_), plot_struct);
                    apply_plot_labels(axs_(i), xname, yname, [titlename num2str(i)]);
                    set(axs_(i), 'YScale', 'log');
                end
            else
                for i = 1:length(axs_)
                    plot_generic_frame_data(axs_(i), x, Ycell{i}, plot_struct);
                    apply_plot_labels(axs_(i), xname, yname, [titlename num2str(i)]);
                    set(axs_(i), 'YScale', 'log');
                end
            end
        end
        function plot_posres_vs_frames_statcomp(axs_, base_color, det_, stat_)
            plot_struct = struct('base_color', base_color);
            for i = 1:det_.beads
                plot_generic_frame_data(axs_(i), 1:det_.Frame_end, det_.sim_beadres{i}, plot_struct);
                apply_plot_labels(axs_(i), 'Frames', 'position residual', ['detect_res_vs_frame_bead' num2str(i)]);
                set(axs_(i), 'YScale', 'log');
            end
            for i = 1:det_.beads
                plot_generic_frame_data(axs_(i+det_.beads), 1:det_.Frame_end, stat_.pos_res_bead{i}(:,1:det_.Frame_end), plot_struct);
                apply_plot_labels(axs_(i+det_.beads), 'Frames', 'position residual', ['stat_res_vs_frame_bead' num2str(i)]);
                set(axs_(i+det_.beads), 'YScale', 'log');
            end

        end

        %% error bar plots
        function plot_3bead_event_stats(axs_, colors_, re_)
            event0 = re_.event0;
            beads = event0.specs.beads
            for i = 1:beads
                Frame_bins = max(event0.bead_events_full(i,:))-min(event0.bead_events_full(i,:))+1;
                histogram(axs_((i-1)*beads+1), re_.event0.bead_events_full(i,:), Frame_bins, 'FaceColor', colors_(1,:))
                histogram(axs_((i-1)*beads+2), re_.event0.bead_res(i,:), 'FaceColor', colors_(2,:))
                histogram(axs_((i-1)*beads+3), re_.event0.bead_alpha(i,:), 'FaceColor', colors_(3,:))

                apply_plot_labels(axs_((i-1)*beads+1), 'event frame', 'count', ['bead_' num2str(i) '_collision_frame']);
                apply_plot_labels(axs_((i-1)*beads+2), 'final residual', 'count', ['bead_' num2str(i) '_final_residual']);
                apply_plot_labels(axs_((i-1)*beads+3), 'final alpha', 'count', ['bead_' num2str(i) '_final_alpha']);
            end
        end
        function plot_gen_weights(axs_, base_color, re, indices_)
            axi = 1;
            for i = indices_
                histogram(axs_(axi), re.gen(i).sample_weights, 'FaceColor', base_color)
                apply_plot_labels(axs_(axi), 'weights', 'counts', ['gen_' num2str(i) '_weights']);
                axi=axi+1;
            end
        end
        function plot_gen_netresiduals(axs_, base_color, re, indices_)
            axi = 1;
            for i = indices_
                histogram(axs_(axi), re.gen(i).net_residuals, 'FaceColor', base_color)
                apply_plot_labels(axs_(axi), 'residuals', 'counts', ['gen_' num2str(i) '_net_residuals']);
                axi=axi+1;
            end
        end
        function plot_gen_probabilities(axs_, base_color, re, indices_)
            axi = 1;
            for i = indices_
                histogram(axs_(axi), re.gen(i).pis, 'FaceColor', base_color)
                apply_plot_labels(axs_(axi), '$$\pi$$', 'counts', ['gen_' num2str(i) 'probabilities']);
                axi=axi+1;
            end
        end
        function plot_gen_zeta(axs_, base_color, re, indices_)
            axi = 1;
            for i = indices_
                histogram(axs_(axi), re.gen(i).zetas, 'FaceColor', base_color)
                apply_plot_labels(axs_(axi), '$$\zeta$$', 'counts', ['gen_' num2str(i) 'zeta_values']);
                axi=axi+1;
            end
        end
        function plot_gen_duplication_count(axs_, base_color, re, indices_)
            axi = 1;
            for i = indices_
                lead_bins = max(re.gen(i).lead_dup_count)-min(re.gen(i).lead_dup_count)+1;
                histogram(axs_(axi), re.gen(i).lead_dup_count, lead_bins, 'FaceColor', base_color)
                apply_plot_labels(axs_(axi), 'duplication count', 'counts', ['gen_' num2str(i) '_duplication_count']);
                axi=axi+1;
            end
        end
        function plot_gen_weight_vs_dup(axs_, base_color, re, indices_)
            axi = 1;
            for i = indices_
                box(axs_(axi),'on');
                hold(axs_(axi), 'on');
                plot(axs_(axi),re.gen(i).sample_weights, re.gen(i).lead_dup_count, ' o', 'Color', base_color, 'LineWidth', 2);
                xlim(axs_(axi), [min(re.gen(i).sample_weights), max(re.gen(i).sample_weights)]);
                ylim(axs_(axi), [min(re.gen(i).lead_dup_count), max(re.gen(i).lead_dup_count)]);
                apply_plot_labels(axs_(axi), 'weight', 'duplication count', ['gen_' num2str(i) '_weight_vs_duplication_count']);

                % set(axs_(axi), 'XScale', 'log');
                % set(axs_(axi), 'YScale', 'log');
                axi=axi+1;
            end
        end
        function plot_gen_param_error(axs_, base_color, re, par_true_, indices_)
            axi = 1;
            for i = indices_
                plot_generic_frame_data(axs_(axi), 1:length(par_true_), ((par_true_-re.gen(i).params)./abs(par_true_))', struct('base_color', base_color));

                plot(axs_(axi),1:length(par_true_), (par_true_-re.gen(i).lead_par_w_mean)./abs(par_true_), ' --', 'Color', [0 0 0], 'LineWidth', 4);

                apply_plot_labels(axs_(axi), 'indices', 'param err', ['gen_' num2str(i) '_param_error']);
                axi=axi+1;
            end
        end
        function plot_param_error( ax_, base_color, par_mat_, par_true_)
            plot_generic_frame_data(ax_, 1:length(par_true_), ((par_true_-par_mat_)./abs(par_true_))', struct('base_color', base_color));
            apply_plot_labels(ax_, 'indices', 'param err', 'plot_param_error');
        end
    end
end
function plot_simple(ax_, x, y, plotspecs_)
    box(ax_,'on');
    hold(ax_, 'on');
    plot(ax_,x, y, ' -', 'Color', plotspecs_.base_color, 'LineWidth', 2);
    xlim(ax_, [min(x), max(x)]);
    ylim(ax_, [min(y), max(y)]);
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
