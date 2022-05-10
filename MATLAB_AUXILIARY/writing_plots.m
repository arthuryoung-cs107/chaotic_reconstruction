classdef writing_plots
    properties

    end
    methods
        function obj = writing_plots()
        end
    end
    methods (Static)
        function plot_1bead_swirl(sw_,AYfig_,color_path_, color_col_)
            dish = sw_.dish;
            Frames = sw_.Frames;

            set_movie_dims(AYfig_);

            walld = 5.72;
            wallL = (2/sqrt(3))*walld;
            wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];

            wsca = wallL;
            lims = [-wsca, wsca, -wsca, wsca];

            figdims = AYfig_.get_dims();
            MS = figdims(4)/(4*wsca); %% radius in pixels
            SS = pi*MS*MS; %% area in pixels

            pos = reshape((permute(sw_.pos(:,1:2,:), [3 2 1])), [Frames, 2])-dish(:,2:3);
            contact_pos = pos(sw_.contact_f, :);

            alpha = 0.4;

            hold(AYfig_.ax, 'on');

            plot(AYfig_.ax, wallv(:, 1), wallv(:, 2), 'k -')
            scatter(AYfig_.ax, pos(1, 1), pos(1, 2), 'o', 'filled', 'CData', color_path_, 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);
            for i = 1:Frames-1
                plot(AYfig_.ax, pos(i:i+1, 1), pos(i:i+1, 2), '-', 'Color', [color_path_, alpha], 'LineWidth', 3);
            end
            plot(AYfig_.ax, contact_pos(:, 1), contact_pos(:, 2), 'x', 'Color', [color_col_, alpha], 'LineWidth', 1, 'MarkerSize', MS);
            set(AYfig_.ax,'xcolor','none')
            set(AYfig_.ax,'ycolor','none')
            % scatter(AYfig_.ax, posIC(i, 1), posIC(i, 2), 'x', 'filled', 'CData', color_col_, 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

            axis(AYfig_.ax, lims);
        end
        function plot_ICs(sw_,AYfig_,color_path_)
            dish = sw_.dish;
            Frames = sw_.Frames;

            set_movie_dims(AYfig_);

            walld = 5.72;
            wallL = (2/sqrt(3))*walld;
            wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];

            wsca = wallL;
            lims = [-wsca, wsca, -wsca, wsca];

            figdims = AYfig_.get_dims();
            MS = figdims(4)/(4*wsca); %% radius in pixels
            SS = pi*MS*MS; %% area in pixels

            posIC = reshape(sw_.pos(:,1:2,1),[],2)-dish(1,2:3);

            hold(AYfig_.ax, 'on');

            plot(AYfig_.ax, wallv(:, 1), wallv(:, 2), 'k -')
            for i = 1:sw_.beads
                scatter(AYfig_.ax, posIC(i, 1), posIC(i, 2), 'o', 'filled', 'CData', color_path_, 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);
                pause
            end
            axis(AYfig_.ax, lims);
        end
        function plot_param_error( ax_, base_color, par_mat, par_true_)
            x = 1:7;
            Y=((par_true_(x)-par_mat(x,:))./abs(par_true_(x)))';

            [len_gp, Frames] = size(Y);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_,x, Y(i,:), ' - o', 'Color', [base_color, 0.1], 'LineWidth', 2);
            end
            xlim(ax_, [min(x), max(x)]);
            ylim(ax_, [min(min(Y)), max(max(Y))]);

            names = {'$$K$$';'$$\gamma_b$$';'$$\gamma_f$$';'$$\gamma_w$$';'$$\mu_b$$';'$$\mu_f$$';'$$\mu_w$$'};

            ylabel(ax_, '$$\frac{\hat{x}_i - x_i}{|\hat{x}_i|}$$', 'Interpreter', 'Latex', 'Fontsize', 14)
            % xlabel(ax_, 'parameters', 'Interpreter', 'Latex', 'Fontsize', 14)
            set(ax_,'xtick',x,'xticklabel',names, 'TickLabelInterpreter', 'Latex', 'Fontsize', 14)
        end
    end
end

function set_movie_dims(AYfig_)
    axdims = get(AYfig_.ax, 'Position');
    axdims(3:4) = min(axdims(3:4));
    set(AYfig_.ax, 'Position', axdims);
end
