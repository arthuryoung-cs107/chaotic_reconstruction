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
        function plot_alpha_bead(AYfig_, base_color1, base_color2, base_color3, te, t_vec)

            plot_generic_frame_data(AYfig_.ax_tile(1), t_vec, te.alpha_bead(:,:,1), struct('base_color', base_color1));
            plot_generic_frame_data(AYfig_.ax_tile(2), t_vec, te.alpha_bead(:,:,2), struct('base_color', base_color2));
            plot_generic_frame_data(AYfig_.ax_tile(3), t_vec, te.alpha_bead(:,:,3), struct('base_color', base_color3));

            apply_plot_labels(AYfig_.ax_tile(1), '$$t_f$$', '$$\alpha_1 = \frac{d \mathcal{R}_1^*(t_f)}{d t_f^*}$$', '');
            apply_plot_labels(AYfig_.ax_tile(2), '$$t_f$$', '$$\alpha_2 = \frac{d \mathcal{R}_2^*(t_f)}{d t_f^*}$$', '');
            apply_plot_labels(AYfig_.ax_tile(3), '$$t_f$$', '$$\alpha_3 = \frac{d \mathcal{R}_3^*(t_f)}{d t_f^*}$$', '');

            reslims = [min(te.alpha_bead(:)),max(te.alpha_bead(:))];
            ylim(AYfig_.ax_tile(1), reslims);
            ylim(AYfig_.ax_tile(2), reslims);
            ylim(AYfig_.ax_tile(3), reslims);

            set(AYfig_.ax_tile, 'YScale', 'log');

            txtbx_a = annotation(AYfig_.fig, 'textbox', [0.075 0.125 0.1 0.1], 'String', '\textbf{a)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_b = annotation(AYfig_.fig, 'textbox', [0.4 0.125 0.1 0.1], 'String', '\textbf{b)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_c = annotation(AYfig_.fig, 'textbox', [0.7325 0.125 0.1 0.1], 'String', '\textbf{c)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);

        end
        function plot_alpha_truenoise(AYfig_, base_color1t, base_color1n, base_color2t, base_color2n, t_vec, alpha_matt, alpha_matn, INTres_matt, INTres_matn)

            plot_generic_frame_data(AYfig_.ax_tile(1), t_vec, alpha_matt, struct('base_color', base_color1t));
            plot_generic_frame_data(AYfig_.ax_tile(2), t_vec, alpha_matn, struct('base_color', base_color1n));

            plot_generic_frame_data(AYfig_.ax_tile(3), t_vec, INTres_matt, struct('base_color', base_color2t));
            plot_generic_frame_data(AYfig_.ax_tile(4), t_vec, INTres_matn, struct('base_color', base_color2n));

            apply_plot_labels(AYfig_.ax_tile(1), '$$t_f$$', '$$\alpha_t = \frac{d \mathcal{R}_t^*(t_f)}{d t_f^*}$$', '');
            apply_plot_labels(AYfig_.ax_tile(2), '$$t_f$$', '$$\alpha_n = \frac{d \mathcal{R}_n^*(t_f)}{d t_f^*}$$', '');

            apply_plot_labels(AYfig_.ax_tile(3), '$$t_f$$', '$$\mathcal{R}_t(t_f) = \int_{t_0}^{t_f} (\hat{\mathbf{p}}(t)-\mathbf{p}(t))^T(\hat{\mathbf{p}}(t)-\mathbf{p}(t)) dt$$', '');
            apply_plot_labels(AYfig_.ax_tile(4), '$$t_f$$', '$$\mathcal{R}_n(t_f) = \int_{t_0}^{t_f} (\tilde{\mathbf{p}}(t)-\mathbf{p}(t))^T(\tilde{\mathbf{p}}(t)-\mathbf{p}(t)) dt$$', '');

            AYfig_.ax_tile(2).YLim = AYfig_.ax_tile(1).YLim;

            reslims = [min(min(INTres_matt(:, 2:end))),max(INTres_matn(:))];
            ylim(AYfig_.ax_tile(3), reslims);
            ylim(AYfig_.ax_tile(4), reslims);

            % AYfig_.ax_tile(4).YLim = AYfig_.ax_tile(3).YLim;
            set(AYfig_.ax_tile, 'YScale', 'log');

            txtbx_a = annotation(AYfig_.fig, 'textbox', [0.41 0.525 0.1 0.1], 'String', '\textbf{a)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_b = annotation(AYfig_.fig, 'textbox', [0.9 0.525 0.1 0.1], 'String', '\textbf{b)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_c = annotation(AYfig_.fig, 'textbox', [0.41 0.025 0.1 0.1], 'String', '\textbf{c)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_d = annotation(AYfig_.fig, 'textbox', [0.9 0.025 0.1 0.1], 'String', '\textbf{d)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);

        end
        function plot_stat_residuals(AYfig_, base_color1, base_color2, st)
            axs_ = AYfig_.ax_tile;
            [f_vec, pos_res, cl] = deal(1:st.sw0.Frames, st.pos_res, st.sw0.params(13));

            plot_generic_frame_data(axs_(1), f_vec, cl*cl*pos_res, struct('base_color', base_color1));
            plot(axs_(1), f_vec, mean(cl*cl*pos_res,1), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(axs_(1), f_vec, cl*cl*pos_res(100, :), ' --', 'Color', [0 0 0], 'LineWidth', 2);
            apply_plot_labels(axs_(1), '$$f$$', '$$(\hat{\mathbf{p}}_f - \mathbf{p}_f)^T(\hat{\mathbf{p}}_f - \mathbf{p}_f)$$', '');
            set(axs_(1), 'YScale', 'log');

            plot_generic_frame_data(axs_(2), f_vec, cumsum(cl*cl*pos_res, 2), struct('base_color', base_color2));
            plot(axs_(2), f_vec, mean(cumsum(cl*cl*pos_res,2),1), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(axs_(2), f_vec, cumsum(cl*cl*pos_res(100, :)), ' --', 'Color', [0 0 0], 'LineWidth', 2);
            apply_plot_labels(axs_(2), '$$f$$', '$$\sum_{i=1}^f (\hat{\mathbf{p}}_i - \mathbf{p}_i)^T(\hat{\mathbf{p}}_i - \mathbf{p}_i)$$', '');
            set(axs_(2), 'YScale', 'log');

            txtbx_a = annotation(AYfig_.fig, 'textbox', [0.41 0.15 0.1 0.1], 'String', '\textbf{a)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_b = annotation(AYfig_.fig, 'textbox', [0.9 0.15 0.1 0.1], 'String', '\textbf{b)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
        end
        function set_collision_fig(AYfig_)
            figure(AYfig_.fig.Number)
            AYfig_.ax_sub = gobjects(6, 1);
            H=AYfig_.fig.Position(4);
            W=AYfig_.fig.Position(3);
            h2w=H/W;
            w2h=W/H;

            hs = 0.025;
            hhs = hs/2;
            h=1/3-hhs;
            w=h2w*h;
            ws =(1.0-2*w)/3;

            Ths = 2*hs;
            Ty = 2*h+hs+Ths;

            Tws = 0.70*ws;
            Tw = ws+w-Tws;
            Th=1-Ty-hhs;

            figure(AYfig_.fig.Number);

            AYfig_.ax_sub(1) = subplot('Position', [Tws,Ty,Tw, Th]);
            AYfig_.ax_sub(2) = subplot('Position', [w+2*ws,Ty,Tw, Th]);
            AYfig_.ax_sub(3) = subplot('Position', [ws,h+hs,w, h]);
            AYfig_.ax_sub(4) = subplot('Position', [w+2*ws,h+hs,w, h]);
            AYfig_.ax_sub(5) = subplot('Position', [ws,hhs,w, h]);
            AYfig_.ax_sub(6) = subplot('Position', [w+2*ws,hhs,w, h]);

            txtbx_a = annotation(AYfig_.fig, 'textbox', [0.4 0.65 0.1 0.1], 'String', '\textbf{a)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_b = annotation(AYfig_.fig, 'textbox', [0.875 0.65 0.1 0.1], 'String', '\textbf{b)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_c = annotation(AYfig_.fig, 'textbox', [ws 0.55 0.1 0.1], 'String', '\textbf{c)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_d = annotation(AYfig_.fig, 'textbox', [w+2*ws 0.55 0.1 0.1], 'String', '\textbf{d)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_e = annotation(AYfig_.fig, 'textbox', [ws 0.225 0.1 0.1], 'String', '\textbf{e)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_f = annotation(AYfig_.fig, 'textbox', [w+2*ws 0.225 0.1 0.1], 'String', '\textbf{f)}', 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);

        end
        function plot_collision_frames(AYfig_, st, base_color1, base_color2, base_color3, base_color4)
            axs_ = AYfig_.ax_sub;
            [f_vec, pos_res, cl] = deal(1:st.sw0.Frames, st.pos_res, st.sw0.params(13));

            plot_generic_frame_data(axs_(1), f_vec, cl*cl*pos_res, struct('base_color', base_color1));
            plot(axs_(1), f_vec, mean(cl*cl*pos_res,1), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(axs_(1), f_vec, cl*cl*pos_res(100, :), ' --', 'Color', [0 0 0], 'LineWidth', 2);
            apply_plot_labels(axs_(1), '$$f$$', '$$(\hat{\mathbf{p}}_f - \mathbf{p}_f)^T(\hat{\mathbf{p}}_f - \mathbf{p}_f)$$', '');
            set(axs_(1), 'YScale', 'log');

            plot_generic_frame_data(axs_(2), f_vec, cumsum(cl*cl*pos_res, 2), struct('base_color', base_color2));
            plot(axs_(2), f_vec, mean(cumsum(cl*cl*pos_res,2),1), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(axs_(2), f_vec, cumsum(cl*cl*pos_res(100, :)), ' --', 'Color', [0 0 0], 'LineWidth', 2);
            apply_plot_labels(axs_(2), '$$f$$', '$$\sum_{i=1}^f (\hat{\mathbf{p}}_i - \mathbf{p}_i)^T(\hat{\mathbf{p}}_i - \mathbf{p}_i)$$', '');
            set(axs_(2), 'YScale', 'log');

            sw_ = st.sw0;
            [rk,ik] = mink(sum(st.pos_res(:,9:150),2), 1); %% interesting case, where we have very low initial residual
            sb_ = st.spawn_swirl_i(ik);
            dish = sw_.dish;
            Frames = sw_.Frames;
            % fcontact = [1, sw_.contact_f(2), sw_.contact_f(4), sw_.contact_f(5)];
            fcontact = [1, sw_.contact_f(1), sw_.contact_f(2), sw_.contact_f(5)];

            walld = 5.72;
            wallL = (2/sqrt(3))*walld;
            wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];

            wsca = wallL;
            lims = [-wsca, wsca, -wsca, wsca];

            figscale = AYfig_.fig.Position(4)*AYfig_.ax_sub(3).Position(4);

            MS = figscale/(4*wsca); %% radius in pixels
            SS = pi*MS*MS; %% area in pixels

            for i = 1:4
                axi = AYfig_.ax_sub(i+2);
                hold(axi, 'on');
                plot(axi, wallv(:, 1), wallv(:, 2), 'k -')
                scatter(axi, sw_.pos(:, 1, fcontact(i))-dish(fcontact(i),2), sw_.pos(:, 2, fcontact(i))-dish(fcontact(i),3), 'o', 'filled', 'CData', base_color3, 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);
                scatter(axi, sb_.pos(:, 1, fcontact(i))-dish(fcontact(i),2), sb_.pos(:, 2, fcontact(i))-dish(fcontact(i),3), 'o', 'filled', 'CData', base_color4, 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);

                set(axi,'xcolor','none')
                set(axi,'ycolor','none')
                axis(axi, lims);
            end

            txtbx_f1 = annotation(AYfig_.fig, 'textbox', [0.4 0.55 0.1 0.1], 'String', ['$$f=' num2str(fcontact(1)) '$$'], 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_f2 = annotation(AYfig_.fig, 'textbox', [0.825 0.55 0.1 0.1], 'String', ['$$f=' num2str(fcontact(2)) '$$'], 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_f3 = annotation(AYfig_.fig, 'textbox', [0.4 0.225 0.1 0.1], 'String', ['$$f=' num2str(fcontact(3)) '$$'], 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
            txtbx_f4 = annotation(AYfig_.fig, 'textbox', [0.825 0.225 0.1 0.1], 'String', ['$$f=' num2str(fcontact(4)) '$$'], 'Interpreter', 'Latex', 'LineStyle', 'none', 'FontSize', 16);
        end
        function plot_convergence(ax_, base_color, re_, ind_, par_true_)
            [par_lineage, par_mat] = re_.get_lineage_mat(re_.gen(ind_(end)).rec(1));

            x = 1:7;
            Y=((par_true_(x)-par_mat(x,:))./abs(par_true_(x)))';

            [len_gp, Frames] = size(Y);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_,x, Y(i,:), ' - o', 'Color', [base_color, 0.1+(i-1)*(0.9)/(len_gp-0.999)], 'LineWidth', 2);
            end
            plot(ax_, x, Y(len_gp,:), ' - p', 'Color', [0 0 0], 'LineWidth', 2);
            xlim(ax_, [min(x), max(x)]);
            ylim(ax_, [min(min(Y)), max(max(Y))]);

            names = {'$$K$$';'$$\gamma_b$$';'$$\gamma_f$$';'$$\gamma_w$$';'$$\mu_b$$';'$$\mu_f$$';'$$\mu_w$$'};

            ylabel(ax_, '$$\frac{\hat{x}_i - x_i}{|\hat{x}_i|}$$', 'Interpreter', 'Latex', 'Fontsize', 14)
            % xlabel(ax_, 'parameters', 'Interpreter', 'Latex', 'Fontsize', 14)
            set(ax_,'xtick',x,'xticklabel',names, 'TickLabelInterpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_K_vs_gamma(ax_, base_color, re_, ind_, par_true_)
            [par_lineage, par_mat] = re_.get_lineage_mat(re_.gen(ind_(end)).rec(1));

            x = 1:7;

            Y=((par_true_(x)-par_mat(x,:))./abs(par_true_(x)))';

            [len_gp, Frames] = size(Y);
            box(ax_,'on');
            hold(ax_, 'on');

            for i=1:len_gp
                scatter(ax_,par_mat(1,i), par_mat(2,i), ' o', 'filled', 'CData', base_color, 'LineWidth', 0.5, 'SizeData', 50, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.1+(i-1)*(0.9)/(len_gp-0.999));
            end
            % plot(ax_, x, Y(len_gp,:), ' - p', 'Color', [0 0 0], 'LineWidth', 2);
            % xlim(ax_, [min(x), max(x)]);
            % ylim(ax_, [min(min(Y)), max(max(Y))]);

            % names = {'$$K$$';'$$\gamma_b$$';'$$\gamma_f$$';'$$\gamma_w$$';'$$\mu_b$$';'$$\mu_f$$';'$$\mu_w$$'};

            ylabel(ax_, '$$K$$', 'Interpreter', 'Latex', 'Fontsize', 14)
            xlabel(ax_, '$$\gamma_b$$', 'Interpreter', 'Latex', 'Fontsize', 14)
            % set(ax_,'xtick',x,'xticklabel',names, 'TickLabelInterpreter', 'Latex', 'Fontsize', 14)
        end
    end
end

function plot_generic_frame_data(ax_, x, Y, plotspecs_)
    [len_gp, Frames] = size(Y);
    box(ax_,'on');
    hold(ax_, 'on');
    for i=1:len_gp
        plot(ax_,x, Y(i,:), ' -', 'Color', [plotspecs_.base_color, 0.05], 'LineWidth', 1);
    end
    xlim(ax_, [min(x), max(x)]);
    ylim(ax_, [min(min(Y)), max(max(Y))]);
end
function set_movie_dims(AYfig_)
    axdims = get(AYfig_.ax, 'Position');
    axdims(3:4) = min(axdims(3:4));
    set(AYfig_.ax, 'Position', axdims);
end
function apply_plot_labels(ax_, xname, yname, titlename)
    xlabel(ax_, xname, 'Interpreter', 'Latex', 'Fontsize', 14)
    ylabel(ax_, yname, 'Interpreter', 'Latex', 'Fontsize', 14)
    title(ax_, ['\textbf{', strrep(titlename, '_', ' ') , '}'], 'Interpreter', 'Latex', 'Fontsize', 14)
end
