classdef swirl_group
    properties
        % identity swirl_element
        sw0;
        %{a cell object that contains each swirl's position data and dish data%}
        gp_cell;
        %{the number of swirl instances collected by this swirl group%}
        len_gp;
        %{the frame duration of each swirl in this swirl group%}
        frame_vec;
        %{the parameters associated with each swirl in this swirl group%}
        params_mat;
    end
    methods
        function obj = swirl_group(sw_, gp_cell_, len_gp_, frame_vec_, params_mat_)
            obj.sw0 = sw_;
            obj.gp_cell = gp_cell_;
            obj.len_gp = len_gp_;
            obj.frame_vec = frame_vec_;
            obj.params_mat = params_mat_;
        end
        function swi = spawn_swirl_i(obj,i_)
            Frames = obj.frame_vec(i_);
            len_dish = obj.sw0.len_dish;
            len_pos = obj.sw0.len_pos;
            beads = obj.sw0.beads;
            [pos,dish] = swirl_group.cellrow2posdish(obj.gp_cell,i_,beads,len_pos);

            swi = swirl(pos,dish,Frames,len_dish,len_pos,beads);
            swi.params=obj.params_mat(:,i_);
            swi.len_par=size(obj.params_mat,1);
        end
        function frscores = compute_frscores(obj)
            sw0=obj.sw0;
            len_gp = obj.len_gp;
            gp_cell = obj.gp_cell;

            Frames = sw0.Frames;
            beads = sw0.beads;
            len_pos = sw0.len_pos;
            pos_true = sw0.pos;

            frscores = nan(len_gp,1);

            for i = 1:len_gp
                frscores(i) = swirl.compute_frscore(Frames, 1, pos_true-reshape(gp_cell{i,1},[beads,len_pos,Frames]));
            end
        end
        function plot_frame_error(obj, ax_, color_left, color_right, t_vec, pos_err, I_best, I_truest)
            [len_gp, Frames] = size(pos_err);
            mean_err = mean(pos_err);

            box(ax_,'on');
            yyaxis(ax_, 'left');
            set(ax_, 'YScale', 'log');
            yyaxis(ax_, 'right');
            set(ax_, 'YScale', 'log');
            hold(ax_, 'on');
            for i=1:len_gp
                yyaxis(ax_, 'left');
                plot(ax_,t_vec, pos_err(i,:), ' -', 'Color', [color_left, 0.1], 'LineWidth', 1);
                yyaxis(ax_, 'right');
                plot(ax_,t_vec, cumsum(pos_err(i,:)), ' -', 'Color', [color_right, 0.1], 'LineWidth', 1);
            end
            yyaxis(ax_, 'left');
            plot(ax_, t_vec, mean_err, ' -', 'Color', color_left, 'LineWidth', 2);
            plot(ax_, t_vec, pos_err(I_best, :), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec, pos_err(I_truest, :), ' :', 'Color', [0 0 0], 'LineWidth', 2);
            yyaxis(ax_, 'right');
            plot(ax_, t_vec, cumsum(mean_err), ' -', 'Color', color_right, 'LineWidth', 2);
            plot(ax_, t_vec, cumsum(pos_err(I_best, :)), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec, cumsum(pos_err(I_truest, :)), ' :', 'Color', [0 0 0], 'LineWidth', 2);

            xlabel(ax_, 'time', 'Interpreter', 'Latex', 'Fontsize', 14)
            yyaxis(ax_, 'left');
            ylabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            yyaxis(ax_, 'right');
            ylabel(ax_, 'cumulative position error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_param_error(obj, ax_, base_color, par_err_, par_true, I_best, I_truest)
            par_err = par_err_./abs(par_true);
            [len_par, len_gp] = size(par_err);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_, 1:len_par, par_err(:,i), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(ax_, 1:len_par, par_err(:,I_best), ' -', 'Color', [0 0 0], 'LineWidth', 1);
            plot(ax_, 1:len_par, par_err(:,I_truest), ' :', 'Color', [0 0 0], 'LineWidth', 1);
            xlabel(ax_, 'parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'parameter error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_err_vs_accerr(obj, ax_, base_color, pos_err, I_best, I_truest)
            [len_gp, Frames] = size(pos_err);

            box(ax_,'on');
            % set(ax_, 'YScale', 'log');
            % set(ax_, 'XScale', 'log');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_,cumsum(pos_err(i,:)), pos_err(i,:), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(ax_, cumsum(pos_err(I_best, :)), pos_err(I_best, :), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, cumsum(pos_err(I_truest, :)), pos_err(I_truest, :), ' :', 'Color', [0 0 0], 'LineWidth', 2);

            xlabel(ax_, 'cumulative position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_param_covariance(obj, ax_, base_color, par_cov)
            % [len_gp, len_par] = size(par_cov);
            % box(ax_,'on');
            % hold(ax_, 'on');
            % plot(ax_, 1:len_par, par_cov, ' -', 'Color', base_color, 'LineWidth', 1);
            % xlabel(ax_, 'parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
            % ylabel(ax_, 'cov', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_derr_vs_err(obj, ax_, base_color, del_t, pos_err, I_best, I_truest)
            [len_gp, Frames] = size(pos_err);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_,pos_err(i,2:Frames), (pos_err(i,2:Frames)-pos_err(i,1:(Frames-1)))/del_t, ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            % plot(ax_, pos_err(I_best,2:end), (pos_err(I_best,2:end)-pos_err(I_best,1:end-1))/(del_t), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            % plot(ax_, cumsum(pos_err(I_truest,2:end), (pos_err(I_truest,2:end)-pos_err(I_truest,1:end-1))/(del_t), ' :', 'Color', [0 0 0], 'LineWidth', 2);

            xlabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'Dt position error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_pos_vs_param_error(obj, ax_, base_color, par_err, par_true, pos_err_acc, I_best, I_truest)
            [len_par,len_gp] = size(par_err);
            par_err_absum = sum(abs(par_err)./abs(par_true),1);
            box(ax_,'on');
            hold(ax_, 'on');
            scatter(ax_, par_err_absum, pos_err_acc, ' o', 'MarkerFaceColor', base_color, 'SizeData', 50);
            scatter(ax_, par_err_absum(I_best), pos_err_acc(I_best), 'p', 'MarkerFaceColor', [0 0 0], 'SizeData', 80);
            scatter(ax_, par_err_absum(I_truest), pos_err_acc(I_truest), 'h', 'MarkerFaceColor', [0 0 0], 'SizeData', 80);
            xlabel('param error', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel('position error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
    end
    methods (Static)
        function [pos, dish] = cellrow2posdish(gp_cell_, i_, beads, len_pos)
            pos=reshape(gp_cell_{i_,1},[beads,len_pos,size(gp_cell_{i_,1},2)]);
            dish=gp_cell_{i_,2};
        end
        function make_movie_comp(AYfig_in,movie_data,movie_specs)
            Frames = movie_specs.Frames;
            AYfig_in.init_movie(Frames);

            movie_fig = AYfig_in.fig;
            movie_ax = AYfig_in.ax;
            movie_gen = AYfig_in.movie_gen;

            dish = movie_specs.dish;

            pos = cat(1,movie_data{:,1});
            colors = cat(1,movie_data{:,2});

            walld = 5.72;
            wallL = (2/sqrt(3))*walld;
            wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];

            wsca = max(abs([min(dish(:, 2))-wallL, max(dish(:, 2))+walld, min(dish(:, 3))-wallL, max(dish(:, 3))+walld]));
            lims = [-wsca, wsca, -wsca, wsca];

            figdims = AYfig_in.get_dims();
            MS = figdims(4)/(4*wsca); %% radius in pixels
            SS = pi*MS*MS; %% area in pixels

            for i=1:Frames
                plot(movie_ax, wallv(:, 1)+dish(i, 2), wallv(:, 2)+dish(i, 3), 'k -')
                hold(AYfig_in.ax, 'on');
                dots = scatter(movie_ax, pos(:, 1, i), pos(:, 2, i), 'o', 'filled', 'CData', colors, 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);
                txtbx = annotation(movie_fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(i-1), 'LineStyle', 'none', 'FontSize', 16);
                hold(movie_ax, 'off');
                axis(movie_ax, lims);
                drawnow
                movie_gen(i) = getframe(movie_fig);
                delete(txtbx);
            end
            AYfig_in.movie_gen = movie_gen;
        end
    end
end
