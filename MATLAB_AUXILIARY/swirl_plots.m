classdef swirl_plots
    properties

    end
    methods
        function obj = swirl_plots()

        end
    end
    methods (Static)
        function plot_frame_error( ax_, color_left, color_right, t_vec, pos_err, I_best, I_truest, I_leader)
            [len_gp, Frames] = size(pos_err);
            mean_err = mean(pos_err);

            box(ax_,'on');
            yyaxis(ax_, 'left');
            % set(ax_, 'YScale', 'log');
            yyaxis(ax_, 'right');
            % set(ax_, 'YScale', 'log');
            hold(ax_, 'on');
            for i=1:len_gp
                yyaxis(ax_, 'left');
                plot(ax_,t_vec, pos_err(i,:), ' -', 'Color', [color_left, 0.1], 'LineWidth', 1);
                yyaxis(ax_, 'right');
                plot(ax_,t_vec, cumsum(pos_err(i,:)), ' -', 'Color', [color_right, 0.1], 'LineWidth', 1);
            end
            yyaxis(ax_, 'left');
            plot(ax_, t_vec, pos_err(I_best, :), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec, pos_err(I_truest, :), ' :', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec, pos_err(I_leader, :), ' -.', 'Color', [0 0 0], 'LineWidth', 2);
            yyaxis(ax_, 'right');
            plot(ax_, t_vec, cumsum(pos_err(I_best, :)), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec, cumsum(pos_err(I_truest, :)), ' :', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec, cumsum(pos_err(I_leader, :)), ' -.', 'Color', [0 0 0], 'LineWidth', 2);

            xlim(ax_, [min(t_vec), max(t_vec)]);

            xlabel(ax_, 'time', 'Interpreter', 'Latex', 'Fontsize', 14)
            yyaxis(ax_, 'left');
            ylim(ax_, [0, max(max(pos_err))]);
            ylabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            yyaxis(ax_, 'right');
            ylim(ax_, [0, max(sum(pos_err,2))]);
            ylabel(ax_, 'accumulated position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot frame error}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_param_error_leaders( ax_, base_color, par_err_, par_true, I_best, I_truest, I_leaders, nlead)
            par_err = par_err_./abs(par_true);
            [len_par, len_gp] = size(par_err);
            I_leader = I_leaders(1);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=reshape(I_leaders(1:nlead), 1,[])
                plot(ax_, 1:len_par, par_err(:,i), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(ax_, 1:len_par, par_err(:,I_best), ' -', 'Color', [0 0 0], 'LineWidth', 1);
            plot(ax_, 1:len_par, par_err(:,I_truest), ' :', 'Color', [0 0 0], 'LineWidth', 1);
            plot(ax_, 1:len_par, par_err(:,I_leader), ' -.', 'Color', [0 0 0], 'LineWidth', 1);
            xlabel(ax_, 'parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'parameter error', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot param error leaders}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_param_error_best( ax_, base_color, par_err_, par_true, I_bestest, I_truest, I_leader, nbest)
            par_err = par_err_./abs(par_true);
            [len_par, len_gp] = size(par_err);
            I_best = I_bestest(1);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=reshape(I_bestest(1:nbest), 1,[])
                plot(ax_, 1:len_par, par_err(:,i), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(ax_, 1:len_par, par_err(:,I_best), ' -', 'Color', [0 0 0], 'LineWidth', 1);
            plot(ax_, 1:len_par, par_err(:,I_truest), ' :', 'Color', [0 0 0], 'LineWidth', 1);
            plot(ax_, 1:len_par, par_err(:,I_leader), ' -.', 'Color', [0 0 0], 'LineWidth', 1);
            xlabel(ax_, 'parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'parameter error', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot param error best}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_err_vs_accerr( ax_, base_color, pos_err, I_best, I_truest, I_leader)
            [len_gp, Frames] = size(pos_err);

            box(ax_,'on');
            set(ax_, 'YScale', 'log');
            set(ax_, 'XScale', 'log');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_,cumsum(pos_err(i,:)), pos_err(i,:), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(ax_, cumsum(pos_err(I_best, :)), pos_err(I_best, :), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, cumsum(pos_err(I_truest, :)), pos_err(I_truest, :), ' :', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, cumsum(pos_err(I_leader, :)), pos_err(I_leader, :), ' -.', 'Color', [0 0 0], 'LineWidth', 2);
            xlabel(ax_, 'cumulative position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot err vs accerr}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_param_covariance( ax_, base_color, par_cov)

        end
        function plot_derr_vs_err( ax_, base_color, dish, pos_err, I_best, I_truest, I_leader)
            [len_gp, Frames] = size(pos_err);
            del_t = dish(2,1);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_,pos_err(i,2:Frames), (pos_err(i,2:Frames)-pos_err(i,1:(Frames-1)))/del_t, ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(ax_,pos_err(I_best,2:Frames), (pos_err(I_best,2:Frames)-pos_err(I_best,1:(Frames-1)))/del_t, ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_,pos_err(I_truest,2:Frames), (pos_err(I_truest,2:Frames)-pos_err(I_truest,1:(Frames-1)))/del_t, ' :', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_,pos_err(I_leader,2:Frames), (pos_err(I_leader,2:Frames)-pos_err(I_leader,1:(Frames-1)))/del_t, ' -.', 'Color', [0 0 0], 'LineWidth', 2);
            xlabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'Dt position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot derr vs err}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_pos_vs_param_error( ax_, base_color, par_err, par_true, pos_err_acc, I_best, I_truest, I_leader)
            [len_par,len_gp] = size(par_err);
            par_err_absum = sum(abs(par_err)./abs(par_true),1);

            box(ax_,'on');
            set(ax_, 'YScale', 'log');
            set(ax_, 'XScale', 'log');
            hold(ax_, 'on');
            scatter(ax_, par_err_absum, pos_err_acc, ' o', 'MarkerFaceColor', base_color, 'SizeData', 50);
            scatter(ax_, par_err_absum(I_best), pos_err_acc(I_best), 'p', 'MarkerFaceColor', [0 0 0], 'SizeData', 80);
            scatter(ax_, par_err_absum(I_truest), pos_err_acc(I_truest), 'h', 'MarkerFaceColor', [0 0 0], 'SizeData', 80);
            scatter(ax_, par_err_absum(I_leader), pos_err_acc(I_leader), 'd', 'MarkerFaceColor', [0 0 0], 'SizeData', 80);
            xlabel('param error', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel('position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot pos vs param error}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_frame_error_kill( ax_, color_left, color_right, t_vec, pos_err, frscores, I_best, I_truest, I_leader)
            [len_gp, Frames] = size(pos_err);
            mean_err = mean(pos_err);
            frmean_floor = floor(mean(frscores));

            box(ax_,'on');
            yyaxis(ax_, 'left');
            % set(ax_, 'YScale', 'log');
            yyaxis(ax_, 'right');
            % set(ax_, 'YScale', 'log');
            hold(ax_, 'on');
            for i=1:len_gp
                yyaxis(ax_, 'left');
                plot(ax_,t_vec(1:frscores(i)), pos_err(i,1:frscores(i)), ' -', 'Color', [color_left, 0.1], 'LineWidth', 1);
                yyaxis(ax_, 'right');
                plot(ax_,t_vec(1:frscores(i)), cumsum(pos_err(i,1:frscores(i))), ' -', 'Color', [color_right, 0.1], 'LineWidth', 1);
            end
            yyaxis(ax_, 'left');
            plot(ax_, t_vec(1:frmean_floor), mean_err(1:frmean_floor), ' -', 'Color', color_left, 'LineWidth', 2);
            plot(ax_, t_vec(1:frscores(I_best)), pos_err(I_best, 1:frscores(I_best)), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec(1:frscores(I_truest)), pos_err(I_truest, 1:frscores(I_truest)), ' :', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec(1:frscores(I_leader)), pos_err(I_leader, 1:frscores(I_leader)), ' -.', 'Color', [0 0 0], 'LineWidth', 2);
            yyaxis(ax_, 'right');
            plot(ax_, t_vec(1:frmean_floor), cumsum(mean_err(1:frmean_floor)), ' -', 'Color', color_right, 'LineWidth', 2);
            plot(ax_, t_vec(1:frscores(I_best)), cumsum(pos_err(I_best, 1:frscores(I_best))), ' -', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec(1:frscores(I_truest)), cumsum(pos_err(I_truest, 1:frscores(I_truest))), ' :', 'Color', [0 0 0], 'LineWidth', 2);
            plot(ax_, t_vec(1:frscores(I_leader)), cumsum(pos_err(I_leader, 1:frscores(I_leader))), ' -.', 'Color', [0 0 0], 'LineWidth', 2);

            xlim(ax_, [min(t_vec), max(t_vec)]);

            xlabel(ax_, 'time', 'Interpreter', 'Latex', 'Fontsize', 14)
            yyaxis(ax_, 'left');
            ylim(ax_, [0, max(max(pos_err))]);
            ylabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            yyaxis(ax_, 'right');
            ylim(ax_, [0, max(sum(pos_err,2))]);
            ylabel(ax_, 'cumulative position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot frame error kill}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_frscore_poserr( ax_, color_left, color_right, pos_err, frscores, I_best, I_truest, I_leader)
            [len_gp, Frames] = size(pos_err);
            nleaders = 100;

            pos_err_acc = sum(pos_err,2);
            pos_err_acc_kill = nan(len_gp, 1);
            pos_err_final_kill = nan(len_gp, 1);
            for i=1:len_gp
                pos_err_acc_kill(i) = sum(pos_err(i, 1:frscores(i)))/frscores(i);
                pos_err_final_kill(i) = pos_err(i, frscores(i))/frscores(i);
            end

            box(ax_,'on');
            hold(ax_, 'on');

            plot(ax_,frscores, pos_err_acc, ' o', 'Color', [color_left, 0.1], 'LineWidth', 1);
            plot(ax_,frscores(I_best), pos_err_acc(I_best), ' *', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerSize', 8);
            plot(ax_,frscores(I_truest), pos_err_acc(I_truest), ' x', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerSize', 8);
            plot(ax_,frscores(I_leader), pos_err_acc(I_leader), ' +', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerSize', 8);

            ylabel(ax_, 'cumulative position error', 'Interpreter', 'Latex', 'Fontsize', 14)
            xlabel(ax_, 'frame score', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot frscore poserr}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_acc_weights( ax_, color_left, color_right, pos_err, frscores, I_best, I_truest, I_leader)
            [len_gp, Frames] = size(pos_err);
            % nleaders = 100;
            % nleaders = 500;
            nleaders = 250;

            [frleaders, I_leaders] = maxk(frscores, nleaders);

            pos_err_acc = sum(pos_err,2);
            pos_err_acc_f1 = nan(len_gp, 1);
            pos_err_acc_f2 = nan(len_gp, 1);
            for i=1:len_gp
                pos_err_acc_f1(i) = sum(pos_err(i, 1:frscores(i)))/frscores(i)/(frscores(i)-min(frleaders)+1);
                pos_err_acc_f2(i) = sum(pos_err(i, 1:frscores(i)))/(frscores(i)+frscores(i)-min(frleaders));
            end

            box(ax_,'on');
            hold(ax_, 'on');

            plot(ax_,frscores(I_leaders), pos_err_acc_f1(I_leaders), ' o', 'Color', [color_left, 0.1], 'LineWidth', 1);

            ylabel(ax_, 'reweighted 1', 'Interpreter', 'Latex', 'Fontsize', 14)
            xlabel(ax_, 'frame score', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot acc weights}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_acc_weights_hist( ax_, color_left, color_right, pos_err, frscores, I_best, I_truest, I_leader)

            [len_gp, Frames] = size(pos_err);
            % nleaders = 100;
            % nleaders = 500;
            nleaders = 250;

            [frleaders, I_leaders] = maxk(frscores, nleaders);

            pos_err_acc = sum(pos_err,2);
            pos_err_acc_f1 = nan(len_gp, 1);
            pos_err_acc_f2 = nan(len_gp, 1);
            for i=1:len_gp
                pos_err_acc_f1(i) = sum(pos_err(i, 1:frscores(i)))/frscores(i)/(frscores(i)-min(frleaders)+1);
                pos_err_acc_f2(i) = sum(pos_err(i, 1:frscores(i)))/(frscores(i)+frscores(i)-min(frleaders));
            end

            lambda = 1/(mean(pos_err_acc_f1(I_leaders)));

            box(ax_,'on');
            hold(ax_, 'on');

            yyaxis(ax_, 'left');
            ylabel(ax_, 'frequency', 'Interpreter', 'Latex', 'Fontsize', 14)
            histogram(ax_, pos_err_acc_f1(I_leaders));

            yyaxis(ax_, 'right');
            fplot(@(z) lambda*exp(-lambda.*z), [0, max(pos_err_acc_f1(I_leaders))], 'Color', color_right);

            xlabel(ax_, 'z', 'Interpreter', 'Latex', 'Fontsize', 14)
            title(ax_, '\textbf{plot acc weights hist}', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
    end
end
