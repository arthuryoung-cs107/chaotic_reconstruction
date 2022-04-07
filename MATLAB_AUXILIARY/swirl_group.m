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
        function [par_err, pos_err, pos_err_acc] = compute_error(obj)
            [len_gp, Frames, len_par, par0, pos0, pars, pos, tensdims] = deal(obj.len_gp, obj.sw0.Frames, obj.sw0.len_par, obj.sw0.params, obj.sw0.pos(:, 1:2, :), obj.params_mat, obj.gp_cell(:,1), size(obj.sw0.pos));

            par_err = nan(len_gp, len_par);
            pos_err = nan(len_gp, Frames-1);
            pos_err_acc = nan(len_gp, 1);
            for i=1:len_gp
                par_err(i,:) = (par0-pars(:,i))./abs(par0);
                pos_err(i,:) = compute_frame_error(pos0,pos{i},tensdims);
                pos_err_acc(i) = sum(pos_err(i,:));
            end
        end
        function [I_best, I_truest, par_cov] = compute_statistics(obj, par_err, pos_err, pos_err_acc)
            len_gp = obj.len_gp;
            [pos_err_best, I_best] = mink(pos_err_acc,len_gp);
            [par_err_truest, I_truest] = mink(sum(abs(par_err)'),len_gp);
            par_cov = sum(pos_err_acc.*par_err)./mean(par_err);
        end
        function plot_frame_error(obj, ax_, base_color, pos_err, I_best, I_truest)
            [len_gp, Frames] = size(pos_err);
            mean_err = mean(pos_err);

            set(ax_, 'YScale', 'log');
            box(ax_,'on');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_,1:Frames, pos_err(i,:), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(ax_, 1:Frames, mean_err, ' -', 'Color', base_color, 'LineWidth', 2);
            plot(ax_, 1:Frames, pos_err(I_best, :), ' -', 'Color', [0 0 0], 'LineWidth', 1);
            plot(ax_, 1:Frames, pos_err(I_truest, :), ' :', 'Color', [0 0 0], 'LineWidth', 1);
            xlabel(ax_, 'Frame', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_param_error(obj, ax_, base_color, par_err, I_best, I_truest)
            [len_gp, len_par] = size(par_err);
            box(ax_,'on');
            hold(ax_, 'on');
            for i=1:len_gp
                plot(ax_, 1:len_par, par_err(i, :), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(ax_, 1:len_par, par_err(I_best,:), ' -', 'Color', [0 0 0], 'LineWidth', 1);
            plot(ax_, 1:len_par, par_err(I_truest,:), ' :', 'Color', [0 0 0], 'LineWidth', 1);
            xlabel(ax_, 'parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
            ylabel(ax_, 'position error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_param_covariance(obj, fig_in, base_color)
            figure(fig_in.Number)
            for i=1:obj.noise_len-1
            plot(1:obj.par_len, obj.par_cov(i, :), ' -', 'Color', [base_color, 0.1], 'LineWidth', 1);
            end
            plot(1:obj.par_len, mean(obj.par_cov), ' -', 'Color', base_color, 'LineWidth', 2);
            plot(1:obj.par_len, obj.par_cov(obj.I_best(1), :), ' -', 'Color', [0 0 0 ], 'LineWidth', 1);
            plot(1:obj.par_len, obj.par_cov(obj.I_truest(1), :), ' :', 'Color', [0 0 0 ], 'LineWidth', 1);
            % xlabel('parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
            % ylabel('param error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function plot_param_pos_error(obj, fig_in, base_color)
            figure(fig_in.Number)
            plot(sum(abs(obj.par_err)'), obj.pos_err_sum, ' o', 'Color', base_color, 'LineWidth', 2);
            plot(sum(abs(obj.par_err(obj.I_best(1), :))), obj.pos_err_sum(obj.I_best(1)), ' p', 'Color', [0 0 0], 'LineWidth', 2);
            plot(sum(abs(obj.par_err(obj.I_truest(1), :))), obj.pos_err_sum(obj.I_truest(1)), ' h', 'Color', [0 0 0], 'LineWidth', 2);
            % xlabel('param error', 'Interpreter', 'Latex', 'Fontsize', 14)
            % ylabel('position error', 'Interpreter', 'Latex', 'Fontsize', 14)
        end
        function make_moviei(obj, AYfig_in, i )
            obj.sw(i).make_movie(AYfig_in);
        end
    end
end

function pos_err_out = compute_frame_error(pos0,posi,dims)
    [beads, len_pos, Frames] = deal(dims(1),dims(2),dims(3));
    pos_err_out = nan(1,Frames-1);
    posi_tens = reshape(posi,dims);
    diff = pos0-posi_tens(:,1:2,:);
    for i=1:(Frames-1)
        pos_err_out(i) = sum(sqrt(sum((diff(:,:,i).*diff(:,:,i))')./sum((pos0(:,:,i).*pos0(:,:,i))')));
    end
end
