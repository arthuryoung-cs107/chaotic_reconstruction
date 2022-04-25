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
    end
    methods (Static)
        function INT_mat = compute_INT(t_vec_, x_mat_)
            [M,N] = size(x_mat_);
            INT_mat = cumsum([zeros(M,1), 0.5*(x_mat_(:,1:(N-1))+x_mat_(:,2:end)).*reshape(t_vec_(2:end)-t_vec_(1:(N-1)), 1, [])] , 2);
        end
        function D_mat = compute_D(t_vec_, x_mat_)
            [M,N] = size(x_mat_);
            del_t = (t_vec_(2:end)-t_vec_(1:(N-1)));
            D_mat = [(x_mat_(:,2)-x_mat_(:,1))/del_t(1), (x_mat_(:,3:end)-x_mat_(:,1:(N-2)))./reshape(del_t(1:(N-2))+del_t(2:(N-1)) ,1,[]), (x_mat_(:,N)-x_mat_(:,N-1))/del_t(N-1)];
        end
        function frame_stats = compute_residual_time_stats(INTpos_res, pos_res, Dpos_res, frames_)
            if (nargin==4)
                frame_stats = [mean(INTpos_res(:,frames_)); std(INTpos_res(:,frames_)); mean(pos_res(:, frames_)); std(pos_res(:, frames_)); mean(Dpos_res(:,frames_)); std(Dpos_res(:,frames_))];
            else
                frame_stats = [mean(INTpos_res); std(INTpos_res); mean(pos_res); std(pos_res); mean(Dpos_res); std(Dpos_res)];
            end
        end
        function alpha_mat = compute_alpha_logD(t_vec_, x_mat_)
            alpha_mat = swirl_group.compute_D(log(t_vec_), log(x_mat_));
            % [M,N] = size(x_mat_);
            % tlog = log(t_vec_);
            % xlog = log(x_mat_);
            % alpha_mat = nan(M,N);
            % for i = 1:M
            %     alpha_mat(i,:) = approx_deriv_weighted_central(tlog,xlog(i,:));
            % end
        end
        function [pos, dish] = cellrow2posdish(gp_cell_, i_, beads, len_pos)
            pos=reshape(gp_cell_{i_,1},[beads,len_pos,size(gp_cell_{i_,1},2)]);
            dish=gp_cell_{i_,2};
        end
        function make_movie_comp(AYfig_in,movie_data,movie_specs,watch_)
            Frame_vec = movie_specs.Frame_vec;

            if nargin==4
                AYfig_in.init_movie(length(Frame_vec), watch_);
            else
                AYfig_in.init_movie(length(Frame_vec));
            end

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

            f = 1;
            for i=reshape(Frame_vec, 1, [])
                plot(movie_ax, wallv(:, 1)+dish(i, 2), wallv(:, 2)+dish(i, 3), 'k -')
                hold(AYfig_in.ax, 'on');
                dots = scatter(movie_ax, pos(:, 1, i), pos(:, 2, i), 'o', 'filled', 'CData', colors, 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.5);
                txtbx = annotation(movie_fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(i-1), 'LineStyle', 'none', 'FontSize', 16);
                hold(movie_ax, 'off');
                axis(movie_ax, lims);
                drawnow
                movie_gen(f) = getframe(movie_fig);
                delete(txtbx);
                f = f+1;
            end
            txtbx = annotation(movie_fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(i-1), 'LineStyle', 'none', 'FontSize', 16);
            AYfig_in.movie_gen = movie_gen;
        end
    end
end

function prime_vec = approx_deriv_weighted_central(t_in, x_in, k_)
    n = length(t_in);
    prime_vec = nan(size(x_in));

    if (nargin==3)
        k=k_;
    else
        k=3;
    end

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
