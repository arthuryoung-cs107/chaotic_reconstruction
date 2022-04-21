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
        function [pos, dish] = cellrow2posdish(gp_cell_, i_, beads, len_pos)
            pos=reshape(gp_cell_{i_,1},[beads,len_pos,size(gp_cell_{i_,1},2)]);
            dish=gp_cell_{i_,2};
        end
        function make_movie_comp(AYfig_in,movie_data,movie_specs,watch_)
            frames = movie_specs.Frames;
            if nargin==4
                AYfig_in.init_movie(frames, watch_);
            else
                AYfig_in.init_movie(frames);
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

            for i=1:frames
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
            txtbx = annotation(movie_fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(frames-1), 'LineStyle', 'none', 'FontSize', 16);
            AYfig_in.movie_gen = movie_gen;
        end
    end
end
