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
        function compute_error()
          % obj.par_err = nan(obj.noise_len-1, obj.par_len);
          % obj.pos_err = nan(obj.noise_len-1, obj.sw(1).Frames-1);
          % for i=1:(obj.noise_len-1)
          %   obj.par_err(i, :) = (obj.pars(:, 1)-obj.pars(:, i+1))./abs(obj.pars(:, 1));
          %   obj.pos_err(i, :) = obj.sw(1).comp_pos_err(obj.sw(i+1));
          % end
          % [obj.par_err_truest, obj.I_truest] = mink(sum(abs(obj.par_err)'), obj.noise_len-1);
          % obj.pos_err_sum = sum(obj.pos_err');
          % [obj.pos_err_best, obj.I_best ] = mink(obj.pos_err_sum, obj.noise_len-1);
          % obj.par_cov = (abs(obj.par_err)).*(obj.pos_err_sum)';
        end
        function make_moviei(obj, AYfig_in, i )
            obj.sw(i).make_movie(AYfig_in);
        end
        function make_movieij(obj, AYfig_in, i, j)
            obj.sw(i).make_movie_comp(AYfig_in, obj.sw(j));
        end
        function make_movieijk(obj, AYfig_in, i, j, k)
            obj.sw(i).make_movie_comp2(AYfig_in, obj.sw(j), obj.sw(k));
        end
    end
end
