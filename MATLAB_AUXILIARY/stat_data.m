classdef stat_data < handle
    properties
        % friend class
        % swirl_group.m;

        swgp;
        swtrue;

        dat_dir_name;
        exp_name;
        dat_name;

        len_par;
        noise_len;

        % shared with swirl_group
        %{a cell object that contains each swirl's position data and dish data%}
        gp_cell;
        %{the parameters associated with each swirl in this swirl group%}
        params_mat;

        %% must be initialized
        par_err;
        pos_err;
        pos_err_acc;

        I_best;
        I_truest;
        par_cov;

        frscores;
      end
    methods
        function obj = stat_data(dat_dir_name_, exp_name_, dat_name_, noise_len_)
            params_mat = AYdata.aysml_read([dat_dir_name_ exp_name_ dat_name_ '_pars']);
            [len_par,noise_len] = size(params_mat);
            if (nargin==4)
                noise_len=noise_len_;
            else
                noise_len=noise_len-1;
            end
            params_true = params_mat(:, 1);
            params_mat = params_mat(:, 2:(noise_len+1));

            rysml = ODR_data.get_rysml(dat_dir_name_, exp_name_, [dat_name_ '.0']);
            swtrue = ODR_data.construct_swirl(dat_dir_name_, exp_name_, [dat_name_ '.0']);
            swtrue.len_par = len_par;
            swtrue.params = params_true;

            gp_cell = cell([noise_len, 2]);
            [Frames, beads, len_pos] = deal(rysml.Frames, rysml.beads, rysml.len_pos);
            mat_dims = [beads*len_pos, Frames];
            for i=1:noise_len
                rysml.dat_name = [dat_name_ '.' num2str(i)];
                rydat_i = ODR_data.get_rydat(rysml);
                [gp_cell{i, :}] = deal(reshape(rydat_i.pos, mat_dims), rydat_i.dish);
            end

            %% set stat data variables
            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;

            obj.len_par = len_par;
            obj.noise_len = noise_len;

            obj.swtrue = swtrue;

            obj.gp_cell = gp_cell;
            obj.params_mat = params_mat;
        end

        function swgp_out = spawn_swirlgroup(obj)
            swgp_out = swirl_group(obj.swtrue, obj.gp_cell, obj.noise_len, obj.swtrue.Frames*ones(obj.noise_len, 1), obj.params_mat);
        end
        function swout = spawn_swindex(obj, id_)
            swout = ODR_data(obj.dat_dir_name, obj.exp_name, [obj.dat_name '.' num2str(i)]);
        end
        function plot_gen_scores(obj, ax_)
            frscores = obj.frscores;
            fcount = histcounts(frscores, max(frscores)-min(frscores)+1);
            generation.plot_gen_scores(ax_, frscores, fcount);
            

        end
        function make_POVray_inputsi(obj, i_, pov_dir_)
            obj.sw(i_).make_POVray_inputs(pov_dir_, obj.dat_name);
        end
    end
    methods (Static)
        function stat = init_statistics(stat_, stgp_)
            stat = stat_;
            stgp = stgp_;

            [stat.par_err, stat.pos_err, stat.pos_err_acc] = stgp.compute_error();
            [stat.I_best, stat.I_truest, stat.par_cov] = stgp.compute_statistics(stat.par_err, stat.pos_err, stat.pos_err_acc);
        end
        function stat = init_stat_generation(stat_, stgp_)
            stat = stat_;
            stgp = stgp_;
            stat.frscores = stgp.compute_frscores();
        end
    end
end
