classdef stat_data < swirl_group
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        noise_len;

        %% init_statistics
        par_err;
        pos_err;
        pos_err_acc;
        par_mean;

        I_best;
        I_truest;
        par_cov;

        %% compute_frscores
        frscores;
      end
    methods
        function obj = stat_data(dat_dir_name_, exp_name_, dat_name_, noise_len_)
            %% in C/C++ memory, the parameters are stored as [p0; p1; ... p_last]. AYlinalg invokes as [p0, ..., p_last]
            name = [dat_dir_name_ exp_name_ dat_name_ '_pars'];
            dims = dlmread([name '.aysml']);
            len_par = dims(1);
            id = fopen([name '.aydat']);
            params_true=fread(id,[len_par,1],'double=>double');
            if (nargin==4)
                noise_len=noise_len_;
            else
                noise_len=dims(2)-1;
            end
            params_mat=fread(id,[len_par,noise_len],'double=>double');
            fclose(id);

            rysml = ODR_data.get_rysml(dat_dir_name_, exp_name_, [dat_name_ '.0']);
            [Frames, beads, len_dish, len_pos] = deal(rysml.Frames, rysml.beads, rysml.len_dish, rysml.len_pos);
            mat_dims = [beads*len_pos, Frames];
            swtrue = ODR_data.construct_swirl(dat_dir_name_, exp_name_, [dat_name_ '.0']);
            swtrue.len_par = len_par;
            swtrue.params = params_true;

            dish_mat = nan(len_dish*Frames,noise_len);
            pos_mat = nan(len_pos*beads*Frames, noise_len);
            for i=1:noise_len
                id = fopen([rysml.dat_dir rysml.exp_name [dat_name_ '.' num2str(i)] '.rydat']);
                data = reshape(fread(id,'double'), [len_dish+beads*len_pos, Frames]);
                fclose(id);
                dish_mat(:,i) = reshape(data(1:len_dish, :), [len_dish*Frames, 1]);
                pos_mat(:,i) = reshape(permute(reshape(data(len_dish+1:end, :), [len_pos, beads, Frames]),[2 1 3]),[len_pos*beads*Frames, 1]);
            end
            gp_cell = cell([noise_len, 2]);
            gp_cell(:,1) = reshape(cellfun(@(pos) reshape(pos, [len_pos*beads, Frames]), num2cell(pos_mat,1), 'UniformOutput',false), [noise_len, 1]);
            gp_cell(:,2) = reshape(cellfun(@(dish) reshape(dish, [len_dish, Frames]), num2cell(dish_mat,1),   'UniformOutput',false), [noise_len, 1]);

            obj@swirl_group(swtrue, gp_cell, noise_len, Frames*ones(noise_len, 1), params_mat);

            obj.par_err = -params_mat+params_true;
            obj.pos_err = cell2mat(cellfun(@(pos) frame_pos_err(swtrue.pos(:,1:2,:),reshape(pos,[beads,len_pos,Frames])), gp_cell(:,1), 'UniformOutput', false));
            obj.pos_err_acc = sum(obj.pos_err, 2);
            obj.par_mean = mean(params_mat, 2);
            [~,obj.I_best]=mink(obj.pos_err_acc, noise_len);
            [~,obj.I_truest]=mink(sum(abs(obj.par_err)./abs(params_true), 1), noise_len);
            obj.par_cov=cov(params_mat');

            %% set stat data variables
            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.noise_len = noise_len;
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
        function stat = init_stat_generation(stat_, stgp_)
            stat = stat_;
            stgp = stgp_;
            stat.frscores = stgp.compute_frscores();
        end
    end
end

function pos_err_out = frame_pos_err(pos0,posi)
    diff = pos0-posi(:,1:2,:);
    pos_err_out = sum(sqrt(squeeze(sum(diff.*diff, 2)))./sqrt(squeeze(sum(pos0.*pos0, 2))), 1);
end
