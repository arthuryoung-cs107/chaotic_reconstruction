classdef race_data
  properties
    dat_dir_name;
    exp_name;
    dat_name;
    gen_count;

    race_id;

    F;
    lambda;

    specs;
    specs_len = 7;

    gen;
  end
  methods
    function obj = race_data(dat_dir_name_, exp_name_, dat_name_, race_id_)
        [obj.dat_dir_name obj.exp_name obj.dat_name] = deal(dat_dir_name_,exp_name_,dat_name_);
        dat = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.rc' num2str(race_id_) '_end.rcdat']);
        obj.specs = fread( dat,[1, obj.specs_len], 'int=>int');
        fclose(dat);
        obj.gen_count = obj.specs(1);
        obj.race_id = race_id_;
        obj.gen = generation.empty(obj.gen_count, 0);
        for i=1:obj.gen_count
            obj.gen(i) = generation(obj.dat_dir_name, obj.exp_name, obj.dat_name, i, obj.race_id);
        end
    end
    function ancestors_zero = get_gen_ancestors_zero(obj, gen)
        leader_count = gen.leader_count;
        ancestors_zero = record.empty(leader_count, 0);
        for i = 1:leader_count
            ancestors_zero(i) = obj.get_ancestor_zero(gen.records(i));
        end
    end
    function anc_0 = get_ancestor_zero(obj, rec)
        if rec.parent_count == 0
            anc_0=rec;
        else
            anc_0 = obj.get_ancestor_zero(obj.gen(rec.parent_gen+1).records(rec.parent_global_index+1)); % get parent
        end
    end
    function lineage = get_lineage_records(obj, rec)
        lin_count = rec.parent_count+1;
        lineage = record.empty(lin_count, 0);
        lineage(1) = rec;

        parent_gen = rec.parent_gen;
        parent_index = rec.parent_global_index;
        for i = 2:lin_count
            parent = obj.gen(parent_gen+1).records(parent_index+1);
            lineage(i) = parent;
            parent_gen = parent.parent_gen;
            parent_index = parent.parent_global_index;
        end
        lineage = flip(lineage);
    end
    function lin_mat = get_lineage_mat(obj, rec)
        lin_count = rec.parent_count+1;
        lin_mat = nan(lin_count,3);
        lin_mat(1,:) = [rec.global_index, rec.gen, rec.frscore];
        parent_gen = rec.parent_gen;
        parent_index = rec.parent_global_index;
        for i = 2:lin_count
            parent = obj.gen(parent_gen+1).records(parent_index+1);
            lin_mat(i,:) = [parent.global_index, parent.gen, parent.frscore];
            parent_gen = parent.parent_gen;
            parent_index = parent.parent_global_index;
        end
        lin_mat = flip(lin_mat);
    end
    function lineage_cell = get_gen_lineage_cell(obj, gen, k_)
        if (nargin == 3)
            k = k_;
        else
            k = gen.leader_count;
        end
        lineage_cell = cell([k, 1]);
        [~,I_k] = maxk(gen.frscores, k);
        j = 1;
        for i = reshape(I_k,1,[])
            rec = gen.records(i);
            lin_count = rec.parent_count+1;
            lin_mat = nan(lin_count,3);
            lin_mat(1,:) = [rec.global_index, rec.gen, rec.frscore];
            parent_gen = rec.parent_gen;
            parent_index = rec.parent_global_index;
            for i = 2:lin_count
                parent = obj.gen(parent_gen+1).records(parent_index+1);
                lin_mat(i,:) = [parent.global_index, parent.gen, parent.frscore];
                parent_gen = parent.parent_gen;
                parent_index = parent.parent_global_index;
            end
            lineage_cell{j} = flip(lin_mat);
            j=j+1;
        end
    end
    function plot_final_lineage(obj, AYfig_, stat)
        red5 = [1 0 0];
        blue5 = [0 0 1];
        green4 = [0.1000 0.4440 0.2440];
        ax_ = AYfig_.ax_tile;

        k = 10;

        swtrue = stat.sw0;
        len_par = swtrue.len_par-2;
        par_true = swtrue.params(3:end);
        par_err_stat = stat.par_err(3:end,:)./abs(par_true);


        gen_count = obj.gen_count;
        gen_final = obj.gen(gen_count);
        gen_zero = obj.gen(1);
        rec_best = gen_final.get_best_record;
        ancestors_zero =  obj.get_gen_ancestors_zero(gen_final);
        best_lineage = obj.get_lineage_mat(rec_best);
        top_k_lineage = obj.get_gen_lineage_cell(gen_final,k);
        par_err_zero = (-gen_zero.param_mat+par_true)./abs(par_true);
        par_err_final = (-gen_final.param_mat+par_true)./abs(par_true);

        gen_vec = nan(gen_count, 1);
        index_vec = nan(gen_count, 1);
        parent_count_vec = nan(gen_count, 1);
        for i = 1:gen_count
            gen_vec(i) = ancestors_zero(i).gen;
            index_vec(i) = ancestors_zero(i).global_index;
            parent_count_vec(i) = gen_final.records(i).parent_count;
        end





        hold(ax_(1), 'on')
        for i=1:stat.len_gp
            plot(ax_(1), 1:len_par, par_err_stat(:,i), ' -', 'Color', [green4, 0.1], 'LineWidth', 1);
        end
        plot(ax_(1), 1:len_par, par_err_stat(:,stat.I_leader(1)), ' -', 'Color', [0 0 0], 'LineWidth', 1);
        xlabel(ax_(1), 'parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylabel(ax_(1), '$$(\hat{P}_i - P_i)/|\hat{P}_i|$$', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylim(ax_(1), [-4 1])

        hold(ax_(2), 'on')
        for i=1:gen_final.leader_count
            plot(ax_(2), 1:len_par, par_err_zero(:,i), ' -', 'Color', [green4, 0.1], 'LineWidth', 1);
        end
        plot(ax_(2), 1:len_par, par_err_zero(:,gen_zero.bleader_index+1), ' -', 'Color', blue5, 'LineWidth', 1);
        xlabel(ax_(2), 'parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylabel(ax_(2), '$$(\hat{P}_i - P_i)/|\hat{P}_i|$$', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylim(ax_(2), [-4 1])

        hold(ax_(3), 'on')
        for i=1:gen_final.leader_count
            plot(ax_(3), 1:len_par, par_err_final(:,i), ' -', 'Color', [green4, 0.1], 'LineWidth', 1);
        end
        plot(ax_(3), 1:len_par, par_err_final(:,gen_final.bleader_index+1), ' -', 'Color', red5, 'LineWidth', 1);
        xlabel(ax_(3), 'parameter index', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylabel(ax_(3), '$$(\hat{P}_i - P_i)/|\hat{P}_i|$$', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylim(ax_(3), [-4 1])

        hold(ax_(4), 'on')
        for i = 1:k
            top_i = top_k_lineage{i};
            plot(ax_(4), top_i(:,2), top_i(:,1), ' - o', 'Color', [0 0 0], 'LineWidth', 1);
        end
        plot(ax_(4), best_lineage(:,2), best_lineage(:,1), ' - p', 'Color', red5, 'LineWidth', 2);
        xlabel(ax_(4), 'generation', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylabel(ax_(4), 'pool index', 'Interpreter', 'Latex', 'Fontsize', 14)

        histogram(ax_(5), index_vec);
        xlabel(ax_(5), 'Global index', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylabel(ax_(5), 'frequency', 'Interpreter', 'Latex', 'Fontsize', 14)

        histogram(ax_(6), parent_count_vec);
        xlabel(ax_(6), 'Gen final parent count', 'Interpreter', 'Latex', 'Fontsize', 14)
        ylabel(ax_(6), 'frequency', 'Interpreter', 'Latex', 'Fontsize', 14)
    end
    function plot_diagnosticsi(obj, AYfig_, i_)
        obj.gen(i_).plot_diagnostics(AYfig_, obj.F, obj.lambda);
    end
    function sw_out = spawn_swbest(obj)
        swbname = ['swirl_best_rc' num2str(obj.race_id) '.odr/'];
        sw_out = ODR_data.construct_swirl([obj.dat_dir_name obj.exp_name], swbname, obj.dat_name);
        sw_out.params = AYdata.aysml_read([obj.dat_dir_name obj.exp_name swbname obj.dat_name '.sparam_best']);
        sw_out.len_par = length(sw_out.params);
    end
  end
end

function lin_mat = linrec2linmat(lineage_)
    lin_count = size(lineage_,2);
    lin_mat = nan(lin_count, 2);
    for i = 1:lin_count
        rec = lineage_(i);
        lin_mat(i,:) = [rec.global_index, rec.gen];
    end
end
