classdef relay_generation
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        gen_count;
        relay_id;

        gen_specs;
        lead_dup_count;
        sample_weights;
        lead_par_w_mean;
        lead_par_w_var;

        net_residuals;
        pis;
        zetas;

        rec;
        leader_data;
        params;
    end
    methods
        function obj = relay_generation()

        end
        function obj = read_generation_data(obj, dat_dir_name_, exp_name_, dat_name_, relay_id_, gen_count_, specs)
            dat = fopen([dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_gen' num2str(gen_count_) '.redat']);
            header = fread(dat,[1, 2], 'int=>int');
            int_vec = fread(dat,[1, header(1)], 'int=>int');
            double_vec = fread(dat,[1, header(2)], 'double=>double');

            gen_specs = struct('gen_count', int_vec(1), ...
                               'leader_count',int_vec(2), ...
                               'pool_success_count',int_vec(3), ...
                               'pool_candidates',int_vec(4), ...
                               'best_leader',int_vec(5), ...
                               'worst_leader',int_vec(6), ...
                               'repl_count',int_vec(7), ...
                               'dup_count',int_vec(8), ...
                               'dup_unique',int_vec(9), ...
                               'resample_count',int_vec(10), ...
                               'residual_best', double_vec(1), ...
                               'residual_worst', double_vec(2), ...
                               'root_res_best', double_vec(3), ...
                               'root_res_worst', double_vec(4), ...
                               'pos_res_global', double_vec(5), ...
                               'lambda', double_vec(6), ...
                               'beta', double_vec(7), ...
                               'w_sum', double_vec(8), ...
                               'w_max', double_vec(9), ...
                               't_wheels', double_vec(10));

            lead_dup_count = fread(dat, [1,gen_specs.leader_count], 'int=>int');
            sample_weights = fread(dat, [1,gen_specs.leader_count], 'double=>double');
            lead_par_w_mean = fread(dat, [specs.param_len, 1], 'double=>double');
            lead_par_w_var = fread(dat, [specs.param_len, 1], 'double=>double');

            rec = relay_record.empty(gen_specs.leader_count, 0);
            leader_data = struct('int_data', nan(specs.record_int_len, specs.nlead), ...
                                'double_data', nan(specs.record_double_len, specs.nlead), ...
                                'double_chunk', nan(specs.beads*(specs.record_double_chunk_len-1), specs.nlead));
            params = nan(specs.param_len, specs.nlead);
            for i = 1:specs.nlead
                leader_data.int_data(:,i) = fread(dat, [specs.record_int_len,1], 'int=>int');
                leader_data.double_data(:, i) = fread(dat, [specs.record_double_len,1], 'double=>double');
                leader_data.double_chunk(:,i) = fread(dat, [specs.beads*(specs.record_double_chunk_len-1),1], 'double=>double');
                params(:,i) = fread(dat, [specs.param_len,1], 'double=>double');
                rec(i) = relay_record(params(:,i));
            end
            fclose(dat);

            %%% assignments

            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.relay_id = relay_id_;

            obj.gen_count = gen_count_;

            obj.gen_specs=gen_specs;
            obj.lead_dup_count = lead_dup_count;
            obj.sample_weights = sample_weights;
            obj.lead_par_w_mean = lead_par_w_mean;
            obj.lead_par_w_var = lead_par_w_var;

            obj.rec = rec;
            obj.leader_data = leader_data;
            obj.params = params;
        end
    end
end
