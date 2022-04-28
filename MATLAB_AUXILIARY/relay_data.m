classdef relay_data
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        relay_id;

        specs;
        event0;

        gen;

        test;
    end
    methods
        function obj = relay_data(dat_dir_name_, exp_name_, dat_name_, relay_id_)
            dat_start = fopen([dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_startup.redat']);
            header_start = fread(dat_start,[1, 2], 'int=>int');
            int_params_start = fread(dat_start,[1, header_start(1)], 'int=>int');
            double_params_start = fread(dat_start,[1, header_start(2)], 'double=>double');
            fclose(dat_start);

            dat_end = fopen([dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_end.redat']);
            header_end = fread(dat_end,[1, 2], 'int=>int');
            int_params_end = fread(dat_end,[1, header_end(1)], 'int=>int');
            double_params_end = fread(dat_end,[1, header_end(2)], 'double=>double');
            fclose(dat_end);

            specs = struct('gen_max', int_params_start(1), ...
            'nlead', int_params_start(2), ...
            'npool', int_params_start(3), ...
            'param_len', int_params_start(4), ...
            'beads', int_params_start(5), ...
            'Frames', int_params_start(6), ...
            'record_int_len', int_params_start(7), ...
            'record_double_len', int_params_start(8), ...
            'record_int_chunk_len', int_params_start(9), ...
            'record_double_chunk_len', int_params_start(10),  ...
            'gen_last', int_params_end(1));

            event0 = event(dat_dir_name_, exp_name_, dat_name_, relay_id_, 0, specs);

            gen = relay_generation.empty(specs.gen_last, 0);

            for i = 1:specs.gen_last
                gen(i) = relay_generation(dat_dir_name_, exp_name_, dat_name_, relay_id_, i, specs);
            end

            %%% assignments

            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.relay_id = relay_id_;

            obj.specs = specs;
            obj.event0 = event0;
            obj.gen = gen;
        end
        function write_relay_test(obj, params_, test_)
            header=size(params_);
            file_id = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.re' num2str(obj.relay_id) '_test' num2str(test_) '.redat'], 'w+');
            fwrite(file_id, header, 'int');
            fwrite(file_id, params_(:), 'double');
            fclose(file_id);
        end
        function test_out = read_relay_test(obj, test_, relay_id_)
            dat_test_in_name = [obj.dat_dir_name obj.exp_name obj.dat_name '.re' num2str(relay_id_) '_test' num2str(test_) '.redat']
            dat_test_in = fopen(dat_test_in_name);
            header_in = fread(dat_test_in, [1,2], 'int=>int');
            [param_len,npool] = deal(header_in(1),header_in(2));
            params = fread(dat_test_in,[param_len,npool], 'double=>double');
            fclose(dat_test_in);


            test_directory = [obj.dat_dir_name obj.exp_name obj.dat_name '.re' num2str(relay_id_) '_test' num2str(test_) '_results/'];

            dat_test_specs = fopen([test_directory, 'results_specs.redat']);
            results_specs = fread(dat_test_specs, [1,3], 'int=>int');
            [Frame_end, beads, dof] = deal(results_specs(1), results_specs(2), results_specs(3));
            fclose(dat_test_specs);

            accres_vec = nan(npool, 1);
            [sim_pos_mat, pos_mat, ref_pos_mat] = deal(nan(Frame_end*beads*dof, npool));
            [pos_res_mat, alpha_mat, INTmat] = deal(nan(Frame_end*beads, npool));
            for i = 1:npool
                dat = fopen([test_directory 'par' num2str(i-1) '.redat']);
                header = fread(dat, [1,2], 'int=>int');
                sim_pos_mat(:,i) = fread(dat, [Frame_end*beads*dof,1], 'double=>double');
                pos_mat(:,i) = fread(dat, [Frame_end*beads*dof,1], 'double=>double');
                ref_pos_mat(:,i) = fread(dat, [Frame_end*beads*dof,1], 'double=>double');
                pos_res_mat(:,i) = fread(dat, [Frame_end*beads,1], 'double=>double');
                alpha_mat(:,i) = fread(dat, [Frame_end*beads,1], 'double=>double');
                INTmat(:,i) = fread(dat, [Frame_end*beads,1], 'double=>double');
                fclose(dat);
            end

            test_out = struct('param_len', param_len,'npool', npool,'Frame_end', Frame_end,'beads', beads,'dof', dof, 'sim_pos_mat', sim_pos_mat, 'pos_mat', pos_mat, 'ref_pos_mat', ref_pos_mat, 'pos_res_mat', pos_res_mat, 'alpha_mat', alpha_mat, 'INTmat', INTmat);
        end
    end
end
