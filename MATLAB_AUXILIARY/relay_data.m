classdef relay_data
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        relay_id;

        specs;
        event0;

        gen_ind;

        gen;

        gen_end_data;
    end
    methods
        function obj = relay_data(dat_dir_name_, exp_name_, dat_name_, relay_id_)
            dat_start = fopen([dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_startup.redat']);
            header_start = fread(dat_start,[1, 2], 'int=>int');
            int_params_start = fread(dat_start,[1, header_start(1)], 'int=>int');
            double_params_start = fread(dat_start,[1, header_start(2)], 'double=>double');
            fclose(dat_start);

            specs = struct('gen_max', int_params_start(1), ...
            'nlead', int_params_start(2), ...
            'npool', int_params_start(3), ...
            'param_len', int_params_start(4), ...
            'beads', int_params_start(5), ...
            'Frames', int_params_start(6), ...
            'record_int_len', int_params_start(7), ...
            'record_double_len', int_params_start(8), ...
            'record_int_chunk_len', int_params_start(9), ...
            'record_double_chunk_len', int_params_start(10));

            event0 = event(dat_dir_name_, exp_name_, dat_name_, relay_id_, 0, specs);

            %%% assignments

            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.relay_id = relay_id_;

            obj.specs = specs;
            obj.event0 = event0;

        end
        function [gen_ind,gen,gen_end_struct] = read_gen_end(obj, gen_ind_)
            dat_end = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.re' num2str(obj.relay_id) '_end.redat']);
            header_end = fread(dat_end,[1, 2], 'int=>int');
            int_params_end = fread(dat_end,[1, header_end(1)], 'int=>int');
            double_params_end = fread(dat_end,[1, header_end(2)], 'double=>double');
            fclose(dat_end);

            gen_end_struct = struct('int_params', int_params_end, 'double_params', double_params_end);
            gen_last = int_params_end(1);
            if (nargin==1)
                gen_ind_set = 1:gen_last;
            else
                gen_ind_set = gen_ind_;
            end
            [gen_ind, gen] = obj.read_generations(gen_ind_set,gen_last);
        end
        function [gen_ind,gen] = read_generations(obj,gen_ind_, gen_last_)
            gen = relay_generation.empty(gen_last_, 0);
            gen_ind = reshape(gen_ind_, 1, []);
            for i = 1:gen_last_;
                gen(i) = relay_generation();
            end
            for i = gen_ind
                gen(i) = gen(i).read_generation_data(obj.dat_dir_name, obj.exp_name, obj.dat_name, obj.relay_id, i, obj.specs);
            end
        end
        function [lin_mat, params_out] = get_lineage_mat(obj, rec)
            lin_count = rec.parent_count;
            lin_mat = nan(lin_count,2);
            lin_mat(1,:) = [rec.global_index, rec.gen];
            parent_gen = rec.parent_gen;
            parent_index = rec.parent_global_index;

            if (nargout==1)
                for i = 2:lin_count
                    parent = obj.gen(parent_gen+1).rec(parent_index+1);
                    lin_mat(i,:) = [parent.global_index, parent.gen];
                    parent_gen = parent.parent_gen;
                    parent_index = parent.parent_global_index;
                end
            else
                len_par = length(rec.params);
                params_out = nan(len_par,lin_count);
                params_out(:, 1) = rec.params;
                for i = 2:lin_count
                    parent = obj.gen(parent_gen+1).rec(parent_index+1);
                    lin_mat(i,:) = [parent.global_index, parent.gen];
                    params_out(:,i) = parent.params;
                    parent_gen = parent.parent_gen;
                    parent_index = parent.parent_global_index;
                    if parent_gen<1
                        break
                    end
                end
                params_out = flip(params_out,2);
            end
            lin_mat = flip(lin_mat);
        end
        function write_relay_test(obj, params_, relay_id_, test_)
            header=size(params_);
            file_id = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.re' num2str(relay_id_) '_test' num2str(test_) '.redat'], 'w+');
            fwrite(file_id, header, 'int');
            fwrite(file_id, params_(:), 'double');
            fclose(file_id);
        end
        function test_out = read_relay_test(obj, swtrue_, test_id_, relay_id_)
            test_out = detection_group(swtrue_, obj.dat_dir_name, obj.exp_name, obj.dat_name, test_id_, relay_id_);
        end
    end
end
