classdef event
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        relay_id;
        event_count;

        specs;

        records;
        params;
        post_params;

        global_indices;

        residuals;

        bead_events_full;

        bead_res;
        bead_alpha;

        post_data; 
    end
    methods
        function obj = event(dat_dir_name_, exp_name_, dat_name_, relay_id_, event_count_, specs)

            pool_data = struct('int_data', nan(specs.record_int_len, specs.npool), ...
                            'double_data', nan(specs.record_double_len, specs.npool), ...
                            'int_chunk_data', nan(specs.beads*specs.record_int_chunk_len, specs.npool), ...
                            'double_chunk_data', nan(specs.beads*specs.record_double_chunk_len, specs.npool));
            records = relay_record.empty(specs.npool, 0);
            params = nan(specs.param_len, specs.npool);
            post_params = nan(specs.param_len, specs.npool);

            dat_pre = fopen([dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_event' num2str(event_count_) '.redat']);
            event_frame_count = fread(dat_pre,[1, specs.beads*specs.Frames], 'int=>int');
            for i = 1:specs.npool
                pool_data.int_data(:,i) = fread(dat_pre,[specs.record_int_len,1], 'int=>int');
                pool_data.double_data(:,i) = fread(dat_pre,[specs.record_double_len,1], 'double=>double');
                pool_data.int_chunk_data(:,i) = fread(dat_pre,[specs.beads*specs.record_int_chunk_len,1], 'int=>int');
                pool_data.double_chunk_data(:,i) = fread(dat_pre,[specs.beads*specs.record_double_chunk_len,1], 'double=>double');

                params(:,i) = fread(dat_pre,[specs.param_len,1], 'double=>double');

                records(i) = relay_record(params(:,i));
            end
            fclose(dat_pre);

            dat_post = fopen([dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_postevent' num2str(event_count_) '.redat']);
            head = fread(dat_post,[1,2], 'int=>int');
            post_data = struct('int_vec', fread(dat_post,[head(1),1], 'int=>int'), ...
                            'double_vec', fread(dat_post,[head(2),1], 'double=>double'), ...
                            'end', fread(dat_post,[specs.beads,1], 'int=>int'));
            for i = 1:specs.npool
                post_params(:,i) = fread(dat_post,[specs.param_len,1], 'double=>double');
                records(i).post_params = post_params(:,i);
            end
            fclose(dat_post);

            %%% assignments

            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.relay_id = relay_id_;
            obj.event_count = event_count_;
            obj.specs = specs;

            obj.records = records;
            obj.params = params;
            obj.post_params = post_params;

            obj.global_indices = pool_data.int_data(1,:);

            obj.residuals = pool_data.double_data(1,:);

            obj.bead_events_full = pool_data.int_chunk_data;

            obj.bead_res = pool_data.double_chunk_data(1:specs.beads, :);
            obj.bead_alpha = pool_data.double_chunk_data(specs.beads:end, :);

            obj.post_data = post_data;
        end
    end
end
