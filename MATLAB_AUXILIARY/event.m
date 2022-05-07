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

        global_indices;

        residuals;

        bead_events_full;

        bead_res;
        bead_alpha;

    end
    methods
        function obj = event(dat_dir_name_, exp_name_, dat_name_, relay_id_, event_count_, specs)

            if event_count_==0
                rec_len = specs.npool;
            else
                rec_len = specs.nlead;
            end

            rec_data = struct('int_data', nan(specs.record_int_len, rec_len), ...
                            'double_data', nan(specs.record_double_len, rec_len), ...
                            'int_chunk_data', nan(specs.beads*specs.record_int_chunk_len, rec_len), ...
                            'double_chunk_data', nan(specs.beads*specs.record_double_chunk_len, rec_len));

            records = relay_record.empty(rec_len, 0);
            params = nan(specs.param_len, rec_len);

            dat = fopen([dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_event' num2str(event_count_) '.redat']);
            head = fread(dat,[1,2], 'int=>int');
            event_data = struct('int_vec', fread(dat,[head(1),1], 'int=>int'), ...
                            'double_vec', fread(dat,[head(2),1], 'double=>double'), ...
                            'event_frame_count', fread(dat,[specs.beads*specs.Frames,1], 'int=>int'), ...
                            'event_frames', fread(dat,[specs.beads,1], 'int=>int'));
            for i = 1:rec_len
                rec_data.int_data(:,i) = fread(dat,[specs.record_int_len,1], 'int=>int');
                rec_data.double_data(:,i) = fread(dat,[specs.record_double_len,1], 'double=>double');
                rec_data.int_chunk_data(:,i) = fread(dat,[specs.beads*specs.record_int_chunk_len,1], 'int=>int');
                rec_data.double_chunk_data(:,i) = fread(dat,[specs.beads*specs.record_double_chunk_len,1], 'double=>double');

                params(:,i) = fread(dat,[specs.param_len,1], 'double=>double');

                records(i) = relay_record(params(:,i));
            end
            fclose(dat);

            %%% assignments

            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.relay_id = relay_id_;
            obj.event_count = event_count_;
            obj.specs = specs;

            obj.records = records;
            obj.params = params;
        end
    end
end
