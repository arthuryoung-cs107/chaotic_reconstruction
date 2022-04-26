classdef event
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        relay_id;
        event_count;

        beads;
        Frames;
    end
    methods
        function obj = event(dat_dir_name_, exp_name_, dat_name_, relay_id_, event_count_, beads_, Frames_, record_int_len_, record_double_len_, record_int_chunk_count_, record_double_chunk_count_)
            dat_pre = fopen([dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_event' num2str(event_count_) '.redat']);
            event_frame_count = fread(dat_pre,[1, beads_*Frames_], 'int=>int');


            %%% assignments

            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.relay_id = relay_id_;
        end
    end
end
