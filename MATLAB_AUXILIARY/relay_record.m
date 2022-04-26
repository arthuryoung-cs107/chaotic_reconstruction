classdef relay_record
    properties
        global_index;
        gen;
        parent_gen;
        parent_count;
        parent_global_index;
        dup_count;

        res;
        event_res;
        weight;

        event_pos;

        params;
        res_data;
        alpha_data;

        post_params;
    end
    methods
        function obj = relay_record(params)
            obj.params = params;
        end
    end
end
