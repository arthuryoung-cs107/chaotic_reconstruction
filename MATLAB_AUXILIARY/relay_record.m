classdef relay_record
    properties
        global_index;
        gen;
        parent_gen;
        parent_count;
        parent_global_index;
        dup_count;

        residual;
        root_residual;
        net_residual;
        net_smooth_residual;
        weight;
        px;
        zeta;

        event_positions;

        smooth_residual;
        stiff_residual;
        alpha_data;

        params;
    end
    methods
        function obj = relay_record(params, int_data, double_data, double_chunk)
            obj.params = params;
            if nargin>1
                obj.global_index=int_data(1);
                obj.gen=int_data(2);
                obj.parent_gen=int_data(3);
                obj.parent_count=int_data(4);
                obj.parent_global_index=int_data(5);
                obj.dup_count=int_data(6);

                obj.residual=double_data(1);
                obj.root_residual=double_data(2);
                obj.net_residual=double_data(3);
                obj.net_smooth_residual=double_data(4);
                obj.weight=double_data(5);
                obj.px=double_data(6);
                obj.zeta=double_data(7);                
            end
        end
    end
end
