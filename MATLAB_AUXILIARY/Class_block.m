classdef Class_block < MH_genetic_it_data
    properties
        ilen;
        dlen;
        ints;
        dubs;

        leaders;
    end
    methods
        function obj = Class_block(file_,ulen_)
            [ilen,dlen,ints,dubs]=MH_genetic_data.read_int_double_chunk(file_);
            leaders=record_block(file_,ulen_);

            %% assignments
            obj@MH_genetic_it_data(ints,dubs);

            obj.ilen=ilen;
            obj.dlen=dlen;
            obj.ints=ints;
            obj.dubs=dubs;
            obj.leaders=leaders;
        end
    end
end
