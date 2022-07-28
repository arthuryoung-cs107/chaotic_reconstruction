classdef generation_block < MH_genetic_it_data
    properties (Constant)
        ilen_full_genetic_old=13;
        pool_flag_index=14;
    end
    properties
        ilen;
        dlen;
        ints;
        dubs;

        u_mean;
        u_wmean;
        u_var;
        u_wvar;

        pool;
    end
    methods
        function obj = generation_block(file_,ulen_)
            [ilen,dlen,ints,dubs]=MH_genetic_data.read_int_double_chunk(file_);
            u_mean=fread(file_,ulen_,'double=>double');
            u_wmean=fread(file_,ulen_,'double=>double');
            u_var=fread(file_,ulen_,'double=>double');
            u_wvar=fread(file_,ulen_,'double=>double');

            pool=generation_block.check_pool_data(file_,ulen_,ilen,ints);

            %% assignments
            obj@MH_genetic_it_data(ints,dubs);

            obj.ilen=ilen;
            obj.dlen=dlen;
            obj.ints=ints;
            obj.dubs=dubs;
            obj.u_mean=u_mean;
            obj.u_wmean=u_wmean;
            obj.u_var=u_var;
            obj.u_wvar=u_wvar;
            obj.pool=pool;
        end
    end
    methods (Static)
        function pool_out = check_pool_data(file_,ulen_,ilen_,ints_)
            if (ilen_>generation_block.ilen_full_genetic_old) %% check if we have an old version
                if (ints_(generation_block.pool_flag_index)) %% check if we should read for a pool
                    pool_out=record_block(file_,ulen_);
                else
                    pool_out=0;
                end
            else
                pool_out=0;
            end
        end
    end
end
