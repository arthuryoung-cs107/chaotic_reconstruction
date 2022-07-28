classdef MH_genetic_data
    properties
        test_dir_name;
        MHpar;

        %% constant system values
        ulen;
        nbeads;
        Frames;
        nlead;
        npool;
        dt_sim;
        t_phys;
        sigma;
    end
    methods
        function obj = MH_genetic_data(nbeads_, par_id_, MH_id_)
            test_dir_name = ['../dat_dir/swirl' num2str(nbeads_) '.odr/pts.' num2str(par_id_) '.MH' num2str(MH_id_) '_results/'];

            startspecs_file = fopen([test_dir_name 'startspecs.mhdat']);
            MHpar = MH_genetic_data.read_MH_params(startspecs_file);
            gen_header = fread(startspecs_file,3,'int=>int');
            ints=fread(startspecs_file,gen_header(2),'int=>int');
            dubs=fread(startspecs_file,gen_header(3),'double=>double');
            fclose(startspecs_file);

            %% assignments
            obj.test_dir_name = test_dir_name;
            obj.MHpar = MHpar;

            obj.ulen=MHpar.ulen;
            obj.nbeads=MHpar.nbeads;
            obj.Frames=MHpar.Frames;
            obj.nlead=MHpar.nlead;
            obj.npool=MHpar.npool;
            obj.dt_sim=MHpar.dt_sim;
            obj.t_phys=MHpar.t_phys;
            obj.sigma=MHpar.sigma;
        end
    end

    methods (Static)
        function MHpar_out = read_MH_params(file_)
            header=fread(file_, [1 3], 'int=>int');
            [hlen,ilen,dlen]=deal(header(1),header(2),header(3));
            ints=fread(file_,ilen,'int=>int');
            dubs=fread(file_,dlen,'double=>double');
            MHpar_out = struct( 'ulen', ints(1), ...
                                'nbeads', ints(2), ...
                                'Frames', ints(3), ...
                                'nlead', ints(4), ...
                                'npool', ints(5), ...
                                'dt_sim', dubs(1), ...
                                't_phys', dubs(2), ...
                                'sigma', dubs(3));
        end
        function event_block_out = read_event_block_i(obj, i_, test_dir_name_, MHpar_)
            if nargin==2
                test_dir_name=obj.test_dir_name;
                nbeads = obj.MHpar.nbeads;
                ulen = obj.MHpar.ulen;
            else
                test_dir_name=test_dir_name_;
                nbeads = MHpar_.nbeads;
                ulen = MHpar_.ulen;
            end

            file = fopen([test_dir_name 'event_block' num2str(i_) '.mhdat']);
            event_block_out = event_block(file,ulen,nbeads);
            fclose(file);
        end
        function Class_out = read_Class_i(obj, i_, test_dir_name_, MHpar_)
            if nargin==2
                test_dir_name=obj.test_dir_name;
                ulen = obj.MHpar.ulen;
            else
                test_dir_name=test_dir_name_;
                ulen = MHpar_.ulen;
            end

            file = fopen([test_dir_name 'Class' num2str(i_) '.mhdat']);
            Class_out = Class_block(file,ulen);
            fclose(file);
        end
        function gen_out = read_gen_i(obj, i_, test_dir_name_, MHpar_)
            if nargin==2
                test_dir_name=obj.test_dir_name;
                ulen = obj.MHpar.ulen;
            else
                test_dir_name=test_dir_name_;
                ulen = MHpar_.ulen;
            end

            file = fopen([test_dir_name 'gen' num2str(i_) '.mhdat']);
            gen_out = generation_block(file,ulen);
            fclose(file);
        end
        function [ilen,dlen,ints,dubs] = read_int_double_chunk(file_)
            header=fread(file_,3,'int=>int');
            [ilen,dlen] = deal(header(2),header(3));
            ints=fread(file_,ilen,'int=>int');
            dubs=fread(file_,dlen,'double=>double');
        end
    end
end
