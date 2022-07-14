classdef MH_genetic_data
    properties
        test_dir_name;
        MHpar;

        event_block0;
        Class0;
        gen0;
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

            obj.test_dir_name = test_dir_name;
            obj.MHpar = MHpar;
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
            header = fread(file,3,'int=>int');
            [hlen,npool,stev_latest]=deal(header(1),header(2),header(3));

            nb_x_F = nbeads*stev_latest;

            event_block_out = struct(   'stev_comp',fread(file,nbeads,'int=>int'), ...
                                        'stev_ordered',fread(file,nbeads,'int=>int'), ...
                                        'comps_ordered',fread(file,nbeads,'int=>int'), ...
                                        'nev_state_comp',fread(file,nb_x_F,'int=>int'), ...
                                        'nobs_state_comp',fread(file,nb_x_F,'int=>int'), ...
                                        'rho2stable_comp',fread(file,nbeads,'double=>double'), ...
                                        'rho2unstable_comp',fread(file,nbeads,'double=>double'), ...
                                        'mur2_state_comp',fread(file,nb_x_F,'double=>double'), ...
                                        'stdr2_state_comp',fread(file,nb_x_F,'double=>double'), ...
                                        'mualpha_state_comp',fread(file,nb_x_F,'double=>double'), ...
                                        'stdalpha_state_comp',fread(file,nb_x_F,'double=>double'), ...
                                        'pool',record_block(file,ulen));
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
            header = fread(file,3,'int=>int');
            [hlen,ilen,dlen] = deal(header(1),header(2),header(3));
            Class_out = MH_genetic_data.read_int_double_chunk(file,ilen,dlen);

            Class_out.leaders = record_block(file,ulen);
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
            header = fread(file,3,'int=>int');
            [hlen,ilen,dlen] = deal(header(1),header(2),header(3));
            gen_out = MH_genetic_data.read_int_double_chunk(file,ilen,dlen);

            gen_out.u_mean=fread(file,ulen,'double=>double');
            gen_out.u_wmean=fread(file,ulen,'double=>double');
            gen_out.u_var=fread(file,ulen,'double=>double');
            gen_out.u_wvar=fread(file,ulen,'double=>double');
            % if (feof(file))
            %     gen_out.full_diagnostics=false;
            % else
            %     gen_out.full_diagnostics=true;
            %     gen_out.pool=record_block(file,ulen);
            % end
            fclose(file);
        end
        function chunk_data_out = read_int_double_chunk(file_,ilen_,dlen_)
            chunk_data_out = struct(    'ilen', ilen_, ...
                                        'dlen', dlen_, ...
                                        'ints', fread(file_,ilen_,'int=>int'), ...
                                        'dubs', fread(file_,dlen_,'double=>double'));
        end
    end
end
