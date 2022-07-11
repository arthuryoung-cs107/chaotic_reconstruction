classdef MH_doctor_data
    properties
        test_dir_name;
        MHpar;

        Frames_test;

        TEST_refp;

        event_block0;

        test_block;
    end
    methods
        function obj = MH_doctor_data(nbeads_,par_id_,MH_id_,test_id_)
            test_dir_name = ['../dat_dir/swirl' num2str(nbeads_) '.odr/pts.' num2str(par_id_) '.MH' num2str(MH_id_) '_test' num2str(test_id_) '_results/'];

            startspecs_file = fopen([test_dir_name 'startspecs.mhdat']);
            MHpar = MH_genetic_data.read_MH_params(startspecs_file);
            % read MH_genetic stuff
            gen_header = fread(startspecs_file,3,'int=>int');
            ints=fread(startspecs_file,gen_header(2),'int=>int');
            dubs=fread(startspecs_file,gen_header(3),'double=>double');
            % read MH_doctor stuff
            doc_header = fread(startspecs_file,2,'int=>int');

            dof=MHpar.nbeads*2;
            Frames_test=doc_header(2);

            TEST_refp=fread(startspecs_file,dof*Frames_test,'double=>double');
            fclose(startspecs_file);

            event_block0=MH_genetic_data.read_event_block_i(0,0,test_dir_name,MHpar);

            obj.test_dir_name=test_dir_name;
            obj.MHpar=MHpar;
            obj.Frames_test=Frames_test;
            obj.TEST_refp=TEST_refp;
            obj.event_block0=event_block0;
        end
        function test_block_out = read_test_block(obj)
            % test_block_out = test_block(obj);
            test_block_out = test_block(obj,obj.event_block0.pool);
        end
    end
end
