classdef test_block
    properties
        Frames_test;
        MHpar;

        psim_mat;
        r2_state_comp;
        alpha_state_comp;
        TEST_p;
        TEST_INTr2;

        records;
    end
    methods
        function obj = test_block(doc_,records_)
            MHpar=doc_.MHpar;
            Frames_test=doc_.Frames_test;
            nbeads=MHpar.nbeads;
            npool=MHpar.npool;
            ulen=MHpar.ulen;
            dof=nbeads*2;

            psim_mat=nan(dof*Frames_test,npool);
            r2_state_comp=nan(nbeads*Frames_test,npool);
            alpha_state_comp=nan(nbeads*Frames_test,npool);
            TEST_p=nan(dof*Frames_test,npool);
            TEST_INTr2=nan(nbeads*Frames_test,npool);

            if (nargin==1)
                ilen=size(doc_.event_block0.pool.int_mat,1);
                dlen=size(doc_.event_block0.pool.dub_mat,1);
                ichunk_len=size(doc_.event_block0.pool.ichunk_mat,1);
                dchunk_len=size(doc_.event_block0.pool.dchunk_mat,1);

                int_mat=nan(ilen,npool);
                dub_mat=nan(dlen,npool);
                ichunk_mat=nan(ichunk_len,npool);
                dchunk_mat=nan(dchunk_len,npool);
                u_mat=nan(ulen,npool);
                for i = 1:npool
                    file_i = fopen([doc_.test_dir_name 'par' num2str(i-1) '.mhdat']);
                    psim_mat(:,i)=fread(file_i,dof*Frames_test,'double=>double');
                    r2_state_comp(:,i)=fread(file_i,nbeads*Frames_test,'double=>double');
                    alpha_state_comp(:,i)=fread(file_i,nbeads*Frames_test,'double=>double');
                    TEST_p(:,i)=fread(file_i,dof*Frames_test,'double=>double');
                    TEST_INTr2(:,i)=fread(file_i,nbeads*Frames_test,'double=>double');

                    head=fread(file_i,6,'int=>int');

                    int_mat(:,i)=fread(file_i,ilen,'int=>int');
                    dub_mat(:,i)=fread(file_i,dlen,'double=>double');
                    ichunk_mat(:,i)=fread(file_i,ichunk_len,'int=>int');
                    dchunk_mat(:,i)=fread(file_i,dchunk_len,'double=>double');
                    u_mat(:,i)=fread(file_i,ulen,'double=>double');

                    fclose(file_i);
                end
                records = record_block(0,ulen,int_mat,dub_mat,ichunk_mat,dchunk_mat,u_mat);
            else
                for i = 1:npool
                    file_i = fopen([doc_.test_dir_name 'par' num2str(i-1) '.mhdat']);
                    psim_mat(:,i)=fread(file_i,dof*Frames_test,'double=>double');
                    r2_state_comp(:,i)=fread(file_i,nbeads*Frames_test,'double=>double');
                    alpha_state_comp(:,i)=fread(file_i,nbeads*Frames_test,'double=>double');
                    TEST_p(:,i)=fread(file_i,dof*Frames_test,'double=>double');
                    TEST_INTr2(:,i)=fread(file_i,nbeads*Frames_test,'double=>double');
                    fclose(file_i);
                end
                records=records_;
            end
            obj.Frames_test=Frames_test;
            obj.MHpar=MHpar;

            obj.psim_mat=psim_mat;
            obj.r2_state_comp=r2_state_comp;
            obj.alpha_state_comp=alpha_state_comp;
            obj.TEST_p=TEST_p;
            obj.TEST_INTr2=TEST_INTr2;

            obj.records=records;
        end
    end
end
