classdef record_block
    properties
        int_mat;
        dub_mat;
        ichunk_mat;
        dchunk_mat;
    end
    methods
        function obj = record_block(file_)
            head=fread(file_,6,'int=>int');
            [hlen,len,ilen,dlen,ichunk_len,dchunk_len]=deal(head(1),head(2),head(3),head(4),head(5),head(6));
            
            int_mat=nan(ilen,len);
            dub_mat=nan(dlen,len);
            ichunk_mat=nan(ichunk_len,len);
            dchunk_mat=nan(dchunk_len,len);
            for i = 1:len
                int_mat(:,i)=fread(file_,[ilen 1],'int=>int');
                dub_mat(:,i)=fread(file_,[dlen 1],'double=>double');
                ichunk_mat(:,i)=fread(file_,[ichunk_len 1],'int=>int');
                dchunk_mat(:,i)=fread(file_,[dchunk_len 1],'double=>double');
            end
            obj.int_mat=int_mat;
            obj.dub_mat=dub_mat;
            obj.ichunk_mat=ichunk_mat;
            obj.dchunk_mat=dchunk_mat;
        end
    end
end
