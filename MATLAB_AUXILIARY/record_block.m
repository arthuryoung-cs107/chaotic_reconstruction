classdef record_block
    properties
        ilen;
        dlen;
        ichunk_len;
        dchunk_len;

        ulen;
        int_mat;
        dub_mat;
        ichunk_mat;
        dchunk_mat;
        u_mat;
    end
    methods
        function obj = record_block(file_, ulen_,int_mat_,dub_mat_,ichunk_mat_,dchunk_mat_,u_mat_)
            if (nargin==2)
                head=fread(file_,6,'int=>int')
                [hlen,len,ilen,dlen,ichunk_len,dchunk_len]=deal(head(1),head(2),head(3),head(4),head(5),head(6))
                ulen_

                int_mat=nan(ilen,len);
                dub_mat=nan(dlen,len);
                ichunk_mat=nan(ichunk_len,len);
                dchunk_mat=nan(dchunk_len,len);
                u_mat=nan(ulen_,len);
                for i = 1:len
                    int_vec=fread(file_,[ilen 1],'int=>int')
                    dub_vec=fread(file_,[dlen 1],'double=>double')
                    ichunk_vec=fread(file_,[ichunk_len 1],'int=>int')
                    dchunk_vec=fread(file_,[dchunk_len 1],'double=>double')
                    % int_mat(:,i)=fread(file_,[ilen 1],'int=>int');
                    % dub_mat(:,i)=fread(file_,[dlen 1],'double=>double');
                    % ichunk_mat(:,i)=fread(file_,[ichunk_len 1],'int=>int');
                    % dchunk_mat(:,i)=fread(file_,[dchunk_len 1],'double=>double');
                    % u_mat(:,i)=fread(file_,[ulen_ 1],'double=>double');
                    pause
                end
            else
                ilen=size(int_mat_,1);
                dlen=size(dub_mat_,1);
                ichunk_len=size(ichunk_mat_,1);
                dchunk_len=size(dchunk_mat_,1);
                int_mat=int_mat_;
                dub_mat=dub_mat_;
                ichunk_mat=ichunk_mat_;
                dchunk_mat=dchunk_mat_;
                u_mat=u_mat_;
            end
            obj.ilen=ilen;
            obj.dlen=dlen;
            obj.ichunk_len=ichunk_len;
            obj.dchunk_len=dchunk_len;
            obj.ulen=ulen_;

            obj.int_mat=int_mat;
            obj.dub_mat=dub_mat;
            obj.ichunk_mat=ichunk_mat;
            obj.dchunk_mat=dchunk_mat;
            obj.u_mat=u_mat;
        end
    end
end
