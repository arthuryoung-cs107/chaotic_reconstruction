classdef event_block
    properties
        nbeads;
        npool;
        stev_latest;

        stev_comp;
        stev_ordered;
        comps_ordered;
        nev_state_comp;
        nobs_state_comp;
        rho2stable_comp;
        rho2unstable_comp;
        mur2_state_comp;
        stdr2_state_comp;
        mualpha_state_comp;
        stdalpha_state_comp;
        pool;

        r2_pool_Framebead;
        alpha_pool_Framebead;
    end
    methods
        function obj = event_block(file_,ulen_,nbeads_)
            hlen = fread(file_,1,'int=>int');
            header = fread(file_,hlen,'int=>int');
            [npool,stev_latest]=deal(header(1),header(2));
            if (hlen==3)
                read_pool_event_data=header(3);
            else
                read_pool_event_data=false;
            end
            nb_x_F = nbeads_*stev_latest;

            %% assignments
            obj.nbeads = nbeads_;
            obj.npool = npool;
            obj.stev_latest = stev_latest;

            obj.stev_comp=fread(file_,nbeads_,'int=>int');
            obj.stev_ordered=fread(file_,nbeads_,'int=>int');
            obj.comps_ordered=fread(file_,nbeads_,'int=>int');
            obj.nev_state_comp=fread(file_,nb_x_F,'int=>int');
            obj.nobs_state_comp=fread(file_,nb_x_F,'int=>int');
            obj.rho2stable_comp=fread(file_,nbeads_,'double=>double');
            obj.rho2unstable_comp=fread(file_,nbeads_,'double=>double');
            obj.mur2_state_comp=fread(file_,nb_x_F,'double=>double');
            obj.stdr2_state_comp=fread(file_,nb_x_F,'double=>double');
            obj.mualpha_state_comp=fread(file_,nb_x_F,'double=>double');
            obj.stdalpha_state_comp=fread(file_,nb_x_F,'double=>double');
            obj.pool=record_block(file_,ulen_);

            if (read_pool_event_data)
                obj.r2_pool_Framebead = fread(file_,npool*nb_x_F,'double=>double');
                obj.alpha_pool_Framebead = fread(file_,npool*nb_x_F,'double=>double');
            end
        end
        function [r2_p_F_b_out alpha_p_F_b_out] = get_pool_Frame_bead(obj)
            nbeads=obj.nbeads;
            Frames=obj.stev_latest;
            npool=obj.npool;

            r2_tens=reshape(obj.r2_pool_Framebead,[nbeads,Frames,npool]);
            r2_p_F_b_out = permute(r2_tens,[3,2,1]); %% bead,Frame,pool->pool,Frame,bead
            if (nargout==2)
                alpha_tens=reshape(obj.r2_pool_Framebead,[nbeads,Frames,npool]);
                alpha_p_F_b_out=permute(alpha_tens,[3,2,1]); %% bead,Frame,pool->pool,Frame,bead
            end
        end
    end
end
