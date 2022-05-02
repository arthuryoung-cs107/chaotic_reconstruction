classdef detection_group
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        relay_id;
        test_id;

        swtrue;

        param_len;
        npool;
        params;

        Frame_end;
        beads;
        dof;
        Frame_full;
        ref_pos_mat;

        earlylate_mat;
        acc_res_mat;
        sim_pos_mat;
        pos_mat;
        pos_res_mat;
        alpha_mat;
        INTmat;
        pos_res_re_mat;
        alpha_re_mat;
        kill_frame_mat;
        alpha_kill_mat;

        sim_beadres;
        beadINTres;
      end
    methods
        function obj = detection_group(sw_, dat_dir_name_, exp_name_, dat_name_, test_id_, relay_id_)
            dat_test_in_name = [dat_dir_name_ exp_name_ dat_name_ '.re' num2str(relay_id_) '_test' num2str(test_id_) '.redat'];
            dat_test_in = fopen(dat_test_in_name);
            header_in = fread(dat_test_in, [1,2], 'int=>int');
            [param_len,npool] = deal(header_in(1),header_in(2));
            params = fread(dat_test_in,[param_len,npool], 'double=>double');
            fclose(dat_test_in);

            test_directory = [dat_dir_name_ exp_name_ dat_name_  '.re' num2str(relay_id_) '_test' num2str(test_id_) '_results/'];

            dat_test_specs = fopen([test_directory, 'results_specs.redat']);
            results_specs = fread(dat_test_specs, [1,4], 'int=>int');
            [Frame_end, beads, dof, Frame_full] = deal(results_specs(1), results_specs(2), results_specs(3), results_specs(4));
            ref_pos_mat=fread(dat_test_specs, [beads*dof,Frame_end], 'double=>double');
            fclose(dat_test_specs);

            earlylate_mat = nan(2,npool);
            acc_res_mat = nan(2,npool);
            [sim_pos_mat, pos_mat] = deal(nan(Frame_end*beads*dof, npool));
            [pos_res_mat, alpha_mat, INTmat] = deal(nan(Frame_end*beads, npool));
            [pos_res_re_mat, alpha_re_mat] = deal(nan(Frame_full*beads,npool));
            [kill_frame_mat, alpha_kill_mat] = deal(nan(beads, npool));

            for i = 1:npool
                dat = fopen([test_directory 'par' num2str(i-1) '.redat']);
                header = fread(dat, [1,2], 'int=>int');
                earlylate_mat(:,i) = fread(dat, [header(2),1], 'int=>int');
                acc_res_mat(:,i) = fread(dat, [header(1),1], 'double=>double');
                sim_pos_mat(:,i) = fread(dat, [Frame_end*beads*dof,1], 'double=>double');
                pos_mat(:,i) = fread(dat, [Frame_end*beads*dof,1], 'double=>double');
                pos_res_mat(:,i) = fread(dat, [Frame_end*beads,1], 'double=>double');
                alpha_mat(:,i) = fread(dat, [Frame_end*beads,1], 'double=>double');
                INTmat(:,i) = fread(dat, [Frame_end*beads,1], 'double=>double');
                pos_res_re_mat(:,i) = fread(dat, [Frame_full*beads, 1], 'double=>double');
                alpha_re_mat(:,i) = fread(dat, [Frame_full*beads, 1], 'double=>double');
                kill_frame_mat(:,i) = fread(dat, [beads, 1], 'int=>int');
                alpha_kill_mat(:,i) = fread(dat, [beads, 1], 'double=>double');
                fclose(dat);
            end

            %% set properties
            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.relay_id = relay_id_;
            obj.test_id = test_id_;

            obj.param_len=param_len;
            obj.npool=npool;
            obj.params = params;

            obj.Frame_end=Frame_end;
            obj.beads=beads;
            obj.dof=dof;
            obj.Frame_full=Frame_full;
            obj.ref_pos_mat=ref_pos_mat;

            obj.earlylate_mat=earlylate_mat;
            obj.acc_res_mat=acc_res_mat;
            obj.sim_pos_mat=sim_pos_mat;
            obj.pos_mat=pos_mat;
            obj.pos_res_mat=pos_res_mat;
            obj.alpha_mat=alpha_mat;
            obj.INTmat=INTmat;
            obj.pos_res_re_mat=pos_res_re_mat;
            obj.alpha_re_mat=alpha_re_mat;
            obj.kill_frame_mat=kill_frame_mat;
            obj.alpha_kill_mat=alpha_kill_mat;

            obj.swtrue = sw_;
        end
        function cell_out = frame_posmat2poolcell(obj, mat_in)
            cell_out = cell([obj.npool, 1]);
            for i = 1:obj.npool
                cell_out{i} = reshape(mat_in(:,i), [obj.dof, obj.beads, obj.Frame_end]);
            end
        end
        function cell_out = frame_posmat2beadcell(obj, mat_in)
            cell_out = cell([obj.beads, 1]);
            [cell_out{:}] = deal(nan(obj.dof, obj.Frame_end, obj.npool));
            for i = 1:obj.npool
                tens_i = reshape(mat_in(:,i), [obj.dof, obj.beads, obj.Frame_end]);
                for j = 1:obj.beads
                    cell_out{j}(:,:,i) = reshape(tens_i(:,j,:), [obj.dof,obj.Frame_end]);
                end
            end
        end
        function cell_out = frame_resmat2beadcell(obj, mat_in)
            cell_out = cell([obj.beads, 1]);
            [cell_out{:}] = deal(nan(obj.Frame_end, obj.npool));
            mat_1 = reshape(mat_in, [obj.beads, obj.Frame_end, obj.npool]);
            for i = 1:obj.beads
                cell_out{i} = (reshape(mat_1(i,:,:), [obj.Frame_end, obj.npool]))';
            end
        end
        function beadres_out = frame_simbeadres(obj, pos_true_, bcell_in_)
            [beads, Frames, npool] = deal(obj.beads, obj.Frame_end, obj.npool);
            if (nargin==3)
                bcell_=bcell_in_;
            else
                bcell_=obj.frame_posmat2beadcell(obj.sim_pos_mat);
            end

            beadres_out = cell([beads,1]);
            [beadres_out{:}] = deal(nan(npool, Frames));

            for i=1:beads
                truei = reshape(pos_true_(i,:,:), [2,Frames]);
                for j = 1:npool
                    diffj = truei-bcell_{i}(:,:,j);
                    beadres_out{i}(j,:) = sum(diffj.*diffj, 1);
                end
            end
        end
        function beadres_out = beadres_matcomp(obj)
            [beads, Frames, npool] = deal(obj.beads, obj.Frame_end, obj.npool);
            beadres_out = cell([beads,1]);
            [beadres_out{:}] = deal(nan(npool, Frames));

            resmat=obj.ref_pos_mat(:)-obj.pos_mat;
            res = reshape(sum(reshape(resmat.*resmat, [2, beads, Frames, npool]),1),[beads,Frames,npool]);
            for i = 1:beads
                beadres_out{i} = (reshape(res(i,:,:),[Frames, npool]))';
            end
        end
        function [refsimpos_out, simpos_out] = datpos2simpos(obj, dish_in_, pcell_, i_)
            dish_ = dish_in_(:,2:3);
            [refsimpos_out, simpos_out] = deal(nan(obj.beads, 2, obj.Frame_end));
            refposi = permute(reshape(obj.ref_pos_mat, [2,obj.beads,obj.Frame_end]), [2 1 3]);
            posi = permute(pcell_{i_},[2 1 3]);
            [cxim,cyim,cl] = deal(obj.params(9,i_), obj.params(10,i_), obj.params(11,i_));
            for i = 1:obj.Frame_end
                refsimpos_out(:,:,i) = [(refposi(:,1,i)-cxim)/(cl)+dish_(i,1), (refposi(:,2,i)-cyim)/(cl)+dish_(i,2)];
                simpos_out(:,:,i) = [(posi(:,1,i)-cxim)/(cl)+dish_(i,1), (posi(:,2,i)-cyim)/(cl)+dish_(i,2)];
            end
        end
    end
    methods (Static)

    end
end
