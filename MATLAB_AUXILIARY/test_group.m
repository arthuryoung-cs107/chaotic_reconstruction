classdef test_group
    properties
        rad=0.5;
        mass=1.0;
        param_len=12;
        true_len=14;

        int_chunk_count=1;
        double_chunk_count=3;

        ref_pos;
        npool;

        params_mat;

        int_chunk;
        double_chunk;
        results;
        pos_res;

        alpha_bead;
    end
    methods
        function obj = test_group(nbeads_,par_id_,relay_id_,test_id_)
            name = ['../dat_dir/swirl' num2str(nbeads_) '.odr/pts.' num2str(par_id_) '.re' num2str(relay_id_) '_test' num2str(test_id_)];

            param_file = fopen([name  '.redat']);
            param_head = fread(param_file, [1, 2], 'int=>int');
            param_mat = fread(param_file, param_head, 'double=>double');
            fclose(param_file);

            dir_name = [name '_results/'];

            specs_file = fopen([dir_name 'results_specs.redat']);
            specs_head = fread(specs_file, [1 4], 'int=>int');
            [Frame_end, nbeads, dof, Frames] = deal(specs_head(1), specs_head(2), specs_head(3), specs_head(4));
            ref_pos_raw = fread(specs_file, [dof*nbeads,Frame_end], 'double=>double');
            fclose(specs_file);

            [param_len, npool] = deal(param_head(1), param_head(2));
            ref_pos = reshape(ref_pos_raw, [dof nbeads Frame_end]);

            int_chunk=nan(obj.int_chunk_count*nbeads, npool);
            double_chunk=nan(obj.double_chunk_count*nbeads, npool);

            results(npool) = struct('sim_pos', nan(dof*nbeads, Frame_end), ...
                                    'test_pos', nan(dof*nbeads, Frame_end), ...
                                    'test_posres', nan(nbeads, Frame_end), ...
                                    'test_alpha', nan(nbeads, Frame_end), ...
                                    'test_intposres', nan(nbeads, Frame_end), ...
                                    'pos_res', nan(nbeads, Frames), ...
                                    'alpha', nan(nbeads, Frames));

            pos_res = nan(npool,Frame_end);
            alpha_bead = nan(nbeads,Frame_end,npool);

            for i = 1:npool
                pool_file = fopen([dir_name 'par' num2str(i-1) '.redat']);
                pool_head = fread(pool_file, [1 2], 'int=>int');
                int_chunk(:,i) = fread(pool_file, pool_head(1), 'int=>int');
                double_chunk(:,i) = fread(pool_file, pool_head(2), 'double=>double');
                params = fread(pool_file, param_len, 'double=>double');

                results(i).sim_pos = fread(pool_file, [dof*nbeads,Frame_end], 'double=>double');
                results(i).test_pos = fread(pool_file, [dof*nbeads,Frame_end], 'double=>double');
                results(i).test_posres = fread(pool_file, [nbeads,Frame_end], 'double=>double');
                results(i).test_alpha = fread(pool_file, [nbeads,Frame_end], 'double=>double');
                results(i).test_intposres = fread(pool_file, [nbeads,Frame_end], 'double=>double');
                results(i).pos_res = fread(pool_file, [nbeads,Frame_end], 'double=>double');
                results(i).alpha = fread(pool_file, [nbeads,Frame_end], 'double=>double');

                pos_res(i,:) = sum(results(i).test_posres);

                alpha_bead(:,:,i) = results(i).test_alpha;

                fclose(pool_file);
            end

            obj.params_mat=param_mat; 

            obj.npool=npool;
            obj.ref_pos=ref_pos;

            obj.int_chunk = int_chunk;
            obj.double_chunk = double_chunk;
            obj.results = results;
            obj.pos_res = pos_res;
            obj.alpha_bead = permute(alpha_bead,[3 2 1]);
        end
    end
end
