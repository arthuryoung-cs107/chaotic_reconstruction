classdef relay_generation
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        gen_count;
        relay_id;


    end
    methods
        function obj = relay_generation(dat_dir_name_, exp_name_, dat_name_, gen_count_, relay_id_)

            dat = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.rc' num2str(obj.race_id) '_gen' num2str(obj.gen_count) '.rcdat']);
            obj.specs = fread( dat,[1, obj.specs_len], 'int=>int');
            obj.gen_correct = obj.gen_count == obj.specs(1);
            obj.leader_count = obj.specs(2);
            obj.par_len = obj.specs(5);
            obj.wleader_index = obj.specs(6);
            obj.bleader_index = obj.specs(7);

            obj.sample_weights = fread( dat,[1, obj.leader_count], 'double=>double');

            obj.records = record.empty(obj.leader_count, 0);
            obj.frscores = nan(obj.leader_count, 1);
            obj.l2scores = nan(obj.leader_count, 1);
            obj.dup_vec = nan(obj.leader_count, 1);
            obj.zvec = nan(obj.leader_count, 1);
            obj.param_mat = nan(obj.par_len,obj.leader_count);
            for i=1:obj.leader_count
                int_vec = fread( dat,[1, obj.specs(3)], 'int=>int');
                double_vec = fread( dat,[1, obj.specs(4)], 'double=>double');
                param_it = fread( dat,[1, obj.par_len], 'double=>double');
                obj.records(i) = record(int_vec, double_vec, param_it);
                obj.frscores(i) = obj.records(i).frscore;
                obj.l2scores(i) = obj.records(i).l2score;
                obj.dup_vec(i) = obj.records(i).dup_count;
                obj.param_mat(:,i) = param_it';
            end
            fclose(dat);

            frmin = min(obj.frscores);
            obj.zvec = (1./(obj.frscores-frmin+1)).*(1./obj.frscores).*obj.l2scores;
            obj.lambda_z = 1/(mean(obj.zvec));
            obj.w = obj.lambda_z*exp(-obj.lambda_z*obj.zvec);

        end
        function best = get_best_record(obj)
            best = obj.records(obj.bleader_index+1);
        end
    end
    methods (Static)
        function plot_gen_scores(ax, frscores, fcount)
            fworst = min(frscores);
            fbest = max(frscores);
            frame_bins = fworst:fbest;
            P = fcount/(sum(fcount));
            lambda = dot(frame_bins,P);

            [fmode_count, i_mode] = max(fcount);
            fmode = fworst+i_mode-1;
            fcount_h = fcount(i_mode:end);
            frame_bins_h = fmode:fbest;
            P_h = fcount_h/(sum(fcount_h));
            lambda_h = dot(frame_bins_h,P_h);

            yyaxis(ax, 'left');
            histogram(ax, frscores, fbest-fworst+1);
            yyaxis(ax, 'right');
            plot(ax, frame_bins, poisspdf(frame_bins, lambda),' o -', 'Color', [0 0 1], 'LineWidth', 2);
            hold(ax, 'on')
            plot(ax, frame_bins, poisspdf(frame_bins, lambda_h),' o -', 'Color', [1 0 0], 'LineWidth', 2);
            hold(ax, 'off')
        end
    end
end
