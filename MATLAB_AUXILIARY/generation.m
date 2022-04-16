classdef generation
    properties
        dat_dir_name;
        exp_name;
        dat_name;
        gen_count;

        race_id;

        specs_len = 7;
        specs;
        gen_correct;
        leader_count;
        par_len;
        wleader_index;
        bleader_index;

        sample_weights;

        records;

        frscores;
        l2scores;
        dup_vec;

        zvec;
        lambda_z;
        w;
        param_mat;
    end
    methods
        function obj = generation(dat_dir_name_, exp_name_, dat_name_, gen_count_, race_id_)
            obj.dat_dir_name = dat_dir_name_;
            obj.exp_name = exp_name_;
            obj.dat_name = dat_name_;
            obj.gen_count = gen_count_;
            obj.race_id = race_id_;

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
        function plot_diagnostics(obj, AYfig_, F_, lambda_)
            green4 = [0.1000 0.4440 0.2440];
            blue5 = [0 0 1];
            red5 = [1 0 0];
            orange1 = [1.0 0.3529 0];

            ax_ = AYfig_.ax_tile;
            frscores = obj.frscores;
            zvec = obj.zvec;
            lambda_z = obj.lambda_z;
            leader_count = obj.leader_count;
            dup_vec = obj.dup_vec;
            w = obj.w;
            weights = obj.sample_weights;

            [fr_sort, I_fr] = mink(frscores,leader_count);
            [z_sort, I_z] = mink(zvec, leader_count);

            dup_zsort = dup_vec(I_z);
            w_zsort = w(I_z);
            weights_zsort = weights(I_z);

            xlabel(ax_(1), 'frame score', 'Interpreter', 'Latex', 'Fontsize', 14)
                yyaxis(ax_(1), 'left');
                histogram(ax_(1), frscores);
                ylabel(ax_(1), 'frequency', 'Interpreter', 'Latex', 'Fontsize', 14)

                yyaxis(ax_(1), 'right');
                plot(ax_(1),frscores, zvec, ' o', 'Color', red5, 'LineWidth', 1);
                ylabel(ax_(1), 'z', 'Interpreter', 'Latex', 'Fontsize', 14)

            xlabel(ax_(2), 'z', 'Interpreter', 'Latex', 'Fontsize', 14)
                yyaxis(ax_(2), 'left');
                histogram(ax_(2), zvec);
                ylabel(ax_(2), 'frequency', 'Interpreter', 'Latex', 'Fontsize', 14)

                yyaxis(ax_(2), 'right');
                % fplot(ax_(2), @(z) lambda_z*exp(-lambda_z.*z), [min(zvec),max(zvec)], 'Color', red5, 'LineWidth', 5);
                plot(ax_(2),zvec, w/sum(w), ' o', 'Color', red5, 'LineWidth', 1);
                hold(ax_(2), 'on')

                plot(ax_(2),zvec, weights, ' o', 'Color', orange1, 'LineWidth', 1);
                ylabel(ax_(2), '$$\lambda_z \mathrm{exp}(-\lambda_z z)$$', 'Interpreter', 'Latex', 'Fontsize', 14)
                hold(ax_(2), 'off')

            xlabel(ax_(3), 'z', 'Interpreter', 'Latex', 'Fontsize', 14)
                % yyaxis(ax_(3), 'left');
                % histogram(ax_(3), 'BinEdges', linspace(0.5, length(dup_vec)+0.5, length(dup_vec)+1), 'BinCounts', dup_zsort);
                % plot(ax_(3),I_z,dup_zsort,' o', 'Color', blue5, 'LineWidth', 1);
                % ylabel(ax_(3), 'frequency', 'Interpreter', 'Latex', 'Fontsize', 14)

                yyaxis(ax_(3), 'right');
                % plot(ax_(3),1:leader_count, weights_zsort, ' o', 'Color', orange1, 'LineWidth', 1);
                plot(ax_(3),1:leader_count, weights, ' o', 'Color', orange1, 'LineWidth', 1);
                ylabel(ax_(3), 'weights', 'Interpreter', 'Latex', 'Fontsize', 14)
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
