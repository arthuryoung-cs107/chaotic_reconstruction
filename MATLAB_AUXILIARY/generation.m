classdef generation < handle
  properties
    dat_dir_name;
    exp_name;
    dat_name;
    gen_count;

    specs_len = 7;
    specs;
    gen_correct;
    leader_count;
    par_len;
    wleader_index;
    bleader_index;

    records;

    frscores;
    param_mat;
  end
  methods
    function obj = generation(dat_dir_name_, exp_name_, dat_name_, gen_count_)
      obj.dat_dir_name = dat_dir_name_;
      obj.exp_name = exp_name_;
      obj.dat_name = dat_name_;
      obj.gen_count = gen_count_;

      dat = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.gen' num2str(obj.gen_count) '.rcdat']);
      obj.specs = fread( dat,[1, obj.specs_len], 'int=>int');
      obj.gen_correct = obj.gen_count == obj.specs(1);
      obj.leader_count = obj.specs(2);
      obj.par_len = obj.specs(5);
      obj.wleader_index = obj.specs(6);
      obj.bleader_index = obj.specs(7);

      obj.records = record.empty(obj.leader_count, 0);
      obj.frscores = nan(obj.leader_count, 1);
      obj.param_mat = nan(obj.leader_count, obj.par_len);
      for i=1:obj.leader_count
        int_vec = fread( dat,[1, obj.specs(3)], 'int=>int');
        double_vec = fread( dat,[1, obj.specs(4)], 'double=>double');
        param_it = fread( dat,[1, obj.par_len], 'double=>double');
        obj.records(i) = record(int_vec, double_vec, param_it);
        obj.frscores(i) = obj.records(i).frscore;
        obj.param_mat(i, :) = param_it;
      end
      fclose(dat);
    end
    function plot_diagnostics(obj, AYfig_, F_, lambda_)
      fbest = max(obj.frscores);
      fworst = min(obj.frscores);
      fmu = mean(obj.frscores);
      frame_bins = fworst:fbest;

      f = mink(obj.frscores, obj.leader_count);
      w = exp(lambda_*obj.frscores/F_); %% raw frame score weights
      w = w/sum(w);

      fcount = histcounts(obj.frscores, max(obj.frscores)-min(obj.frscores)+1);
      lambda = dot(fcount, fworst:fbest)/(length(obj.frscores)) - fworst;

      w_old = @(f) exp(lambda_*f/F_);
      w_new = @(f) (poisspdf(f-fworst, lambda)).^-1;

      pi_old = double(fcount).*w_old(frame_bins)/sum(w_old(frame_bins));
      pi_new = double(fcount).*w_new(frame_bins)/sum(w_new(frame_bins));

      yyaxis(AYfig_.ax_tile(1), 'left');
      histogram(AYfig_.ax_tile(1), obj.frscores, max(obj.frscores)-min(obj.frscores)+1);
      yyaxis(AYfig_.ax_tile(1), 'right');
      plot(AYfig_.ax_tile(1), frame_bins, poisspdf(frame_bins-fworst, lambda),' o -', 'Color', [1 0 0], 'LineWidth', 2);

      plot(AYfig_.ax_tile(2), frame_bins, w_old(frame_bins)/sum(w_old(frame_bins)), ' o -', 'Color', [0 0 1], 'LineWidth', 2);
      hold(AYfig_.ax_tile(2), 'on');
      plot(AYfig_.ax_tile(2), frame_bins, w_new(frame_bins)/sum(w_new(frame_bins)), ' o -', 'Color', [1 0 0], 'LineWidth', 2);
      hold(AYfig_.ax_tile(2), 'off');
      set(AYfig_.ax_tile(2), 'YScale', 'log')

      plot(AYfig_.ax_tile(3), frame_bins, pi_new, ' o -', 'Color', [1 0 0], 'LineWidth', 2);
      hold(AYfig_.ax_tile(3), 'on');
      % plot(AYfig_.ax_tile(3), frame_bins, pi_old, ' o -', 'Color', [0 0 1], 'LineWidth', 2);
      hold(AYfig_.ax_tile(3), 'off');


    end
  end
end
