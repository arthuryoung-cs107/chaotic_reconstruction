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
      s = svd(obj.param_mat(:, 1:8)); %% hard coded for now, these are the varying parameters

      fbest = max(obj.frscores);
      fworst = min(obj.frscores);
      fmu = mean(obj.frscores);

      [fhat, ij] = mink(obj.frscores-fworst, obj.leader_count);
      lambda_pi = 1/(mean(fhat));
      pi_fhat = lambda_pi*exp(-lambda_pi*fhat);
      lambda_poi = mean(fhat);
      pi_fhat_poi = poisspdf(fhat, lambda_poi);


      [f, ik] = mink(obj.frscores, obj.leader_count);
      w = exp(lambda_*obj.frscores/F_); %% raw frame score weights
      w = w/sum(w);
      frame_bins = fworst:fbest;
      w2 = exp(lambda_pi*(obj.frscores-fmu));
      w2 = w2/sum(w2);

      w_acc = nan(size(frame_bins));
      wnew_acc = nan(size(frame_bins));
      for i=1:length(frame_bins)
        w_acc(i) = sum(w(obj.frscores == frame_bins(i)));
        wnew_acc(i) = sum(w2(obj.frscores == frame_bins(i)));
      end
      w_acc = w_acc/sum(w_acc);
      wnew_acc = wnew_acc/sum(wnew_acc);

      w_old = @(f) (lambda_*exp(lambda_*f/F_))/(F_*(exp(lambda_/F_*fbest)-exp(lambda_/F_*fworst)));
      w_new = @(f) (lambda_pi*exp(lambda_pi*(f-fmu)))/(exp(lambda_pi*(fbest-fmu))-exp(lambda_pi*(fworst-fmu)));

      yyaxis(AYfig_.ax_tile(1), 'left');
      histogram(AYfig_.ax_tile(1), obj.frscores, max(obj.frscores)-min(obj.frscores)+1);
      yyaxis(AYfig_.ax_tile(1), 'right');
      % plot(AYfig_.ax_tile(1), obj.frscores(ij), pi_fhat,' o -', 'Color', [1 0 0], 'LineWidth', 2);
      plot(AYfig_.ax_tile(1), obj.frscores(ij), pi_fhat_poi,' o -', 'Color', [1 0 0], 'LineWidth', 2);

      plot(AYfig_.ax_tile(2), 1:8, s/s(1), ' o -', 'Color', [0 0 0], 'LineWidth', 2);
      set(AYfig_.ax_tile(2), 'YScale', 'log')

      yyaxis(AYfig_.ax_tile(3), 'left');
      fplot(@(f) w_old(f), [fworst, fbest], '-', 'Color', [1 0 0], 'LineWidth', 2);
      hold(AYfig_.ax_tile(3), 'on');
      fplot(@(f) w_new(f), [fworst, fbest], ':', 'Color', [1 0 0], 'LineWidth', 2);
      hold(AYfig_.ax_tile(3), 'off');

      yyaxis(AYfig_.ax_tile(3), 'right');
      % plot(AYfig_.ax_tile(3), frame_bins, w_acc,' o -', 'Color', [0 0 1], 'LineWidth', 2);
      hold(AYfig_.ax_tile(3), 'on');
      plot(AYfig_.ax_tile(3), frame_bins, wnew_acc,'. :', 'Color', [0 0 1], 'LineWidth', 2);
      hold(AYfig_.ax_tile(3), 'off');

    end
  end
end
