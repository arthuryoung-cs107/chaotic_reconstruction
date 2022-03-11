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

      for i=1:obj.leader_count
        int_vec = fread( dat,[1, obj.specs(3)], 'int=>int');
        double_vec = fread( dat,[1, obj.specs(4)], 'double=>double');
        param_it = fread( dat,[1, obj.par_len], 'double=>double');
        obj.records(i) = record(int_vec, double_vec, param_it);
      end
      fclose(dat);
    end
  end
end
