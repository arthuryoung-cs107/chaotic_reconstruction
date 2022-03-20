classdef race_data < handle
  properties
    dat_dir_name;
    exp_name;
    dat_name;
    gen_count;

    specs;
    specs_len = 7;

    gen;
  end
  methods
    function obj = race_data(dat_dir_name_, exp_name_, dat_name_)
      obj.dat_dir_name = dat_dir_name_;
      obj.exp_name = exp_name_;
      obj.dat_name = dat_name_;

      dat = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.end.rcdat']);
      obj.specs = fread( dat,[1, obj.specs_len], 'int=>int');
      fclose(dat);
      obj.gen_count = obj.specs(1);
      obj.gen = generation.empty(obj.gen_count, 0);
      for i=1:obj.gen_count
        obj.gen(i) = generation(obj.dat_dir_name, obj.exp_name, obj.dat_name, i);
      end
    end
    function frscore_histogrami(obj, AYfig_in, i_)
      obj.gen(i_).frscore_histogram(AYfig_in);
    end
    function sigmas_ploti(obj, AYfig_in, i_)
      obj.gen(i_).sigmas_plot(AYfig_in);
    end
  end
end
