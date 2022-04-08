classdef race_data < handle
  properties
    dat_dir_name;
    exp_name;
    dat_name;
    gen_count;

    F;
    lambda;

    specs;
    specs_len = 7;

    gen;
  end
  methods
    function obj = race_data(dat_dir_name_, exp_name_, dat_name_)
        [obj.dat_dir_name obj.exp_name obj.dat_name] = deal(dat_dir_name_,exp_name_,dat_name_);
        dat = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.end.rcdat']);
        obj.specs = fread( dat,[1, obj.specs_len], 'int=>int');
        fclose(dat);
        obj.gen_count = obj.specs(1);
        obj.gen = generation.empty(obj.gen_count, 0);
        for i=1:obj.gen_count
            obj.gen(i) = generation(obj.dat_dir_name, obj.exp_name, obj.dat_name, i);
        end
    end
    function plot_diagnosticsi(obj, AYfig_, i_)
        obj.gen(i_).plot_diagnostics(AYfig_, obj.F, obj.lambda);
    end
    function sw_out = spawn_swbest(obj)
        sw_out = ODR_data.construct_swirl([obj.dat_dir_name obj.exp_name], 'swirl_best.odr/', obj.dat_name);
        sw_out.params = AYdata.aysml_read([obj.dat_dir_name obj.exp_name 'swirl_best.odr/' obj.dat_name '.sparam_best']);
        sw_out.len_par = length(sw_out.params);
    end
  end
end
