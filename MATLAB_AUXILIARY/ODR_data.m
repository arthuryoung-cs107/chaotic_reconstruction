classdef ODR_data < swirl
  properties
    dat_dir_name;
    exp_name;
    dat_name;

    rysml; 
  end

  methods
    function obj = ODR_data(dat_dir_name_, exp_name_, dat_name_)
      obj@swirl();
      obj.dat_dir_name = dat_dir_name_;
      obj.exp_name = exp_name_;
      obj.dat_name = dat_name_;
      obj.rysml = dlmread([obj.dat_dir_name obj.exp_name obj.dat_name '.rysml']);
      obj.Frames = obj.rysml(1, 2);
      obj.len_specs = obj.rysml(2, 3);
      obj.len_dat = obj.rysml(3, 3);
      obj.P = obj.rysml(3, 2);

      obj.specs = nan(obj.Frames, obj.len_specs);
      obj.data = nan(obj.P, obj.len_dat, obj.Frames);
      if (obj.rysml(1,1)==1)
        for f=0:(obj.Frames-1)
          id = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.' num2str(f) '.aydat']);
          obj.specs(f+1, :) = (fread( id,[1, obj.len_specs], 'float64=>float64'));
          obj.data(:, :, f+1) = (fread( id,[obj.len_dat, obj.P], 'float64=>float64'))';
          fclose(id);
        end
      elseif (obj.rysml(1,1)==0)
        id = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.rydat']);
        for f=0:(obj.Frames-1)
          obj.specs(f+1, :) = (fread( id,[1, obj.len_specs], 'float64=>float64'));
          obj.data(:, :, f+1) = (fread( id,[obj.len_dat, obj.P], 'float64=>float64'))';
        end
        fclose(id);
      end
    end
    function load_filin(obj)
      obj.filin = filter_inputs([obj.dat_dir_name obj.exp_name obj.dat_name]);
    end
    function gen_out = read_geni(obj, gen_count_)
      gen_out = generation(obj.dat_dir_name, obj.exp_name, obj.dat_name, gen_count_);
    end
    function gen_out = read_race_data(obj)
      gen_out = race_data(obj.dat_dir_name, obj.exp_name, obj.dat_name);
    end
    function make_POVray_inputs(obj, pov_dir_, dat_name_)
      save_dir = [pov_dir_ obj.exp_name];
      status = mkdir(save_dir);
      for i=1:obj.Frames
        id = fopen([save_dir dat_name_ '.' num2str(i-1)], 'w');
        fprintf(id, '%g %g %g %g\n', obj.specs(i, :));
        fprintf(id, '%g %g %g %g %g %g %g\n', obj.data(:, :, i)');
        fclose(id);
      end
    end
    function sw_out = spawn_best(obj)
      sw_out = ODR_data([obj.dat_dir_name obj.exp_name], 'swirl_best.odr/', obj.dat_name);
    end
  end

end
