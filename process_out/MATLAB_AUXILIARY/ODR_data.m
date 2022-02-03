classdef ODR_data < handle
  properties
    dat_dir_name;
    exp_name;
    dat_name;

    Frames;
    len_specs;
    len_dat;
    P;

    rysml;
    specs;
    data;
  end

  methods
    function obj = ODR_data(dat_dir_name_, exp_name_, dat_name_)
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

      for f=0:(obj.Frames-1)
        id = fopen([obj.dat_dir_name obj.exp_name obj.dat_name '.' num2str(f) '.aydat']);
        obj.specs(f+1, :) = (fread( id,[1, obj.len_specs], 'float64=>float64'));
        obj.data(:, :, f+1) = (fread( id,[obj.len_dat, obj.P], 'float64=>float64'))';
        fclose(id);
      end
    end
  end

end
