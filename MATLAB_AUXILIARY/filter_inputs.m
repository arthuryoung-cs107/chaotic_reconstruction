classdef filter_inputs < handle
  properties
    fisml_name;
    fisml;

    t;
    pos;
    theta;
    vidspecs;

    t_phys;
    cx_im;
    cy_im;
    cl_im;
  end

  methods
    function obj = filter_inputs(fisml_name_)
      obj.fisml_name = [fisml_name_, '.fisml'];
      obj.fisml = dlmread(obj.fisml_name);

      id = fopen([fisml_name_, '.filin']);
      obj.t = (fread( id,[1, obj.fisml(2, 3)], 'float64=>float64'));
      obj.pos = nan(obj.fisml(3, 2), obj.fisml(3, 3), obj.fisml(3, 1));
      for i=1:obj.fisml(3, 1)
        obj.pos(:, :, i) = (fread( id,[obj.fisml(3, 3), obj.fisml(3, 2)], 'float64=>float64'))';
      end
      obj.theta = (fread( id,[1, obj.fisml(4, 3)], 'float64=>float64'));
      obj.vidspecs = (fread( id,[1, obj.fisml(5, 3)], 'float64=>float64'));
      fclose(id);

      [obj.t_phys, obj.cx_im, obj.cy_im, obj.cl_im] = deal(obj.vidspecs(1), obj.vidspecs(2), obj.vidspecs(3), obj.vidspecs(4));


    end
  end

end
