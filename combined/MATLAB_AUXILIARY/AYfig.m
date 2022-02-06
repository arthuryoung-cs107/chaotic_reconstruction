classdef AYfig < handle
  properties
    fig;
    props_in_;
  end
  methods(Static)
    function obj = fluid(props_in_)
      obj.fig = figure;
      obj.props_in = props_in_;
      for i=1:size(props_in_, 1)
        obj.fig.set(props_in_{i, 1}, props_in_{i, 2});
      end
    end
    function fig_out = figure(props_in_)
      fig_out = figure;
      for i=1:size(props_in_, 1)
        fig_out.set(props_in_{i, 1}, props_in_{i, 2});
      end
    end
    function struct_out = specs_gen(name_in_, pos_in_)
      struct_out = {'Name', name_in_; 'Renderer', 'painters'; 'Position', pos_in_;};
    end
  end
end
