function dat_return = aysml_read(name)
  dims = dlmread([name '.aysml']);
  if (size(dims, 1) == 1)
    if (size(dims, 2) == 2)
      m = dims(1);
      n = dims(2);
      id = fopen([name '.aydat']);
      dat_return = (fread( id,[m, n], 'float64=>float64'));
      fclose(id);
    elseif ((size(dims, 2) == 4)&&(dims(1)== 1))      
      m = dims(2);
      n = dims(3);
      w = dims(4);
      dat_return = nan(m, n, w);
      id = fopen([name '.aytens']);
      for i=1:w
        dat_return(:, :, i) = (fread( id,[m, n], 'float64=>float64'));
      end
      fclose(id);
    end
  end
end
