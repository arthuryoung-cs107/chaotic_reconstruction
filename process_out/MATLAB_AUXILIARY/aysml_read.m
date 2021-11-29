function mat_return = aysml_read(name)
  dims = dlmread([name '.aysml']);
  m = dims(1);
  n = dims(2);
  id = fopen([name '.aydat']);
  mat_return = (fread( id,[m, n], 'float64=>float64'));
  fclose(id);
end
