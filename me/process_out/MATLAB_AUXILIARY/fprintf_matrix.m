function fprintf_matrix(mat, name)
  file_id = fopen([name, '.aydat'], 'w+');
  fwrite(file_id, mat(:), 'double');
  fclose(file_id);
  file_id2 = fopen([name, '.aysml'], 'w+');
  fprintf(file_id2, '%d %d', size(mat, 1), size(mat, 2));
  fclose(file_id2);
end
