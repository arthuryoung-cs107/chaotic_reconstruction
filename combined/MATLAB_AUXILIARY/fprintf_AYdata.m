function fprintf_matrix(mat, name)
  if (length(size(mat)) == 2)
    file_id = fopen([name, '.aydat'], 'w+');
    fwrite(file_id, mat(:), 'double');
    fclose(file_id);
    file_id2 = fopen([name, '.aysml'], 'w+');
    fprintf(file_id2, '%d %d', size(mat, 1), size(mat, 2));
    fclose(file_id2);
  elseif length(size(mat) == 3)
    file_id = fopen([name, '.aytens'], 'w+');
    for i=1:size(mat, 3)
      fwrite(file_id, reshape(mat(:, :, i), [size(mat, 1)*size(mat, 2), 1]), 'double');
    end
    fclose(file_id);
    file_id2 = fopen([name, '.aysml'], 'w+');
    fprintf(file_id2, '1 %d %d %d', size(mat, 1), size(mat, 2), size(mat, 3));
    fclose(file_id2);
  end
end
