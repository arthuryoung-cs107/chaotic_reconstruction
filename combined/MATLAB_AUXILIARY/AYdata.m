classdef AYdata < handle
  properties
    len;
  end
  methods(Static)
    function obj = AYdata(len_)
      obj.len = len_;
    end
    function dat_return = aysml_read(name)
      dims = dlmread([name, '.aysml']);
      if (size(dims, 1) == 1)
        switch size(dims, 2)
          case 2 %% matrix or vector case
            m = dims(1);
            n = dims(2);
            id = fopen([name '.aydat']);
            dat_return = (fread( id,[m, n], 'float64=>float64'));
            fclose(id);

          case 3 %% AYsym case
            m = dims(2);
            n = m;
            dat_return = zeros(m, n);
            id = fopen([name '.aydat']);
            for i=1:n
              row = (fread( id,[n-i+1, 1], 'float64=>float64'));
              dat_return(i, i) = row(1);
              dat_return(i, i+1:n) = row(2:length(row));
              dat_return(i+1:n, i) = row(2:length(row));
            end
            fclose(id);
          case 4 %% tensor case
            m = dims(2);
            n = dims(3);
            w = dims(4);
            dat_return = nan(m, n, w);
            if (dims(1)== 1) %% one file to source tensor
              id = fopen([name '.aydat']);
              for i=1:w
                dat_return(:, :, i) = (fread( id,[m, n], 'float64=>float64'));
              end
              fclose(id);
            elseif (dims(1)== 0) %% a file for every page of tensor
              for i=1:w
                id = fopen([name '.' num2str(i-1) '.aydat']);
                dat_return(:, :, i) = (fread( id,[m, n], 'float64=>float64'));
                fclose(id);
              end
            end
        end
      else
        fprintf('aysml_read: Failed. Implement the data structure case');
        dat_return = 0;
      end
    end
    function aydat_write(mat, name)
      switch length(size(mat))
        case 2
          file_id = fopen([name, '.aydat'], 'w+');
          fwrite(file_id, mat(:), 'double');
          fclose(file_id);
          file_id2 = fopen([name, '.aysml'], 'w+');
          fprintf(file_id2, '%d %d', size(mat, 1), size(mat, 2));
          fclose(file_id2);
        case 1
          file_id = fopen([name, '.aydat'], 'w+');
          fwrite(file_id, mat(:), 'double');
          fclose(file_id);
          file_id2 = fopen([name, '.aysml'], 'w+');
          fprintf(file_id2, '%d %d', length(mat), 1);
          fclose(file_id2);
        case 3
          file_id = fopen([name, '.aydat'], 'w+');
          for i=1:size(mat, 3)
            fwrite(file_id, reshape(mat(:, :, i), [size(mat, 1)*size(mat, 2), 1]), 'double');
          end
          fclose(file_id);
          file_id2 = fopen([name, '.aysml'], 'w+');
          fprintf(file_id2, '1 %d %d %d', size(mat, 1), size(mat, 2), size(mat, 3));
          fclose(file_id2);
      end
    end
    function aydat_write_split_tens(tens, name)
      if (length(size(tens))==3)
        for i=1:size(tens, 3)
          file_id = fopen([name '.' num2str(i-1) '.aydat'], 'w+');
          fwrite(file_id, reshape(tens(:, :, i), [size(tens, 1)*size(tens, 2), 1]), 'double');
          fclose(file_id);
        end
        file_id2 = fopen([name, '.aysml'], 'w+');
        fprintf(file_id2, '0 %d %d %d', size(tens, 1), size(tens, 2), size(tens, 3));
        fclose(file_id2);
      else
        fprintf('AYdata: aydat_write_split_tens failed, did not pass in a tensor');
      end
    end
  end
end
