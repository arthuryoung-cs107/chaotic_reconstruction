clear
close all

mat1 = aysml_read('../dat_dir/mat1')
fprintf_AYdata(mat1, '../dat_dir/mat1_back');

tens1 = aysml_read('../dat_dir/tens1')
fprintf_AYdata(tens1, '../dat_dir/tens1_back');

sym3 = aysml_read('../dat_dir/sym3')
