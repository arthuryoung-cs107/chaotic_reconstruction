clear
close all

sym1 = aysml_read('../dat_dir/sym1')
L = chol(sym1)

b = (1:size(sym1, 2))'
x = sym1\b
