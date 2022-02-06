clear
close all
dat_dir_name = '../dat_dir/';
circ6_name = 'circ6.odr/';
circ6_bintest_name = 'circ6_bintest.odr/';
dat_name = 'pts';
data_directory = circ6_name;

circ6 = ODR_data(dat_dir_name, circ6_name, dat_name);
circ6_bintest = ODR_data(dat_dir_name, circ6_bintest_name, dat_name);
