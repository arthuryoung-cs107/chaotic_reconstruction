clear
close all
dat_dir_name = '../dat_dir/';
circ6_name = 'circ6.odr/';
circ6_swrl_name = 'circ6_swrl.odr/';
dat_name = 'pts';

circ6 = ODR_data(dat_dir_name, circ6_name, dat_name);
circ6_swrl = ODR_data(dat_dir_name, circ6_swrl_name, dat_name);
circ6_swrl.load_filin();
