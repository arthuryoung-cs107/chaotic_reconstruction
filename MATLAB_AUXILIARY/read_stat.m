function stat = read_stat(nbeads, test, ran_id)
  dat_dir_name = '../dat_dir/';
  dat_name = 'rand';
  exp_name = ['stat', num2str(nbeads), '_' test num2str(ran_id) '.odr/'];
  stat = stat_data(dat_dir_name, exp_name, dat_name, 0);
end
