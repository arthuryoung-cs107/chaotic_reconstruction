function sw = read_swirl(nbeads, par_id)
  proc_name = '../dat_dir/';
  dat_subfix ='pts';
  exp_name = ['swirl' num2str(nbeads) '.odr/'];
  file_name = [dat_subfix '.' num2str(par_id)];
  sw = ODR_data.construct_swirl(proc_name, exp_name, file_name);
end
