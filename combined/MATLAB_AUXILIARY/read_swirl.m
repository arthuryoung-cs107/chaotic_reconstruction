function sw = read_swirl(nbeads, par_id)
  proc_name = '../dat_dir/';
  pov_dir = '../POV_AUXILIARY/';
  dat_name = 'pts';
  exp_name = ['swirl' num2str(nbeads) '.odr/'];
  file_name = [dat_name '.' num2str(par_id)];
  sw = ODR_data(proc_name, exp_name, file_name);
end
