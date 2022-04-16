function odr = read_race(nbeads, par_id, race_id)
  proc_name = '../dat_dir/';
  dat_name = 'pts';
  exp_name = ['swirl' num2str(nbeads) '.odr/'];
  file_name = [dat_name '.' num2str(par_id)];
  odr = race_data(proc_name, exp_name, file_name, race_id);
end
