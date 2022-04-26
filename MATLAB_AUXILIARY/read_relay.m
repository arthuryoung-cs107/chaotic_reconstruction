function relay = read_relay(nbeads, par_id, relay_id)
    proc_name = '../dat_dir/';
    dat_name = 'pts';
    exp_name = ['swirl' num2str(nbeads) '.odr/'];
    file_name = [dat_name '.' num2str(par_id)];
    relay = relay_data(proc_name, exp_name, file_name, relay_id);
end
