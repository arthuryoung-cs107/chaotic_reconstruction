function det = read_detection_group(swtrue, nbeads, par_id, test_id, relay_id)
    proc_name = '../dat_dir/';
    dat_name = 'pts';
    exp_name = ['swirl' num2str(nbeads) '.odr/'];
    file_name = [dat_name '.' num2str(par_id)];
    det = detection_group(swtrue, proc_name, exp_name, file_name, test_id, relay_id);
end
