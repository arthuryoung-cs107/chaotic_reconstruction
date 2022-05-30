function relay = read_relay(nbeads, par_id, relay_id, ended_)
    if (nargin==3)
        ended=false;
    else
        ended = ended_;
    end

    proc_name = '../dat_dir/';
    dat_name = 'pts';
    exp_name = ['swirl' num2str(nbeads) '.odr/'];
    file_name = [dat_name '.' num2str(par_id)];
    relay = relay_data(proc_name, exp_name, file_name, relay_id);

    if (ended)
        [relay.gen_ind, relay.gen, relay.gen_end_data] = relay.read_gen_end; 
    end
end
