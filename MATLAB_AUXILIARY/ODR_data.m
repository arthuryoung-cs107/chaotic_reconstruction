classdef ODR_data < swirl
    properties
        dat_dir_name;
        exp_name;
        dat_name;

        filin;

    end
    methods (Static)
        function sw = construct_swirl(dat_dir_name_, exp_name_, dat_name_)
          rysml = ODR_data.get_rysml(dat_dir_name_, exp_name_, dat_name_);
          rydat = ODR_data.get_rydat(rysml);
          sw = swirl(rydat.pos, rydat.dish, rysml.Frames, rysml.len_dish, rysml.len_pos, rysml.beads);
        end
        function rysml_out = get_rysml(dat_dir_name_,exp_name_,dat_name_)
            rysml = dlmread([dat_dir_name_ exp_name_ dat_name_ '.rysml']);
            rysml_out = struct('dat_dir', dat_dir_name_, 'exp_name',exp_name_,'dat_name',dat_name_,'split',rysml(1,1),'Frames',rysml(1,2),'beads',rysml(3,2),'len_dish',rysml(2, 3),'len_pos',rysml(3,3));
        end
        function rydat_out = get_rydat(rysml)
            dish = nan(rysml.Frames, rysml.len_dish);
            pos = nan(rysml.beads, rysml.len_pos, rysml.Frames);
            if (rysml.split==1) %% if the data is stored across many binary files
                for f=0:(rysml.Frames-1)
                    id = fopen([rysml.dat_dir rysml.exp_name rysml.dat_name '.' num2str(f) '.aydat']);
                    dish(f+1, :) = (fread( id,[1, rysml.len_dish], 'float64=>float64'));
                    pos(:, :, f+1) = (fread( id,[rysml.len_pos, rysml.beads], 'float64=>float64'))';
                    fclose(id);
                end
            elseif (rysml.split==0) %% if the data is consolidated into one binary file
                id = fopen([rysml.dat_dir rysml.exp_name rysml.dat_name '.rydat']);
                for f=0:(rysml.Frames-1)
                    dish(f+1, :) = fread( id,[1, rysml.len_dish], 'float64=>float64');
                    pos(:, :, f+1) = (fread( id,[rysml.len_pos, rysml.beads], 'float64=>float64'))';
                end
                fclose(id);
            end
            rydat_out = struct('pos', pos, 'dish', dish);
        end
    end
    methods
        function obj = ODR_data(dat_dir_name_, exp_name_, dat_name_)
            %% first, construct the swirl that represents this data simulation outputs
            rysml = ODR_data.get_rysml(dat_dir_name_, exp_name_, dat_name_);
            rydat = ODR_data.get_rydat(rysml);

            obj@swirl(rydat.pos, rydat.dish, rysml.Frames, rysml.len_dish, rysml.len_pos, rysml.beads);

            %% ODR member assignments
            obj.dat_dir_name = rysml.dat_dir;
            obj.exp_name = rysml.exp_name;
            obj.dat_name = rysml.dat_name;
        end
        function sw = spawn_swirl(obj)
            sw = swirl(obj.pos, obj.dish, obj.Frames, obj.len_dish, obj.len_pos, obj.beads);
        end
        function load_filin(obj)
            obj.filin = filter_inputs([obj.dat_dir_name obj.exp_name obj.dat_name]);
        end
        function gen_out = read_geni(obj, gen_count_)
            gen_out = generation(obj.dat_dir_name, obj.exp_name, obj.dat_name, gen_count_);
        end
        function gen_out = read_race_data(obj)
            gen_out = race_data(obj.dat_dir_name, obj.exp_name, obj.dat_name);
        end
        function make_POVray_inputs(obj, pov_dir_, dat_name_)
            save_dir = [pov_dir_ obj.exp_name];
            status = mkdir(save_dir);
            for i=1:obj.Frames
                id = fopen([save_dir dat_name_ '.' num2str(i-1)], 'w');
                fprintf(id, '%g %g %g %g\n', obj.dish(i, :));
                fprintf(id, '%g %g %g %g %g %g %g\n', obj.pos(:, :, i)');
                fclose(id);
            end
        end
        function sw_out = spawn_best(obj)
            sw_out = ODR_data([obj.dat_dir_name obj.exp_name], 'swirl_best.odr/', obj.dat_name);
        end
    end
end
