classdef swirl
  properties

    %% data shared with swirl ODR data
    Frames=0; % number of frames in swirl
    beads=0; % number of beads
    len_pos=0; % row length of position matrix; (position, quaternion)
    len_dish=0; % length of dish data in each frame
    len_par=14;
    pos=0; % Note that a swirl IS this tensor
    dish=0; % (t cx cy wall_sca) for each frame
    params; % the parameters associated with the swirl

    contact_f;
    bead_contact_f;
    wall_contact_f;
    
  end
  methods
    function data = swirl(pos_, dish_, Frames_, len_dish_, len_pos_, beads_)
        if(nargin==1)
            data=pos_;
        elseif(nargin<6)
            data.pos=pos_;
            data.beads = size(pos_, 1);
            data.len_pos = size(pos_, 2);
            data.Frames = size(pos_, 3);
            if (nargin==2)
                data.len_dish = size(dish_, 2);
                data.dish=dish_;
            end
        elseif (nargin==6)
            data.pos = pos_(:,:,1:Frames_);
            data.dish = dish_(1:Frames_, :);
            [data.Frames, data.len_dish, data.len_pos, data.beads] = deal(Frames_, len_dish_, len_pos_, beads_);
        end
    end
    function t_vec_out = t_vec(obj)
        t_vec_out = obj.dish(:,1);
    end
    function [out, valid_sub_swirl] = spawn_sub_swirl(obj, oth)
        %% on output, we want to give the sub swirl of the input that matches THIS swirl
        out = 0;
        valid_sub_swirl = false;
        [f,not_lost,fcap,tol] = deal(0, true, min(obj.Frames, oth.Frames), 1);
        diff_full = obj(:, 1:2, :)-oth(:, 1:2, :);
        while ((f<fcap)&&(not_lost))
            f=f+1;
            diff = diff_full(:, :, f);
            not_lost = max(sqrt(sum((diff.*diff)')))<tol;
        end
        f_final = f;
        if (f_final > 1)
            valid_sub_swirl = true;
            out = swirl(oth, oth.dish, f_final, oth.len_dish, oth.len_pos, oth.beads);
        end
    end
    function wall_contact_frames = find_wall_contact_frames(obj)
        [Frames, dish, pos] = deal(obj.Frames, obj.dish, obj.pos);
        wall_contact_frames = 1:Frames;
        contact_detected = zeros(Frames, 1);
        fa = sqrt(0.75);
        r = 0.5;
        d = 5.72;
        walls = [0, 1; fa, 0.5; fa , -0.5]';
        for i=1:Frames
            s_mat = (pos(:, 1:2, i)-dish(i, 2:3))*walls;
            w_mat = s_mat - d*(double(s_mat>0)) + d*(double(s_mat<0));
            contact_detected(i) = sum(abs(w_mat(:))<r)>0;
        end
        wall_contact_frames = wall_contact_frames(logical(contact_detected));
    end
    function [contact_frames, bead_contact_frames,wall_contact_frames] = find_contact_frames(obj)
        [Frames, dish, pos, beads] = deal(obj.Frames, obj.dish, obj.pos, obj.beads);
        full_frames = 1:Frames;
        wall_contact = zeros(Frames, 1);
        bead_contact = zeros(Frames, 1);
        fa = sqrt(0.75);
        r = 0.5;
        tol = 2*r;
        d = 5.72;
        walls = [0, 1; fa, 0.5; fa , -0.5]';
        for i=1:Frames
            posi = pos(:, 1:2, i);
            s_mat = (posi-dish(i, 2:3))*walls;
            w_mat = s_mat - d*(double(s_mat>0)) + d*(double(s_mat<0));
            wall_contact(i) = sum(abs(w_mat(:))<r)>0;

            bead_contact_count=0;
            for j = 1:beads
                beadj = posi(j,:);
                for k = j+1:beads
                    if (norm(beadj-posi(k,:)) <= tol)
                        bead_contact_count=1;
                        break;
                    end
                end
                if (bead_contact_count>0)
                    break;
                end
            end
            bead_contact(i)=bead_contact_count;
        end
        contact_frames = full_frames(logical(wall_contact+bead_contact));
        if (nargout>=2)
            bead_contact_frames = full_frames(logical(bead_contact));
            if (nargout==3)
                wall_contact_frames = full_frames(logical(wall_contact));
            end
        end
    end
    function make_movie(obj, AYfig_in, watch_)
        frames = obj.Frames;
        if nargin==3
            AYfig_in.init_movie(frames, watch_);
        else
            AYfig_in.init_movie(frames);
        end

        walld = 5.72;
        wallL = (2/sqrt(3))*walld;
        wallv = [-wallL/2 -walld; -wallL 0; -wallL/2 walld ; wallL/2 walld; wallL 0; wallL/2 -walld; -wallL/2 -walld];
        wsca = max(abs([min(obj.dish(:, 2))-wallL, max(obj.dish(:, 2))+walld, min(obj.dish(:, 3))-wallL, max(obj.dish(:, 3))+walld]));
        lims = [-wsca, wsca, -wsca, wsca];

        rads = 0.5*ones(obj.beads,1);
        figdims = AYfig_in.get_dims();
        MS = figdims(4)/(4*wsca); %% radius in pixels
        SS = pi*MS*MS; %% area in pixels

        movie_fig = AYfig_in.fig;
        movie_ax = AYfig_in.ax;
        movie_gen = AYfig_in.movie_gen;
        dish = obj.dish;
        pos = obj.pos(:, 1:2, :);

        for i=1:frames
            plot(movie_ax, wallv(:, 1)+dish(i, 2), wallv(:, 2)+dish(i, 3), 'k -')
            hold(AYfig_in.ax, 'on');
            objdots = scatter(movie_ax, pos(:, 1, i), pos(:, 2, i), 'o', 'LineWidth', 1, 'SizeData', SS, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.1000 0.4440 0.2440], 'MarkerFaceAlpha', 0.5);
            txtbx = annotation(movie_fig, 'textbox', [0.8 0.9 0.1 0.1], 'String', num2str(i-1), 'LineStyle', 'none', 'FontSize', 16);
            hold(movie_ax, 'off');
            axis(movie_ax, lims);
            drawnow
            movie_gen(i) = getframe(movie_ax);
            delete(txtbx);
        end
        AYfig_in.movie_gen = movie_gen;
    end
  end
  methods (Static)
    function mat_pos = tens2mat(tens_)
        [beads, len_dat, Frames] = size(tens);
        mat_pos = reshape(tens_, [beads*len_dat, Frames]);
    end
    function tens_pos = mat2tens(mat_,beads,len_dat)
        tens_pos = reshape(mat_,[beads,len_dat,size(mat_,2)]);
    end
    function frscore = compute_frscore(fcap, tol, diff_full)
        f = 0;
        not_lost = true;
        while ((f<fcap)&&(not_lost))
            f=f+1;
            diff = diff_full(:, :, f);
            not_lost = max(sqrt(sum((diff.*diff)')))<tol;
        end
        frscore = f;
    end
  end
end
