classdef debug_MH_plots
    properties (Constant)
        posdim_full = [1 1 1728 1000];
        posdim_toprow = [1 551 1728 460];
        posdim_bottomrow = [0 1 1728 460];

        uerr_range=[-4,1];
    end
    methods
        function obj = debug_MH_plots()
        end
    end
    methods (Static)
        function AYfig_out = make_fig(name_,posdim_)
            AYfig_out = AYfig(AYfig.specs_gen(name_,posdim_),false);
        end
        function AYfig_out = make_fullfig(name_)
            AYfig_out = debug_MH_plots.make_fig(name_,debug_MH_plots.posdim_full);
        end
        function plot_recblock_uerr(ax_,recblock_,utrue_,base_color_)
            x=1:7;
            names = {'$$K$$';'$$\gamma_b$$';'$$\gamma_f$$';'$$\gamma_w$$';'$$\mu_b$$';'$$\mu_f$$';'$$\mu_w$$'};
            udiff_full=recblock_.comp_udiff(utrue_);
            Y=udiff_full(x,:)./abs(utrue_(x));
            box(ax_,'on');
            hold(ax_,'on');
            plot(ax_,x,Y,'- o','Color',[base_color_, 0.1],'LineWidth',2);
            ylabel(ax_, '$$\frac{\hat{x}_i - x_i}{|\hat{x}_i|}$$', 'Interpreter', 'Latex', 'Fontsize', 14)
            set(ax_,'xtick',x,'xticklabel',names, 'TickLabelInterpreter', 'Latex', 'Fontsize', 14)
            ylim(ax_,debug_MH_plots.uerr_range);
        end
        function plot_Frame_bead_data(axs_,Y_,base_color_)
            Y=permute(Y_,[2,1,3]); %% pool,Frame,bead->Frame,pool,bead

            Frames=size(Y,1);
            for i = 1:length(axs_)
                plot(axs_(i),1:Frames,Y(:,:,i),'-','Color',[base_color_,0.1],'LineWidth',2);
            end
            set(axs_,'YScale','log')
        end

        function plot_reclist_uerr(AYfig_,rec_list_,utrue_,base_color_)
            for i = 1:length(rec_list_)
                debug_MH_plots.plot_recblock_uerr(AYfig_.ax_tile(i),rec_list_(i),utrue_,base_color_)
            end
        end
        function plot_Class_diagnostics(AYfig_, mhgen_, Classi_, utrue_)
            run AYfigprops.m
            AYfig_.init_tiles([1,3]);
            debug_MH_plots.plot_recblock_uerr(AYfig_.ax_tile(1),Classi_.leaders,utrue_,green4);
        end
        function plot_gen_diagnostics(AYfig_, mhgen_, geni_, utrue_)
            run AYfigprops.m
            AYfig_.init_tiles([1,3]);
            debug_MH_plots.plot_recblock_uerr(AYfig_.ax_tile(1),geni_.pool,utrue_,green4);
        end
        function plot_event_diagnostics(AYfig_, mhgen_, eventi_, utrue_)
            run AYfigprops.m
            AYfig_.init_tiles([2,3]);
            debug_MH_plots.plot_recblock_uerr(AYfig_.ax_tile(1),eventi_.pool,utrue_,green4);

            [r2_p_F_b, alpha_p_F_b]=eventi_.get_pool_Frame_bead();

            % debug_MH_plots.plot_Frame_bead_data(AYfig_.ax_tile(4:6),r2_p_F_b,red5);
            debug_MH_plots.plot_Frame_bead_data(AYfig_.ax_tile(4:6),alpha_p_F_b,purple1);

            title(AYfig_.ax_tile(1), '$$\vec{u}$$ error', 'Interpreter', 'Latex', 'Fontsize', 14);
        end
    end
end
