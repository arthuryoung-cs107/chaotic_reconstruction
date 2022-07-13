classdef debug_MH_plots
    properties

    end
    methods
        function obj = debug_MH_plots()
        end
    end
    methods (Static)
        function plot_rec_uerr(AYfig_,rec_list_,utrue_,base_color_)
            x=1:7;
            names = {'$$K$$';'$$\gamma_b$$';'$$\gamma_f$$';'$$\gamma_w$$';'$$\mu_b$$';'$$\mu_f$$';'$$\mu_w$$'};

            for i = 1:length(rec_list_)
                udiff_full=rec_list_(i).comp_udiff(utrue_);
                Y=udiff_full(x,:)./abs(utrue_(x));

                ax_=AYfig_.ax_tile(i);
                box(ax_,'on');
                hold(ax_,'on');
                plot(ax_,x,Y,'- o','Color',[base_color_, 0.1],'LineWidth',2);
                ylabel(ax_, '$$\frac{\hat{x}_i - x_i}{|\hat{x}_i|}$$', 'Interpreter', 'Latex', 'Fontsize', 14)
                set(ax_,'xtick',x,'xticklabel',names, 'TickLabelInterpreter', 'Latex', 'Fontsize', 14)
                ylim(ax_,[-4 1]);
            end
        end
    end
end
