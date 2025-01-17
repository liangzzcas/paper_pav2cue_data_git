%% figure 3A
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_merged.mat']);
% cwa = load("E:\Safa_Processed\#paper\#paper_figures\_data\highpass03_wGLM_v2\delivery\components_window_activity_filtered_rebase_rew_merged.mat");
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

plot_info_1 = {"G12",39,30;"G21",27,18;};
plot_trace_clim = {[-0.01,0.02];[-0.025,0.025];};
plot_imagesc_clim = {0.03;0.04;};
plot_info_2 = ["rew1early","rew1late"];
plot_colors = lines(length(plot_info_2));

fig = figure(Position=[40,40,1800,1310]);
tiled = tiledlayout(fig,length(plot_info_2)+1,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,(length(plot_info_2)+1)*size(plot_info_1,1));
for i=1:(length(plot_info_2)+1)*size(plot_info_1,1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axes([pi+size(plot_info_1,1).*[0:length(plot_info_2)]]);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    
    plot_x = (1:sr*4)/sr-1;
    nts = nan(1,length(plot_info_2));
    lgs = gobjects(1,length(plot_info_2));
    for ni=1:length(plot_info_2)
        tmp = cwa.(plot_info_2(ni)).(this_info{1}).activity(:,this_info{2}+3,:); nts(ni) = size(tmp,3);
        lgs(ni)=common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(ni,:)*0.8,sem_color=plot_colors(ni,:));
        imagesc(this_axes(ni+1),plot_x,1:nts(ni),permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
        this_axes(ni+1).YDir = "reverse";
        colormap(this_axes(ni+1),common_functions.redblue());
        colorbar(this_axes(ni+1));
        coord={ct_table{ct_table.mouse_name==this_info{1} & ct_table.ROI_original==this_info{2},["fiber_bottom_AP","fiber_bottom_ML","fiber_bottom_DV"]}};
        title(this_axes(1),this_info{1}+" ROI "+this_info{2}+sprintf(" AP/ML/DV=%0.2f/%0.2f/%0.2f",coord{:}))
    end

    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-0.5,2]);
    ylim(this_axes(1),plot_trace_clim{pi});
    for axi=1:length(plot_info_2)
        ax=this_axes(axi+1);
        xline(ax,0,'--');
        xlim(ax,plot_x([1,end]));
        ylim(this_axes(axi+1),[1,nts(axi)]);
        ylabel(this_axes(axi+1),plot_info_2(axi))
    end
end
hold(axes,"off");
legend(axes(1),lgs,plot_info_2);
saveas(fig,[cf,'fig3A.fig']);
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 3BC
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_rewmerged",{["rew1late","rew1early"]},mouse_names,'fig3C',tabular_sheet_name="fig3B");

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 3DE
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.modality_between_context(cf,"cue_vs_rew",mouse_names,'fig3DE')


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

