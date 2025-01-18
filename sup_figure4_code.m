%% sup figure 4A
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue1LEDomi","cue1late"]},mouse_names,'sup_fig4A',cmap_limit={[-1,1]*0.011,[-1,1]*0.012,[-1,1]*0.011});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 4B
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue2Toneomi","cue2late"]},mouse_names,'sup_fig4B',cmap_limit={[-1,1]*0.011,[-1,1]*0.012,[-1,1]*0.011});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 4C
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue2Toneomi","cue2late"]},mouse_names,'sup_fig4C',tabular_sheet_name="sup_fig4C");


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 4D
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G25","G26","G27"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue1LEDomi","cue1late"]},mouse_names,'sup_fig4D_smoothed',tabular_sheet_name="sup_fig4D",excel_only=true);
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue1LEDomi","cue1late"]},mouse_names,'sup_fig4D',cmap_limit={[-1,1]*0.014,[-1,1]*0.013,[-1,1]*0.013});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 4EF
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.modality_between_context(cf,"cue_vs_cue_extinction",mouse_names,'fig4EF')


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4G
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
% plot_info_1 = {"G23",36,18;"G23",48,18;"G21",41,18;"G21",9,18;"G12",5,30;};
plot_info_1 = {"G23",36,18;"G23",48,18;"G12",5,30;};
plot_info_2 = ["cue2late","cue2LEDomi","cue2Toneomi"];
plot_colors = lines(3);

fig = figure(Position=[100,100,1600,300]);
nr=1;
tiled = tiledlayout(fig,nr,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,nr*size(plot_info_1,1));
for i=1:nr*size(plot_info_1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axes(pi);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    
    plot_x = (1:sr*4)/sr-1;
    tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
    
    tmp = cwa.(plot_info_2(2)).(this_info{1}).activity(:,this_info{2}+3,:);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));

    tmp = cwa.(plot_info_2(3)).(this_info{1}).activity(:,this_info{2}+3,:);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(3,:)*0.8,sem_color=plot_colors(3,:));
   
    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-0.5,2]);
end
hold(axes,"off");
saveas(fig,'sup_fig4G.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

