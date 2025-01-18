%% sup fig 2AB
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue1late","cue1early"]},mouse_names,'sup_fig2B',tabular_sheet_name="sup_fig2A");


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup fig 2CD
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue2late","cue2early"]},mouse_names,'sup_fig2D',tabular_sheet_name="sup_fig2C");


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup fig 2E
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue1late","cue1early"]},mouse_names,'sup_fig2E');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup fig 2F
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue2late","cue2early"]},mouse_names,'sup_fig2F');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup fig 2GH
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.modality_between_context(cf,"cue_vs_cue_learning",mouse_names,'fig2GH')


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup fig 2I
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
CT_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
plot_colors = lines(2);
main_sr = common_functions.get_main_samplerate();

plot_info = {["G12",39],["G21",27]};
fig = figure(Position=[100,100,1200,400]);
tiled = tiledlayout(fig,1,3,TileSpacing="tight");
axes=gobjects(1,3);
for i=1:3
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:length(plot_info)
    this_ax = axes(pi);
    this_info = plot_info{pi};
    sr = main_sr(this_info(1));
    plot_x = (1:sr*4)/sr-1;

    tmp = cwa.cue2early.(this_info(1)).activity(:,str2num(this_info(2))+3,:);
    common_functions.plot_data_single(this_ax,plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
    tmp = cwa.cue2late.(this_info(1)).activity(:,str2num(this_info(2))+3,:);
    common_functions.plot_data_single(this_ax,plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));
    xline(this_ax,0,'--');yline(this_ax,0);
    xlim(this_ax,[-0.5,2]);
end
hold(axes,"off");
saveas(fig,[cf,'sup_fig2I.png'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup fig 2J
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G25","G26","G27"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue2late","cue2early"]},mouse_names,'sup_fig2J_smoothed',tabular_sheet_name="sup_fig2J",excel_only=true);
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue1late","cue1early"]},mouse_names,'sup_fig2J',cmap_limit={[-1,1]*0.018,[-1,1]*0.017,[-1,1]*0.014});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


