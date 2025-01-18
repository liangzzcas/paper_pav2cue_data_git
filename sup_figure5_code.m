%% sup figure 5A
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue1LEDomi","cue1Toneomi"]},mouse_names,'sup_fig5A_smoothed',tabular_sheet_name="sup_fig5A",excel_only=true);
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue1LEDomi","cue1Toneomi"]},mouse_names,'sup_fig5A',cmap_limit={[-1,1]*0.011,[-1,1]*0.012,[-1,1]*0.011});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 5B
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue2Toneomi","cue2LEDomi"]},mouse_names,'sup_fig5B_smoothed',tabular_sheet_name="sup_fig5B",excel_only=true);
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue2Toneomi","cue2LEDomi"]},mouse_names,'sup_fig5B',cmap_limit={[-1,1]*0.011,[-1,1]*0.012,[-1,1]*0.011});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

