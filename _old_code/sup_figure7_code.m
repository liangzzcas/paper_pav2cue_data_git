%% figure 7A-B
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
cwa_all_days = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']).cwa_raw;
cwa_all_days_fields = ["cueOn","cue1","activity"];
session_info = common_functions.get_training_info();
plot_info_1 = {"G23",1,18;"G23",2,18;"G21",1,18;"G21",2,18;};
plot_info_2 = ["cue1late","cue1LEDomi","cue1Toneomi"];
plot_info_3 = {11:15,16:20;11:15,16:20;16:20,11:15;16:20,11:15;};

vel_smooth_info = {"lowess",1/2};

fig = figure(Position=[100,100,1600,900]);
nr=3;
tiled = tiledlayout(fig,nr,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,nr*size(plot_info_1,1));
for i=1:nr*size(plot_info_1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    if pi == 2 || pi == 4
        acc_flag = true;
    else
        acc_flag = false;
    end

    this_axes = axes([pi,pi+size(plot_info_1,1),pi+2*size(plot_info_1,1)]);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    plot_x = (1:sr*4)/sr-1;
    tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2},:); nt1 = size(tmp,3);
    if acc_flag
        vel_smoothed = smoothdata(squeeze(tmp),1,vel_smooth_info{1},round(vel_smooth_info{2}*sr))';
        decel = diff(vel_smoothed,1,2)*sr; decel = cat(2,decel,decel(:,end));
        tmp(:,1,:) = smoothdata(decel,2,vel_smooth_info{1},round(vel_smooth_info{2}*sr))';
        this_clim = [-2.5,2.5];
        xlcolor = [0,0,0];
    else
        tmp = tmp/max(tmp,[],"all","omitmissing");
        this_clim = [0,1];
        xlcolor = [1,1,1];
    end
    imagesc(this_axes(1),plot_x,1:nt1,permute(tmp,[3,1,2]),this_clim);
    this_axes(1).YDir = "reverse";
    if acc_flag
        colormap(this_axes(1),common_functions.redblue());
    end
    colorbar(this_axes(1));

    % get all days for omission
    tmp = arrayfun(@(x) getfield(cwa_all_days.(this_info{1}).(x),cwa_all_days_fields{:}),"file"+plot_info_3{pi,1}, UniformOutput=false);
    tmp = cell2mat(reshape(tmp,1,1,[])); tmp = tmp(:,this_info{2},:);
    if acc_flag
        vel_smoothed = smoothdata(squeeze(tmp),1,vel_smooth_info{1},round(vel_smooth_info{2}*sr))';
        decel = diff(vel_smoothed,1,2)*sr; decel = cat(2,decel,decel(:,end));
        tmp(:,1,:)  = smoothdata(decel,2,vel_smooth_info{1},round(vel_smooth_info{2}*sr))';
        this_clim = [-2.5,2.5];
        xlcolor = [0,0,0];
    else
        tmp = tmp/max(tmp,[],"all","omitmissing");
        this_clim = [0,1];
        xlcolor = [1,1,1];
    end
    nt2 = size(tmp,3);
    imagesc(this_axes(2),plot_x,1:nt2,permute(tmp,[3,1,2]),this_clim);
    this_axes(2).YDir = "reverse";
    if acc_flag
        colormap(this_axes(2),common_functions.redblue());
    end
    colorbar(this_axes(2));
    



    % get all days for omission
    tmp = arrayfun(@(x) getfield(cwa_all_days.(this_info{1}).(x),cwa_all_days_fields{:}),"file"+plot_info_3{pi,2}, UniformOutput=false);
    tmp = cell2mat(reshape(tmp,1,1,[])); tmp = tmp(:,this_info{2},:);
    if acc_flag
        vel_smoothed = smoothdata(squeeze(tmp),1,vel_smooth_info{1},round(vel_smooth_info{2}*sr))';
        decel = diff(vel_smoothed,1,2)*sr; decel = cat(2,decel,decel(:,end));
        tmp(:,1,:) = smoothdata(decel,2,vel_smooth_info{1},round(vel_smooth_info{2}*sr))';
        this_clim = [-2.5,2.5];
        xlcolor = [0,0,0];
    else
        tmp = tmp/max(tmp,[],"all","omitmissing");
        this_clim = [0,1];
        xlcolor = [1,1,1];
    end
    nt3 = size(tmp,3);
    imagesc(this_axes(3),plot_x,1:nt3,permute(tmp,[3,1,2]),this_clim);
    this_axes(3).YDir = "reverse";
    if acc_flag
        colormap(this_axes(3),common_functions.redblue());
    end
    colorbar(this_axes(3));

    for ax=this_axes
        xline(ax,0,'-',Color=xlcolor);
        xlim(ax,plot_x([1,end]));
    end
    ylim(this_axes(1),[1,nt1]);
    ylim(this_axes(2),[1,nt2]);
    ylim(this_axes(3),[1,nt3]);
end
hold(axes,"off");
saveas(fig,'sup_fig7A_B.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7C
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue1LEDomi","cue1Toneomi"]},mouse_names,'sup_fig7C',cmap_limit={[-1,1]*0.011,[-1,1]*0.012,[-1,1]*0.011});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7D
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue2Toneomi","cue2LEDomi"]},mouse_names,'sup_fig7D',cmap_limit={[-1,1]*0.011,[-1,1]*0.012,[-1,1]*0.011});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%