%% figure 1C: GRAB-ACh raw traces
% loading % setting
close all;clear;clc;
cf0 = [pwd,'\'];
cf = [cf0,'raw_data\'];
settings = struct;
settings.paper_v1_data.mouse_name = "G19";
settings.paper_v1_data.behav = load([cf,'\G19\220504\patched_behav\behav_470.mat']);
settings.paper_v1_data.roi_405 = load([cf,'\G19\220504\Data163_405_crop_MC_ROIs.mat']).Fc;
settings.paper_v1_data.roi_470 = load([cf,'\G19\220504\Data163_470_crop_MC_ROIs.mat']).Fc;
settings.paper_v1_data.sr = 18;
settings.paper_v1_data.roi_nums = [49,46,22,43,50,25,26,51,40,38,56,20,14,17];
settings.paper_v1_data.frame_window = 1100:2180;

settings.paper_v1_mute.mouse_name = "G25";
settings.paper_v1_mute.behav = load ([cf,'G25\220926\patched_behav\behav470.mat']);
settings.paper_v1_mute.roi_405 = load([cf,'G25\220926\Data172_405_crop_MC_ROIs.mat']).Fc;
settings.paper_v1_mute.roi_470 = load([cf,'G25\220926\Data172_470_crop_MC_ROIs.mat']).Fc;
settings.paper_v1_mute.sr = 18;
settings.paper_v1_mute.roi_nums = [10,18,6,12,3,17,19,33,21,8,61,42,25,60,30];
settings.paper_v1_mute.frame_window = 3440:4520;

this_setting = settings.paper_v1_data;
% this_setting = settings.paper_v1_mute;

mouse_name = this_setting.mouse_name;
behav = this_setting.behav;
roi_405 = this_setting.roi_405;
roi_470 = this_setting.roi_470;
sr = this_setting.sr;
roi_nums = this_setting.roi_nums;
frame_window = this_setting.frame_window;
vel = common_functions.ball2xy(behav);

seg = vel.angular_velocity_ms(frame_window,1);
plot_x = (frame_window-frame_window(1))/sr;
n_rois = numel(roi_nums);

% actual plot
ax_height = 1.1; % change this to change plot density

unit_height = ax_height/(n_rois);
seg_scaled = seg/max(seg,[],"omitnan")*unit_height - 0.5*unit_height;
min_seg_scaled = min(seg_scaled,[],"omitnan");

reward = behav.reward(frame_window)*(ax_height+0.5*unit_height); reward(reward==0) = min_seg_scaled;
stim_LED = behav.stimulus_led(frame_window)*(ax_height+0.5*unit_height); stim_LED(stim_LED==0) = min_seg_scaled;
stim_sound = behav.stimulus_sound2(frame_window)*(ax_height+0.5*unit_height); stim_sound(stim_sound==0) = min_seg_scaled;

fig = figure('Position',[50 50 560 1032]);
ax1 = axes(fig,"Position",[0.05,0.05,0.9,0.9]);
hold(ax1,"on")
for r = 1:n_rois
    roi_num = roi_nums(r);
    plot(ax1,plot_x,r*unit_height+roi_470(frame_window,roi_num),'-','Color',[0.13,0.65,0.47])
    % plot(ax1,plot_x,r*unit_height+roi_405(frame_window,roi_num),'-','Color',[0.0 0.0 .0])
end
% plot(ax1,plot_x,seg_scaled,'-k','LineWidth',0.1)

% plot(ax1,plot_x,reward,'Color','k')
% plot(ax1,plot_x,stim_LED,'Color','b')
% plot(ax1,plot_x,stim_sound,'Color','r')
% yline(ax1,[-0.5,1:n_rois]*unit_height)
hold(ax1,"off")
% title(ax1,mouse_name)

set(gca,'YColor',"none") % comment this out to show Y axis
set(gca,'YLim',[min_seg_scaled,ax_height+0.5*unit_height])
set(gca,'XLim',plot_x([1,end]))
saveas(fig,[cf0,'fig1C.png'],'png')
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 1E: fiber locations of all mice
% mouse implant coverage map, circle plot containing all mice, each mice with independent color
close all;clear;clc;
cf0 = [pwd,'\'];
cf = [cf0,'raw_data\'];
CT_table = readtable([cf,'CT_across_GXX_mice.xlsx']);

mouse_names = ["G12","G15","G17","G19","G22","G21","G23","G24"];
this_CT_table = CT_table(contains(CT_table{:,"mouse_name"},mouse_names),:);
mouse_data = nan(size(this_CT_table,1),1);
for i = 1:length(mouse_names)
    mouse_data(string(this_CT_table{:,"mouse_name"})==mouse_names(i)) = i;
end

fig = figure(Position=[100,100,1600,900]);
tiled = tiledlayout(fig,1,2);
axs = gobjects(1,2);
for i=1:2
    axs(i) = nexttile(tiled,i);
end
title(axs(1),"implant covearge map of "+strjoin(mouse_names," "))
hold(axs,"on")
ps = common_functions.scatter_3d(axs(1),mouse_data,this_CT_table,skipBubble=~this_CT_table{:,"significance"},...
    colormapOption='jet',colormapBins=length(mouse_names),setuserdata=0,viewAngle=[0,90]);

cb = fig.Children.Children(1);
cb.Ticks = 1/length(mouse_names)/2 + (0:(length(mouse_names)-1))/length(mouse_names);
cb.TickLabels = mouse_names;
axis(axs(1),"vis3d")
ps = circle_UI_function.scatter_3d(axs(2),mouse_data,this_CT_table,skipBubble=~this_CT_table{:,"significance"},...
    colormapOption='jet',colormapBins=length(mouse_names),setuserdata=0,viewAngle=[-90,0]);
axis(axs(2),"vis3d")
hold(axs,"off")
saveas(fig,[cf0,'\fig1E.png'],'png')
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 1F: unpred across weeks
close all;clear;clc;
cf0 = [pwd,'\'];
cf = [pwd,'\raw_data\'];
this_ta = load([cf,'event_aligned_highpass03.mat']).cwa_raw.G19;
plot_roi = 52;
ta_consump_unpred_weekly = struct;

out_ta = structfun(@(x) x.rewOn.unpred.rew,this_ta,UniformOutput=false);
ta_in = out_ta; main_sr = 18;
f_names = string(fields(ta_in)');
n_days = length(f_names);

batches = {1:6,7:12,13:18};

output = struct;
for i = 1:length(batches)
    output.("week"+i) = struct;
    this_activity = [];
    batch_f_names = batches{i};
    for batch_f_name = "file"+batch_f_names
        this_act = ta_in.(batch_f_name).activity;
        this_sr = size(this_act,1)/4;
        if this_sr~= main_sr
            this_act = common_functions.interp_ta(this_act,this_sr,main_sr);
        end
        this_activity = cat(3, this_activity, this_act);
    end
    output.("week"+i).activity = this_activity;
    output.("week"+i).mu = mean(output.("week"+i).activity,3,"omitmissing");
    output.("week"+i).sem = std(output.("week"+i).activity,[],3,"omitmissing")/sqrt(size(output.("week"+i).activity,3));
end

f_names_week = string(fields(output)');
plot_x = 1:main_sr*4;
plot_x = plot_x/main_sr-1;

weekly_color = parula(length(batches));
fig = figure(Position=[100,100,1600,900]); % by session
ax1 = axes(fig);
hold(ax1,"on")
for s_name_i = 1:length(f_names_week)
    s_name = f_names_week(s_name_i);
    common_functions.plot_ta_single_just_mu_sem(ax1,plot_x,output.(s_name),plot_roi+3,...
        plot_color=weekly_color(s_name_i,:)*0.8,sem_color=weekly_color(s_name_i,:));
end
hold(ax1,"off")
xline(ax1,0);yline(ax1,0);
saveas(fig,[cf0,'fig1F.png'],'png')
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 1G
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_time_histogram(cf,"unpred",'fig1G')


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 1H
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_presence_circle(cf,"pav2cue task","unpred",'fig1H');

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 1I
close all;clear;clc;
cf = [pwd,'\'];
event_ta = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']).cwa_raw;
CT_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
consump_window = common_functions.get_comp_timewindow().unpred;

mouse_names = ["G19"];
plot_roi = 52;
output = struct;
for mouse_name = mouse_names
    rois_to_plot = find(CT_table{CT_table{:,"mouse_name"}==mouse_name,"significance"}');
    mouse_data = event_ta.(mouse_name);
    n_days = length(fields(mouse_data));
    n_rois = length(rois_to_plot);
    this_comps = struct;
    this_comps.rois = rois_to_plot;
    this_comps.mus = nan(3,n_rois,n_days);
    this_comps.sems = nan(3,n_rois,n_days);
    this_comps.locs = nan(3,n_rois,n_days);
    this_comps.activity = cell(3,n_rois,n_days);
    for d_i = 1:n_days
        this_data = mouse_data.("file"+d_i).rewOn.unpred.rew_consump;
        n_trials = size(this_data.activity,3);

        this_sr = size(this_data.activity,1)/4; % assume 4 seconds
        this_window = round(consump_window*this_sr);

        [M_pk,I_pk] = max(this_data.mu(this_window(1,1):this_window(1,2),rois_to_plot+3),[],1,"omitmissing"); I_pk = I_pk+this_window(1,1)-1;
        [M_dp,I_dp] = min(this_data.mu(this_window(2,1):this_window(2,2),rois_to_plot+3),[],1,"omitmissing"); I_dp = I_dp+this_window(2,1)-1;
        [M_re,I_re] = max(this_data.mu(this_window(3,1):this_window(3,2),rois_to_plot+3),[],1,"omitmissing"); I_re = I_re+this_window(3,1)-1;
        
        this_act_pk = nan(n_rois,n_trials);
        this_act_dp = nan(n_rois,n_trials);
        this_act_re = nan(n_rois,n_trials);

        for r_i = 1:length(rois_to_plot)
            this_act_pk(r_i,:) = this_data.activity(I_pk(r_i),rois_to_plot(r_i)+3,:);
            this_act_dp(r_i,:) = this_data.activity(I_dp(r_i),rois_to_plot(r_i)+3,:);
            this_act_re(r_i,:) = this_data.activity(I_re(r_i),rois_to_plot(r_i)+3,:);
        end
        
        this_comps.mus(1,:,d_i) = M_pk;
        this_comps.locs(1,:,d_i) = I_pk;
        this_comps.sems(1,:,d_i) = std(this_act_pk,[],2,"omitmissing")./sqrt(n_trials);
        this_comps.activity(1,:,d_i) = num2cell(this_act_pk,2);

        this_comps.mus(2,:,d_i) = M_dp;
        this_comps.locs(2,:,d_i) = I_dp;
        this_comps.sems(2,:,d_i) = std(this_act_dp,[],2,"omitmissing")./sqrt(n_trials);
        this_comps.activity(2,:,d_i) = num2cell(this_act_dp,2);

        this_comps.mus(3,:,d_i) = M_re;
        this_comps.locs(3,:,d_i) = I_re;
        this_comps.sems(3,:,d_i) = std(this_act_re,[],2,"omitmissing")./sqrt(n_trials);
        this_comps.activity(3,:,d_i) = num2cell(this_act_re,2);
    end
    output.(mouse_name) = this_comps;

    % build anova group tag
    anova_mouse_data = output.(mouse_name);
    tmp_act_r = squeeze(anova_mouse_data.activity(:,1,:));
    group_tags = [];
    for d_i = 1:n_days
        group_tags = cat(2,group_tags, repmat("day"+d_i,[1,length(tmp_act_r{1,d_i})]) );
    end
    output.(mouse_name).anova1_group_tags = group_tags;
    output.(mouse_name).anova1_tbl = cell(3,n_rois);
    % build anova data
    for r_i = 1:n_rois
        act_r = squeeze(anova_mouse_data.activity(:,r_i,:));
        for comp_i = 1:3
            act_this_comp = cell2mat(act_r(comp_i,:));
            [p,tbl,stats] = anova1(act_this_comp,group_tags,"off");
            output.(mouse_name).anova1_tbl(comp_i,r_i) = cell({tbl});
        end
    end
end

mouse_names = "G19";
line_colors = lines(3);
comp_text = ["pk","dp","re"];
for mouse_name = mouse_names
    mouse_data = output.(mouse_name);
    n_days = size(mouse_data.mus,3);
    n_rois = length(mouse_data.rois);

    plot_x = 1:n_days;
    title_text = string(mouse_data.rois);

    group_tags = mouse_data.anova1_group_tags;
    % for r_i = 1:n_rois
    for r_i = find(mouse_data.rois==plot_roi)
        mu_r = squeeze(mouse_data.mus(:,r_i,:));
        sem_r = squeeze(mouse_data.sems(:,r_i,:));
        act_r = squeeze(mouse_data.activity(:,r_i,:));

        ps = gobjects(1,3);
        fig = figure(Position=[100,100,600,500]);
        tiled = tiledlayout(fig,1,1);
        sgtitle(tiled,mouse_name+" ROI "+mouse_data.rois(r_i)+" unpred rew components across days")
        axs = gobjects(1,1);
        for i = 1:1
            axs(i) = nexttile(tiled,i);
        end
        hold(axs,"on")
        for comp_i = 1:3
            ps(comp_i) = errorbar(axs(1),plot_x,mu_r(comp_i,:),sem_r(comp_i,:),Color=line_colors(comp_i,:),LineWidth=1.5);
        end
        hold(axs,"off")
        xlim(axs,[0,n_days])
        xlabel(axs,"day")
        ylabel(axs,"DF/F")
        legend(axs,ps,comp_text)
        linkaxes(axs,"xy")

        saveas(fig,[cf,'fig1I.fig'])
        delete(fig)
    end
end