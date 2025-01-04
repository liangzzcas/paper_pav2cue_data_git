%% figure 8C
close all;clear;clc;
cf = [pwd,'\'];
raw_path = [cf,'raw_data\'];
late_d3 = raw_path+["glu924\230219\data22_405_crop_MC_ROIs.mat","glu924\230219\data22_470_crop_MC_ROIs.mat",...
    "glu924\230219\glu924_pav2cue_d8_405_470_2023.02.19_09.32.05_ttlIn1_movie1_405.mat",...
    "glu924\230219\glu924_pav2cue_d8_405_470_2023.02.19_09.32.05_ttlIn1_movie1_470.mat"];
late_d5 = raw_path+["glu924\230221\data28_405_crop_MC_ROIs.mat","glu924\230221\data28_470_crop_MC_ROIs.mat",...
    "glu924\230221\glu_924_pav2cue_d10_405_470_2023.02.21_10.35.51_ttlIn1_movie1_405.mat",...
    "glu924\230221\glu_924_pav2cue_d10_405_470_2023.02.21_10.35.51_ttlIn1_movie1_470.mat"];

late_trial_info = {["d3","ch7",5431-17],["d5","ch7",9072-17]};

late_omi = "late";
this_info = late_trial_info;

plot_x = (1:540)/18;
trials = this_info;
n_trials = length(trials);
ax_height = 0.06; % change this to change plot density
unit_height = ax_height/(n_trials);

fig = figure('Position',[100,100,800,450]);
ax = axes(fig,"Position",[0.05,0.05,0.9,0.9]);
hold(ax,"on")
for r = 1:n_trials
    this_trial = trials{r};
    this_info_name = late_omi+"_"+this_trial(1);
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
    eval("this_paths = "+this_info_name+";");
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
    roi5 = load(this_paths(1));
    roi7 = load(this_paths(2));
    behav = load(this_paths(4));
    starting_frame = str2double(this_trial(3));
    plot(ax,plot_x,r*unit_height+ roi7.Fc(starting_frame:starting_frame+539,1) ,'-','Color',[0.13,0.65,0.47]);
    plot(ax,plot_x,r*unit_height+ roi5.Fc(starting_frame:starting_frame+539,1) ,'-','Color',[0.3,0.3,0.3]);
end
yline(ax,[1:n_trials]*unit_height)
xline(ax,0)
hold(ax,"off")

set(gca,'YColor',"none") % comment this out to show Y axis

set(gca,'YLim',[0.2*unit_height,ax_height+0.8*unit_height])
set(gca,'XLim',plot_x([1,end])-[1/18,0])
saveas(fig,[pwd,'\fig8C.png'],'png')
delete(fig);


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 8D
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\Glu_task_ta.mat']);

p_names_to_plot = ["post","LED_omi"]; % phase used to get trace plot
p_names_to_plot_1 = ["post","LED_omi_complete"];
p_text = ["Post-learning","Extinction 20%"];

main_sr = 18;
plot_x = (1:4*main_sr)/main_sr-1;
plot_colors = lines(2);

lgs=gobjects(1,2);
fig = figure(Position=[100,100,1600,500]);
tiled=tiledlayout(1,3,TileSpacing="tight");
axs = gobjects(1,3);
for i=1:3
    axs(i)=nexttile(tiled,i);
end
hold(axs,"on")

for i=1:length(p_names_to_plot)
    this_act = permute(cwa.glu924.(p_names_to_plot(i)).cueOn.cue1.activity(:,4,:),[1,3,2]);
    tmp = this_act-mean(this_act(1:main_sr,:),"all","omitmissing");
    lgs(i)=common_functions.plot_data_single(axs(1),plot_x,tmp,plot_color=plot_colors(i,:),sem_color=plot_colors(i,:)*0.8);

    tmp_clim = max(abs(prctile(tmp,[2,98],"all")));
    imagesc(axs(i+1),plot_x,1:size(tmp,2),tmp',[-1,1]*tmp_clim);
    colormap(axs(i+1),"redblue")
end
xline(axs(1),0,'k--');yline(axs(1),0);
legend(axs(1),lgs,p_text)

for i=1:length(p_names_to_plot_1)
    this_act = permute(cwa.glu924.(p_names_to_plot_1(i)).cueOn.cue1.activity(:,4,:),[1,3,2]);
    tmp = this_act-mean(this_act(1:main_sr,:),"all","omitmissing");
    tmp_clim = max(abs(prctile(tmp,[2,98],"all")));
    imagesc(axs(i+1),plot_x,1:size(tmp,2),tmp',[-1,1]*0.0045);
    colormap(axs(i+1),"redblue")
    axs(i+1).YDir='reverse';
    ylim(axs(i+1),[1,size(tmp,2)]); xlim(axs(i+1),[-1,3]);
    xline(axs(i+1),0);
    title(axs(i+1),"light "+p_text(i));
end
colorbar(axs(3));

hold(axs,"off")
saveas(fig,[pwd,'\fig8D.png'],'png')
delete(fig);

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 8D and 8E
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\Glu_task_ta.mat']);
mouse_names = string(fields(cwa))';

p_names_to_plot = ["post","LED_omi_complete"]; % phase used to get trace plot
p_text = ["Post-learning","Extinction 20%"];
across_roi_struct = struct; % pre-allocate across ROIs data (personally don't favor combining ROIs)
p_names = p_names_to_plot;
for pi = 1:length(p_names)
    p_name = p_names(pi);
    across_roi_struct.(p_name).cue1 = [];
    across_roi_struct.(p_name).cue2 = [];
    across_roi_struct.(p_name).cue1_session_level = [];
    across_roi_struct.(p_name).cue2_session_level = [];
end

vel_smooth_info = {"lowess",1/2};
main_sr = 18;
plot_x = (1:4*main_sr)/main_sr-1;
plot_color = lines(3);
fiber_text = ["lick","velocity","Fc"];
for mouse_name = mouse_names
    for pi = 1:length(p_names)
        p_name = p_names(pi);
        if ~isfield(cwa.(mouse_name).(p_name),"cueOn") % doesn't have this phase, skip
            continue
        end
        for ci = 1:2
            c = "cue"+ci;
            this_act = cwa.(mouse_name).(p_name).cueOn.(c).activity;
            n_trial = size(this_act,3);

            % replace angular_v by linear_decel
            lin_v = permute(this_act(:,2,:),[1,3,2]);
            smoothed_linv = smoothdata(lin_v,1,vel_smooth_info{1},round(vel_smooth_info{2}*main_sr));
            decel = diff(smoothed_linv,1,1)*main_sr; decel = cat(1,decel,decel(end,:));
            smoothed_decel = smoothdata(decel,1,vel_smooth_info{1},round(vel_smooth_info{2}*main_sr));
            this_act(:,2,:) = smoothed_linv;
            this_act(:,3,:) = smoothed_decel;
            % rebaseline Fc
            tmp = mean(permute(this_act(1:main_sr,4,:),[1,3,2]),"all","omitmissing");
            this_act(:,4,:) = this_act(:,4,:)-tmp;
            % save activity to across_ROI stuct
            across_roi_struct.(p_name).(c) = cat(3,across_roi_struct.(p_name).(c),this_act);
            %   - get session level
            session_n_trial = cwa.(mouse_name).(p_name).cueOn.(c).session_num_trials;
            this_act_session = mat2cell(this_act,size(this_act,1),size(this_act,2),session_n_trial);
            this_act_session = cell2mat(cellfun(@(x) mean(x,3,"omitmissing"),this_act_session,UniformOutput=false));
            across_roi_struct.(p_name).(c+"_session_level") = cat(3,across_roi_struct.(p_name).(c+"_session_level"),this_act_session);
        end
    end
end

level_text = ["session"];
for ave_i=1 % loop thru trial or session level
    fig3 = figure(position=[100,100,1600,900]);
    tiled3 = tiledlayout(fig3,2,3);
    sgtitle(tiled3,"Glu across two fibers "+level_text(ave_i)+" level")
    axs3 = gobjects(1,6);
    for i=1:3
        axs3(i) = nexttile(tiled3,i);
        axs3(i+3) = nexttile(tiled3,i+3);
        title(axs3(i),"cue1 "+fiber_text(i));
        title(axs3(i+3),"cue2 "+fiber_text(i));
    end
    hold(axs3,"on")
    plot_r = [1,2,4];
    for ri=1:3
        r=plot_r(ri);
        ax1 = axs3(ri); ax2 =  axs3(ri+3);
        ps = gobjects(1,2);
        for pi = 1:length(p_names)
            p_name = p_names(pi);
            ps(pi)=common_functions.plot_data_single(ax1,plot_x,across_roi_struct.(p_name).cue1_session_level(:,r,:),plot_color=plot_color(pi,:)*0.8,sem_color=plot_color(pi,:));
            common_functions.plot_data_single(ax2,plot_x,across_roi_struct.(p_name).cue2_session_level(:,r,:),plot_color=plot_color(pi,:)*0.8,sem_color=plot_color(pi,:));
        end
        if ri==1
            legend(ax1,ps,p_text,AutoUpdate="off")
        end
    end
    hold(axs3,"off")
    for ax = axs3
        xline(ax,0,'k--'); yline(ax,0);
    end
    for i=1:3
        linkaxes([axs3(i),axs3(i+3)],"xy");
    end
    saveas(fig3,[pwd,'\fig8E.png'],'png')
    delete(fig3);
end

