%% sup figure 7A
close all;clear;clc;
cf = [pwd,'\'];
cwa_data = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
main_srs = common_functions.get_main_samplerate();
vel_window = [0,0.6]; % velocity window of getting acceleration etc
vel_smooth_info = {"sgolay",1/4};
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
cue_names = ["LED","Tone"];

% build acceleration data
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = cell(1,16);
output_tbl_varnames = reshape(["peak_decel_light_","peak_decel_tone_","vel_change_light_","vel_change_tone_",]+["pre";"post";"ex_10";"ex_20"],1,[]);
%--------------------------------------------------- output file ---------------------------------------------------%
vel_struct = struct;
test_data = struct;
for c_i = 1:2
    cue_text = "cue"+c_i;
    phase_names = [cue_text+"early",cue_text+"late",cue_text+cue_names(c_i)+"omi",cue_text+cue_names(3-c_i)+"omi"];
    test_data.(cue_text) = nan(length(phase_names),8,2);
    for m_i = 1:8
        mouse_name = mouse_names(m_i);
        this_sr = main_srs(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr;
        vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        for p_i = 1:length(phase_names)
            vel = cwa_data.(phase_names(p_i)).(mouse_name).mu(:,2);
            vel_smoothed = smoothdata(vel,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));

            vel_info = common_functions.find_acc_start_end_tbl(vel_smoothed,vel_window_frame,this_sr);
            vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_original = vel;
            vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_smoothed = vel_smoothed;
            vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_info = vel_info;

            test_data.(cue_text)(p_i,m_i,:) = [-vel_info{1,"acc_peak_value"},vel_info{1,"vel_start"}-vel_info{1,"vel_end"}];
        end
    end
    %--------------------------------------------------- output file ---------------------------------------------------%
    output_tbl(1,[1:4,9:12]+4*(c_i-1)) = [num2cell(test_data.(cue_text)(:,:,1)',1),num2cell(test_data.(cue_text)(:,:,2)',1)];
    %--------------------------------------------------- output file ---------------------------------------------------%
end
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = table(output_tbl{:},VariableNames=output_tbl_varnames);
%--------------------------------------------------- output file ---------------------------------------------------%


testf = @kruskalwallis;
test_name = "kruskalwallis";
t_tags = ["decel","diff_v"];
t_texts = ["deceleration","change in velocity"];
test_result = struct;
fig_test = figure(position=[100,100,1600,900]);
tiled_test = tiledlayout(fig_test,2,2);
axs = gobjects(1,4);
for ax_i = 1:4
    axs(ax_i) = nexttile(tiled_test,ax_i);
end
hold(axs,"on")
for c_i = 1:2
    cue_text = "cue"+c_i;
    phase_names = [cue_text+"early",cue_text+"late",cue_text+cue_names(c_i)+"omi",cue_text+cue_names(3-c_i)+"omi",];
    paired_x = repmat((1:length(phase_names))',1,8);
    for t_i = 1:2 % loop across decel vs diff_vel
        this_test_data = test_data.(cue_text)(:,:,t_i)';
        [p,tbl,stats] = testf(this_test_data,[],"off");
        c = multcompare(stats,display="off");
        multi_p = c(:,6)';
        asters = common_functions.p_to_asterisk(multi_p);
        test_result.(t_tags(t_i)) = tbl;
        this_test_data_tag = categorical([repmat(phase_names(1),8,1);repmat(phase_names(2),8,1);repmat(phase_names(3),8,1);repmat(phase_names(4),8,1)],phase_names);
        boxchart(axs(2*c_i-2+t_i),this_test_data_tag,reshape(this_test_data,[],1),MarkerStyle="none");
        plot(axs(2*c_i-2+t_i),paired_x,this_test_data','ko-');
        common_functions.add_significance_to_ax(axs(2*c_i-2+t_i),{[1.1,1.9],[1.1,1.9]+1,[1.1,1.9]+2,[1.1,2.9],[2.1,3.9],[1.1,3.9]},[1,1,1,1.05,1.1,1.15],string(multi_p([1,4,6,2,5,3])));
        title(axs(2*c_i-2+t_i),strjoin([cue_text,t_texts(t_i)]," ")+" "+test_name+"="+p)
    end
end
hold(axs,"off")
saveas(fig_test,[cf,'sup_fig_7A.png']);
close(fig_test)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7B
close all;clear;clc;
cf = [pwd,'\'];
cwa_data = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
main_srs = common_functions.get_main_samplerate();
vel_window = [-0.5,0.5]; % velocity window of getting acceleration etc
vel_smooth_info = {"sgolay",1/4};
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
cue_names = ["LED","Tone"];

% build acceleration data
vel_struct = struct;
test_data = struct;

c_i=1;
cue_text = "rew"+c_i;
phase_names = [cue_text+"early",cue_text+"late",cue_text+cue_names(c_i)+"omi",cue_text+cue_names(3-c_i)+"omi",];
test_data.(cue_text) = nan(length(phase_names),8,2);
for m_i = 1:8
    mouse_name = mouse_names(m_i);
    this_sr = main_srs(mouse_name);
    vel_window_frame = (vel_window+1)*this_sr;
    vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
    for p_i = 1:length(phase_names)
        vel = cwa_data.(phase_names(p_i)).(mouse_name).mu(:,2);
        vel_smoothed = smoothdata(vel,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));

        vel_info = common_functions.find_acc_start_end_tbl(vel_smoothed,vel_window_frame,this_sr);
        vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_original = vel;
        vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_smoothed = vel_smoothed;
        vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_info = vel_info;

        test_data.(cue_text)(p_i,m_i,:) = [-vel_info{1,"acc_peak_value"},vel_info{1,"vel_start"}-vel_info{1,"vel_end"}];
    end
end

testf = @kruskalwallis;
test_name = "kruskalwallis";

t_tags = ["decel","diff_v"];
t_texts = ["deceleration","change in velocity"];
test_result = struct;
fig_test = figure(position=[100,100,1600,900]);
tiled_test = tiledlayout(fig_test,2,2);
axs = gobjects(1,4);
for ax_i = 1:4
    axs(ax_i) = nexttile(tiled_test,ax_i);
end
hold(axs,"on")

c_i = 1;
cue_text = "rew"+c_i;
phase_names = [cue_text+"early",cue_text+"late",cue_text+cue_names(c_i)+"omi",cue_text+cue_names(3-c_i)+"omi",];
paired_x = repmat((1:length(phase_names))',1,8);
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = nan(8,8);
output_tbl_varnames = reshape(["peak_decel_","vel_change_",]+["pre";"post";"LEDomission";"TONEomission"],1,[]);
%--------------------------------------------------- output file ---------------------------------------------------%
for t_i = 1:2 % loop across decel vs diff_vel
    this_test_data = test_data.(cue_text)(:,:,t_i)';
    [p,tbl,stats] = testf(this_test_data,[],"off");
    c = multcompare(stats,display="off");
    multi_p = c(:,6)';
    asters = common_functions.p_to_asterisk(multi_p);
    test_result.(t_tags(t_i)) = tbl;
    this_test_data_tag = categorical([repmat(phase_names(1),8,1);repmat(phase_names(2),8,1);repmat(phase_names(3),8,1);repmat(phase_names(4),8,1)],phase_names);
    boxchart(axs(2*c_i-2+t_i),this_test_data_tag,reshape(this_test_data,[],1),MarkerStyle="none");
    plot(axs(2*c_i-2+t_i),paired_x,this_test_data','ko-');
    common_functions.add_significance_to_ax(axs(2*c_i-2+t_i),{[1.1,1.9],[1.1,1.9]+1,[1.1,1.9]+2,[1.1,2.9],[2.1,3.9],[1.1,3.9]},[1,1,1,1.05,1.1,1.15],string(multi_p([1,4,6,2,5,3])));
    title(axs(2*c_i-2+t_i),strjoin([cue_text,t_texts(t_i)]," ")+" "+test_name+"="+p)
    %--------------------------------------------------- output file ---------------------------------------------------%
    output_tbl(:,[1:4]+4*(t_i-1)) = this_test_data;
    %--------------------------------------------------- output file ---------------------------------------------------%
end
hold(axs,"off")
saveas(fig_test,[cf,'sup_fig7B.png'])
close(fig_test)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7C
close all;clear;clc;
cf = [pwd,'\'];
cwa_data = load([cf,'processed_and_organized_data\ITI_licking_components_window_activity_filtered.mat']);
main_srs = common_functions.get_main_samplerate();
vel_window = [-0.5,0]; % velocity window of getting acceleration etc
vel_smooth_info = {"sgolay",1/4};
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];

cue_names = "ITI";
phase_names = ["early","late"];

% build acceleration data
vel_struct = struct;
test_data = struct;
for c_i = 1
    cue_text = cue_names;
    test_data.(cue_text) = nan(length(phase_names),8,2);
    for m_i = 1:8
        mouse_name = mouse_names(m_i);
        this_sr = main_srs(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr;
        vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        for p_i = 1:length(phase_names)
            vel = cwa_data.(phase_names(p_i)).(mouse_name).mu(:,2);
            vel_smoothed = smoothdata(vel,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));

            vel_info = common_functions.find_acc_start_end_tbl(vel_smoothed,vel_window_frame,this_sr);
            vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_original = vel;
            vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_smoothed = vel_smoothed;
            vel_struct.(cue_text).(mouse_name).(phase_names(p_i)).vel_info = vel_info;

            test_data.(cue_text)(p_i,m_i,:) = [-vel_info{1,"acc_peak_value"},vel_info{1,"vel_start"}-vel_info{1,"vel_end"}];
        end
    end
end

testf = @ranksum;
test_name = "ranksum";
t_tags = ["decel","diff_v"];
t_texts = ["peak deceleration amplitude","velocity changes"];
t_texts_suf = [" (m/s)"," (m/s^2)"];
test_result = struct;
fig_test = figure(position=[100,100,1600,900]);
tiled_test = tiledlayout(fig_test,2,2);
axs = gobjects(1,4);
for ax_i = 1:4
    axs(ax_i) = nexttile(tiled_test,ax_i);
end
hold(axs,"on")
for c_i = 1
    cue_text = cue_names;
    paired_x = repmat((1:length(phase_names))',1,8);
    for t_i = 1:2 % loop across decel vs diff_vel
        this_test_data = test_data.(cue_text)(:,:,t_i)';
        [p,tbl,stats] = testf(this_test_data(:,1),this_test_data(:,2));
        asters = common_functions.p_to_asterisk(p);
        test_result.(t_tags(t_i)) = tbl;
        this_test_data_tag = categorical([repmat(phase_names(1),8,1);repmat(phase_names(2),8,1)],phase_names);
        boxchart(axs(2*c_i-2+t_i),this_test_data_tag,reshape(this_test_data,[],1),MarkerStyle='o');
        plot(axs(2*c_i-2+t_i),paired_x,this_test_data','ko-');
        common_functions.add_significance_to_ax(axs(2*c_i-2+t_i),{[1.1,1.9]},1.05,asters);
        ylabel(axs(2*c_i-2+t_i),t_texts(t_i)+t_texts_suf(t_i));
        title(axs(2*c_i-2+t_i),strjoin([cue_text,t_texts(t_i)]," ")+" "+test_name+"="+p)
    end
end
hold(axs,"off")
saveas(fig_test,[cf,'sup_fig7C.png'])
close(fig_test)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7D
close all;clear;clc;
cf = [pwd,'\'];
results = load([cf,'processed_and_organized_data\decel_signal_iti_highpass03.mat']);
results_bstrp = load([cf,'processed_and_organized_data\decel_signal_iti_btsrp_highpass03']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
mouse_names = string(fields(results)');
mouse_names = mouse_names(contains(mouse_names,"G"));
act_accvel_tag = ["acc","vel"];

% ---------------------- set up special manual plot figure first ---------------------- %
pointer_special_plot = 1;
trace_imagesc_plot = ["G21",27,18;"G21",28,18;];

trace_ylim = [-0.03,0.02;-0.01,0.02;];
imagesc_cmap = [0.03;0.02;];

nr=2;nc=size(trace_imagesc_plot,1);
manual_fig = figure(Position=[41,41,400*nc,800]);
manual_tiled = tiledlayout(manual_fig,nr,nc,TileSpacing="tight");
manual_axs = gobjects(1,nr*nc);
for i=1:nr*nc;manual_axs(i)=nexttile(manual_tiled,i);end
% ---------------------- set up special manual plot figure first ---------------------- %


circle_ct_table = {};
circle_value = [];
min_consecutive = 3;

fiber_stats_table = table(Size=[length(mouse_names),4],VariableTypes=["int16","int16","int16","int16"],VariableNames=["positive_only","negative_only","both","x"],RowNames=mouse_names);
decel_stats_table = table(Size=[length(mouse_names),2],VariableTypes=["double","double"],VariableNames=["peak_decel","vel_change",],RowNames=mouse_names);

vel_dec_fig = figure(Position=[1,41,2560,1323]);
vel_dec_tiled = tiledlayout(2,length(mouse_names),TileSpacing="tight");
vel_dec_axs = gobjects(1,2*length(mouse_names));
for i=1:2*length(mouse_names);vel_dec_axs(i)=nexttile(vel_dec_tiled,i);end
ylabel(vel_dec_axs(1),"vel"); ylabel(vel_dec_axs(1+length(mouse_names)),"dec");
hold(vel_dec_axs,"on");

for mi = 1:length(mouse_names)
    mname = mouse_names(mi);
    plot_x = (-18:18)/18;
    this_table = ct_table(ct_table.mouse_name==mname & ct_table.significance==1,:);
    this_data = results.(mname).all.combined_decel;
    this_bstrp = results_bstrp.(mname).all.bstrp;
    n_trial = size(this_data,3);

    subselect_params = {"decel_peak",1.5,"vel_before",0.6,"vel_after",0.4,"before_after_window",[0.3,0.6]};
    [valid_bit,valid_bits] = subselect_velocity(this_data(:,1:2,:),subselect_params{:});
    subselect_info = sprintf("n=%d <- before(%d),after(%d),peak(%d)",sum(valid_bit),sum(valid_bits(1,:)),sum(valid_bits(2,:)),sum(valid_bits(3,:)));
    fprintf("%s\n",subselect_info);

    % plot vel_dec trace and save to table
    plot_functions.plot_data_single(vel_dec_axs(mi),plot_x,this_data(:,1,valid_bit),plot_color=lines(1)*0.8,sem_color=lines(1));
    plot_functions.plot_data_single(vel_dec_axs(mi+length(mouse_names)),plot_x,this_data(:,2,valid_bit),plot_color=lines(1)*0.8,sem_color=lines(1));
    xline(vel_dec_axs(mi),0);yline(vel_dec_axs(mi),0);xline(vel_dec_axs(mi++length(mouse_names)),0);yline(vel_dec_axs(mi+length(mouse_names)),0);
    title(vel_dec_axs(mi+length(mouse_names)),mname);

    vel_smoothed = mean(this_data(:,2,valid_bit),3,"omitmissing");
    this_sr = 18; vel_window_frame = round(([-0.5,0.5]+1)*this_sr); vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
    vel_info = find_acc_start_end(vel_smoothed,vel_window_frame,this_sr);
    decel_stats_table{mi,:} = [-vel_info{1,"acc_peak_value"},vel_info{1,"vel_start"}-vel_info{1,"vel_end"}];
    
    % --------------------- special plots ----------------------- %
    for session_special_plot = pointer_special_plot
    this_roi_info = trace_imagesc_plot(trace_imagesc_plot(:,1)==mname,:);
    for ri=1:size(this_roi_info,1)
        this_axs = manual_axs([pointer_special_plot+ri-1,pointer_special_plot+ri-1+size(trace_imagesc_plot,1)]);
        this_limit_id = find(all(trace_imagesc_plot == this_roi_info(ri,:),2));
        this_sr = str2double(this_roi_info(ri,3));
        this_r = str2double(this_roi_info(ri,2));
        this_coord = ct_table{ct_table.mouse_name==this_roi_info(ri,1) & ct_table.ROI_original==this_r,["fiber_bottom_AP","fiber_bottom_ML","fiber_bottom_DV"]};
        title(this_axs(1),sprintf("AP/ML/DV = %0.2f/%0.2f/%0.2f",this_coord(1),this_coord(2),this_coord(3)));
        title(this_axs(2),this_roi_info(ri,1)+" ROI"+this_roi_info(ri,2));

        this_plot_x = (-this_sr:this_sr)/this_sr;
        hold(this_axs,"on")
        plot_functions.plot_data_single(this_axs(1),this_plot_x,this_data(:,this_r+2,valid_bit),plot_color=[0,0.447,0.741],sem_color=[0,0.447,0.741]);
        bstrp_plot = permute(this_bstrp(:,this_r,:),[1,3,2]);
        plot(this_axs(1),this_plot_x,bstrp_plot,'k-',LineWidth=1.5);
        xline(this_axs(1),0); yline(this_axs(1),0);
        ylim(this_axs(1),trace_ylim(this_limit_id,:));
        
        %%%%%%
        % tmp = permute(this_data(:,this_r+2,valid_bit),[3,1,2]);
        % tmp = prctile(abs(tmp),95,"all");
        %%%%%%
        tmp =imagesc_cmap(this_limit_id);
        %%%%%%
        imagesc(this_axs(2),this_plot_x,1:sum(valid_bit),permute(this_data(:,this_r+2,valid_bit),[3,1,2]),[-1,1]*tmp);
        this_axs(2).YDir = "reverse";
        colormap(this_axs(2),redblue);
        colorbar(this_axs(2));
        xline(this_axs(2),0,'k-',LineWidth=2);
        xlim(this_axs(2),[-1,1]); ylim(this_axs(2),[1,sum(valid_bit)]);

        hold(this_axs,"off")
        xlabel(this_axs(1),"time since peak deceleration (s)");ylabel(this_axs(1),"\DeltaF/F");
    end
    pointer_special_plot=pointer_special_plot+size(this_roi_info,1);
    end
    % --------------------- special plots ----------------------- %

    
    % for circle plot: just taling max,min which have 2 consecutive significant point
    circle_ct_table = cat(1,circle_ct_table,this_table);
    n_rois = size(this_table,1);
    rois = this_table.ROI_original;

    this_data = results.(mname).all.combined_decel;
    this_bstrp = results_bstrp.(mname).all.bstrp;

    this_values = nan(n_rois,2);
    for ri = 1:n_rois
        r = rois(ri);
        this_mu = mean(this_data(:,r+2,valid_bit),3,"omitmissing");
        this_control = permute(this_bstrp(:,r,:),[1,3,2]);
        tmp = this_mu<this_control(:,1);
        consecutive_bit = get_consecutive_bit(tmp,min_consecutive);
        if any(consecutive_bit)
            this_values(ri,1) = min(this_mu(consecutive_bit));
        end
        tmp = this_mu>this_control(:,2);
        consecutive_bit = get_consecutive_bit(tmp,min_consecutive);
        if any(consecutive_bit)
            this_values(ri,2) = max(this_mu(consecutive_bit));
        end
    end
    % save fiber stats data to table
    tmp = [sum(this_values(:,2)>0),sum(this_values(:,1)<0),sum(this_values(:,2)>0&this_values(:,1)<0),sum(all(isnan(this_values),2))];
    tmp(1:2) = tmp(1:2)-tmp(3);
    fiber_stats_table{mi,:} = tmp;
    circle_value = cat(1,circle_value,this_values);
end
save([cf,'processed_and_organized_data\decel_circle_min_max.mat'],"circle_value");
writetable(fiber_stats_table,[cf,'tabular_data_of_figures_autogenerated.xlsx'],WriteRowNames=true,Sheet="sup_fig7D",WriteMode="overwritesheet");

hold(vel_dec_axs,"off");
linkaxes(vel_dec_axs(1:length(mouse_names)),"xy"); linkaxes(vel_dec_axs([1:length(mouse_names)]+length(mouse_names)),"xy");
saveas(vel_dec_fig,[cf,'sup_fig7E_decel_decel.png'])
delete(vel_dec_fig)

% ---------------------- set up special manual plot figure first ---------------------- %
saveas(manual_fig,[cf,'sup_fig7E_decel_signal.png']);
delete(manual_fig)
% ---------------------- set up special manual plot figure first ---------------------- %

% circle plot for deceleration

fig = figure(position=[204,59,1328,1236]);
tiled = tiledlayout(fig,2,2);
sgtitle(tiled,"deceleration "+strjoin(mouse_names,' '))
axs=gobjects(1,4); for i=1:4;axs(i)=nexttile(tiled,i);end;
hold(axs,"on")
for pni=1:2
    ax1=axs(pni);ax2=axs(pni+2);
    switch pni 
        case 1
            cmapBounds = [-1,1]*0.017;
        case 2
            cmapBounds = [-1,1]*0.018;
    end
    plot_functions.scatter_3d_with_datatip(ax1,circle_value(:,pni)',circle_ct_table,cmapBounds=cmapBounds,...
        viewAngle=[0,90],sorting_tag="abs",colormapOption="redblue");
    plot_functions.scatter_3d_with_datatip(ax2,circle_value(:,pni)',circle_ct_table,cmapBounds=cmapBounds,...
        viewAngle=[-90,0],sorting_tag="abs",colormapOption="redblue");
end
hold(axs,"off")
title(axs(3),"min value");title(axs(4),"max value");
for ax=axs
    axis(ax,"vis3d");
end
saveas(fig,[cf,'sup_fig7D_decel_circle.png'])
delete(fig)

% get task signal for comparison
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
cwa_all_days = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']).cwa_raw;
cwa_all_days_fields = ["cueOn","cue1","activity"];
session_info = common_functions.get_training_info();
plot_info_1 = {"G21",-1,18;"G21",27,18;"G21",28,18;};
plot_imagesc_clim = {-3;0.03;0.02;};
plot_info_2 = ["cue1late",];
plot_info_3 = {16:20,11:15;16:20,11:15;16:20,11:15;};
plot_colors = lines(3);

fig = figure(Position=[100,100,1600,900]);
nr=1;
tiled = tiledlayout(fig,nr,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,nr*size(plot_info_1,1));
for i=1:nr*size(plot_info_1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axes([pi]);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    acc_smooth_info = {"sgolay",round(1/4*sr)};
    plot_x = (1:sr*4)/sr-1;
    
    if this_info{2}+3 ~= 2
        tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
    else
        tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
        early_fc = squeeze(tmp);
        early_fc = smoothdata(early_fc,1,acc_smooth_info{1},acc_smooth_info{2});
        early_fc = diff(early_fc,1,1)*sr;
        early_fc = cat(1,early_fc(1,:),early_fc);
        tmp = smoothdata(early_fc,1,acc_smooth_info{1},acc_smooth_info{2});
    end
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-0.5,2]);
end
hold(axes,"off");
saveas(fig,'sup_fig7D_task_signal.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7F


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7H
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
CT = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
main_sr = common_functions.get_main_samplerate();

phase_names = ["cue1early","cue1late"];
mouse_names = "G23";
plot_rois = [48,13];
vel_smooth_info = {"sgolay",1/4};

% build movement data
mov_output = struct;
mov_validate = struct;
for phase_name_i = 1:length(phase_names)
    phase_name = phase_names(phase_name_i);
    this_cwa = cwa.(phase_name);
    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end

    decel_struct = struct;
    decel_struct.var_names = ["time_start","time_acc_peak","time_end",...,
            "acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"];
    decel_struct.across = [];
    for mouse_name = mouse_names
        linv = squeeze(this_cwa.(mouse_name).activity(:,2,:));
        this_sr = main_sr(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr; vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        this_vel = smoothdata(linv,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
        acc_smooth_info = {vel_smooth_info{1},round(vel_smooth_info{2}*this_sr)};
        [this_vel_data,output_validate] = common_functions.find_acc_start_end(this_vel,vel_window_frame,this_sr,acc_smooth_info=acc_smooth_info);
        decel_struct.(mouse_name) = this_vel_data;
        decel_struct.across = cat(1,decel_struct.across,this_vel_data);
        mov_validate.(phase_name).(mouse_name) = output_validate;
    end
    mov_output.(phase_name) = decel_struct;
end

% plot single trial as imagesc ordered by movement
% --------heatmap plot-------- %
pattern_phase = ["cue1early","cue1late"];
xlim_bound = [-1,3];

for mouse_name = mouse_names
    this_sr = main_sr(mouse_name);
    plot_x = [1:size(cwa.(pattern_phase(2)).(mouse_name).mu,1)]/this_sr-1;
    acc_smooth_info = {vel_smooth_info{1},round(vel_smooth_info{2}*this_sr)};

    r_to_plot = [1,2,3,plot_rois+3];
    title_text = ["lick","lin V","acc","ROI "+string(plot_rois)];
    % special case: plot early vs late colorplot
    cmap_scale = nan(1,4);
    [tmp_a_e,acc_early_index] = sort(-mov_output.(pattern_phase(1)).(mouse_name)(:,4));
    [tmp_a_l,acc_late_index] = sort(-mov_output.(pattern_phase(2)).(mouse_name)(:,4));
    acc_early_index = acc_early_index(~isnan(tmp_a_e));
    acc_late_index = acc_late_index(~isnan(tmp_a_l));
    [tmp_v_e,vel_early_index] = sort(mov_output.(pattern_phase(1)).(mouse_name)(:,6));
    [tmp_v_l,vel_late_index] = sort(mov_output.(pattern_phase(2)).(mouse_name)(:,6));
    vel_early_index = vel_early_index(~isnan(tmp_v_e));
    vel_late_index = vel_late_index(~isnan(tmp_v_l));
    % clear("tmp_a_e","tmp_a_l","tmp_v_e","tmp_v_l");
    nr = 2; nc = 5;
    k = 0; fg = 0;
    for ri = 1:length(r_to_plot)
        if k == 0
            fg = fg + 1;
            k = k + 1;
            fig_1 = figure(Position=[100,100,1600,700]);
            tiled_1 = tiledlayout(nr,nc,Parent=fig_1,TileSpacing="tight");
            sgtitle(tiled_1,mouse_name+" colored by acc")
        end
        r = r_to_plot(ri);

        if r~=3 
            early_fc = squeeze(cwa.(pattern_phase(1)).(mouse_name).activity(:,r,:));
            late_fc = squeeze(cwa.(pattern_phase(2)).(mouse_name).activity(:,r,:));
        else % replace angular velocity with deceleration
            early_fc = squeeze(cwa.(pattern_phase(1)).(mouse_name).activity(:,2,:));
            early_fc = smoothdata(early_fc,1,acc_smooth_info{1},acc_smooth_info{2});
            early_fc = diff(early_fc,1,1)*main_sr(mouse_name);
            early_fc = cat(1,early_fc(1,:),early_fc);
            early_fc = smoothdata(early_fc,1,acc_smooth_info{1},acc_smooth_info{2});

            late_fc = squeeze(cwa.(pattern_phase(2)).(mouse_name).activity(:,2,:));
            late_fc = smoothdata(late_fc,1,acc_smooth_info{1},acc_smooth_info{2});
            late_fc = diff(late_fc,1,1)*main_sr(mouse_name);
            late_fc = cat(1,late_fc(1,:),late_fc);
            late_fc = smoothdata(late_fc,1,acc_smooth_info{1},acc_smooth_info{2});
        end

        axs_1 = [nexttile(tiled_1,k),nexttile(tiled_1,k+5)];

        hold(axs_1,"on")
        tmp_clim = nan(2,1);
        tmp = early_fc(:,acc_early_index)';
        tmp_clim(1) = max(abs(tmp),[],"all");
        ps_e_1 = imagesc(axs_1(1),XData=plot_x,CData=tmp);
        tmp = late_fc(:,acc_late_index)';
        tmp_clim(2) = max(abs(tmp),[],"all");
        ps_l_1 = imagesc(axs_1(2),tmp_clim,XData=plot_x,CData=tmp);
        hold(axs_1,"off")
        % actualy give early late same scale
        tmp_clim([1,2]) = max(tmp_clim([1,2]));
        tmp_clim = cat(2,-tmp_clim,tmp_clim);
        tmp = [axs_1];
        if ri>2
            for ax_i = 1:2
                ax = tmp(ax_i);
                colormap(ax,'redblue')
                ax.CLim = tmp_clim(ax_i,:);
            end
        end
        for ax_i = 1:2
            ax = tmp(ax_i);
            ax.YDir = 'reverse';
            cb = colorbar(ax);
            tk_label = linspace(cb.Limits(1),cb.Limits(2),5);
            tk_label_text = common_functions.decimal_to_str(tk_label,3);
            cb.Ticks = tk_label;
            cb.TickLabels = string(tk_label_text);
            xlim(ax,xlim_bound)
            xline(ax,0)
            yline(ax,0)
        end
        ylim(axs_1(1),[0,length(acc_early_index)]+0.5);
        ylim(axs_1(2),[0,length(acc_late_index)]+0.5);
        linkaxes(axs_1,"x")
        title(axs_1(2),title_text(ri));
        delete(axs_1(1)); % remove unwanted axes
        k = k+1;
        if k == nr * nc + 1 || ri == length(r_to_plot)
            k = 0;
            saveas(fig_1,[cf,'sup_fig7H_upper.png'])
            delete([fig_1])
        end
    end
end

% --------trace plot-------- %
phase_names = ["cue1late"];
mouse_names = "G23";
plot_rois = [48,13];

output_bkup = struct;

cdf_xline_quantile = [0.2,0.8];
vel_thres_quantile = struct;
for phase_name_i = 1:length(phase_names)
    phase_name = phase_names(phase_name_i);
    this_cwa = cwa.(phase_name);
    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end

    decel_struct = struct;
    decel_struct.var_names = ["time_start","time_acc_peak","time_end",...,
            "acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"];
    decel_struct.across = [];
    for mouse_name = ["G17","G19","G21","G22","G23","G24"]
        linv = squeeze(this_cwa.(mouse_name).activity(:,2,:));
        this_sr = main_sr(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr; vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        this_vel = smoothdata(linv,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
        this_vel_data = common_functions.find_acc_start_end(this_vel,vel_window_frame,this_sr);
        decel_struct.(mouse_name) = this_vel_data;
        decel_struct.across = cat(1,decel_struct.across,this_vel_data);
        vel_thres_quantile.(phase_name).(mouse_name) = quantile(this_vel_data(:,4),cdf_xline_quantile);
    end
    vel_thres_quantile.(phase_name).across = quantile(decel_struct.across(:,4),cdf_xline_quantile);

    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end
    high_low_decel_thres = vel_thres_quantile.(phase_name).across;

    % seperate high low decel trials
    seperated_struct = struct;
    seperated_struct.high_decel = struct; seperated_struct.medium_decel = struct; seperated_struct.low_decel = struct;
    for mouse_name = mouse_names
        this_decel = decel_struct.(mouse_name)(:,4);
        low_bit = this_decel >= high_low_decel_thres(2);
        medium_bit = this_decel > high_low_decel_thres(1) & this_decel < high_low_decel_thres(2);
        high_bit = this_decel <= high_low_decel_thres(1);
    
        seperated_struct.low_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,low_bit);
        seperated_struct.medium_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,medium_bit);
        seperated_struct.high_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,high_bit);
    end
    
    % plot seperated trials
    hi_me_lo_color = flipud(parula(3));
    decel_text = string(fields(seperated_struct)');
    for mouse_name = mouse_names
        this_sr = main_sr(mouse_name);
        plot_x = [1:size(seperated_struct.(decel_text(1)).(mouse_name),1)]/this_sr-1;

        r_to_plot = [1,2,3,plot_rois+3];
        title_text = ["lick","lin V","acc","ROI "+string(plot_rois)];
        nr = 1; nc = 5;
        k = 0; fg = 0;
        for ri = 1:length(r_to_plot)
            r = r_to_plot(ri);
            if k == 0
                fg = fg + 1;
                k = k + 1;
                fig = figure(Position=[100,100,1500,300]);
                tiled = tiledlayout(nr,nc,Parent=fig,TileSpacing="tight");
                sgtitle(tiled,mouse_name+" "+phase_name)
            end
            ax = nexttile(tiled,k);
            hold(ax,"on")
            for i = [1,3]
                if r~=3
                    tmp_act_hi = seperated_struct.(decel_text(i)).(mouse_name)(:,r,:);
                else
                    tmp_act_hi = squeeze(seperated_struct.(decel_text(i)).(mouse_name)(:,2,:));
                    tmp_act_hi = diff(tmp_act_hi,1,1);
                    tmp_act_hi = cat(1,tmp_act_hi(1,:),tmp_act_hi);
                    tmp_act_hi = smoothdata(tmp_act_hi,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
                end
                if ~isempty(tmp_act_hi)
                    p = common_functions.plot_data_single(ax,plot_x,tmp_act_hi,plot_color = hi_me_lo_color(i,:)*0.8,sem_color = hi_me_lo_color(i,:));
                end
            end
            hold(ax,"off")
            xline(ax,0)
            yline(ax,0)
            title(ax,title_text(ri))
            k = k+1;
            if k == nr * nc + 1 || ri == length(r_to_plot)
                k = 0;
                saveas(fig,[cf,'sup_fig7H_lower.png'])
                delete(fig)
            end
        end
    end
end



%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7I
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
main_sr = common_functions.get_main_samplerate();

phase_names = ["unpred","unpred"];
mouse_names = ["G19","G24"];
plot_rois = {[49],[23]};
vel_smooth_info = {"sgolay",1/4};

% build movement data
mov_output = struct;
mov_validate = struct;
for phase_name_i = 1:length(phase_names)
    phase_name = phase_names(phase_name_i);
    this_cwa = cwa.(phase_name);
    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end

    decel_struct = struct;
    decel_struct.var_names = ["time_start","time_acc_peak","time_end",...,
            "acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"];
    decel_struct.across = [];
    for mouse_name = mouse_names
        linv = squeeze(this_cwa.(mouse_name).activity(:,2,:));
        this_sr = main_sr(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr; vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        this_vel = smoothdata(linv,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
        acc_smooth_info = {vel_smooth_info{1},round(vel_smooth_info{2}*this_sr)};
        [this_vel_data,output_validate] = common_functions.find_acc_start_end(this_vel,vel_window_frame,this_sr,acc_smooth_info=acc_smooth_info);
        decel_struct.(mouse_name) = this_vel_data;
        decel_struct.across = cat(1,decel_struct.across,this_vel_data);
        mov_validate.(phase_name).(mouse_name) = output_validate;
    end
    mov_output.(phase_name) = decel_struct;
end

% plot single trial as imagesc ordered by movement
% --------heatmap plot-------- %
pattern_phase = ["unpred","unpred"];
xlim_bound = [-1,3];

for mouse_name_i = 1:length(mouse_names)
    mouse_name = mouse_names(mouse_name_i);
    plot_roi = plot_rois{mouse_name_i};
    this_sr = main_sr(mouse_name);
    plot_x = [1:size(cwa.(pattern_phase(2)).(mouse_name).mu,1)]/this_sr-1;
    acc_smooth_info = {vel_smooth_info{1},round(vel_smooth_info{2}*this_sr)};

    r_to_plot = [2,3,plot_roi+3];
    title_text = ["lin V","acc","ROI "+string(plot_roi)];
    % special case: plot early vs late colorplot
    cmap_scale = nan(1,4);
    [tmp_a_e,acc_early_index] = sort(-mov_output.(pattern_phase(1)).(mouse_name)(:,4));
    [tmp_a_l,acc_late_index] = sort(-mov_output.(pattern_phase(2)).(mouse_name)(:,4));
    acc_early_index = acc_early_index(~isnan(tmp_a_e));
    acc_late_index = acc_late_index(~isnan(tmp_a_l));
    [tmp_v_e,vel_early_index] = sort(mov_output.(pattern_phase(1)).(mouse_name)(:,6));
    [tmp_v_l,vel_late_index] = sort(mov_output.(pattern_phase(2)).(mouse_name)(:,6));
    vel_early_index = vel_early_index(~isnan(tmp_v_e));
    vel_late_index = vel_late_index(~isnan(tmp_v_l));
    % clear("tmp_a_e","tmp_a_l","tmp_v_e","tmp_v_l");
    nr = 2; nc = 3;
    k = 0; fg = 0;
    for ri = 1:length(r_to_plot)
        if k == 0
            fg = fg + 1;
            k = k + 1;
            fig_1 = figure(Position=[100,100,1600,700]);
            tiled_1 = tiledlayout(nr,nc,Parent=fig_1,TileSpacing="tight");
            sgtitle(tiled_1,mouse_name+" ordered by acc")
        end
        r = r_to_plot(ri);

        if r~=3 
            early_fc = squeeze(cwa.(pattern_phase(1)).(mouse_name).activity(:,r,:));
            late_fc = squeeze(cwa.(pattern_phase(2)).(mouse_name).activity(:,r,:));
        else % replace angular velocity with deceleration
            early_fc = squeeze(cwa.(pattern_phase(1)).(mouse_name).activity(:,2,:));
            early_fc = smoothdata(early_fc,1,acc_smooth_info{1},acc_smooth_info{2});
            early_fc = diff(early_fc,1,1)*main_sr(mouse_name);
            early_fc = cat(1,early_fc(1,:),early_fc);
            early_fc = smoothdata(early_fc,1,acc_smooth_info{1},acc_smooth_info{2});

            late_fc = squeeze(cwa.(pattern_phase(2)).(mouse_name).activity(:,2,:));
            late_fc = smoothdata(late_fc,1,acc_smooth_info{1},acc_smooth_info{2});
            late_fc = diff(late_fc,1,1)*main_sr(mouse_name);
            late_fc = cat(1,late_fc(1,:),late_fc);
            late_fc = smoothdata(late_fc,1,acc_smooth_info{1},acc_smooth_info{2});
        end

        axs_1 = [nexttile(tiled_1,k),nexttile(tiled_1,k+nc)];

        hold(axs_1,"on")
        tmp_clim = nan(2,1);
        tmp = early_fc(:,acc_early_index)';
        tmp_clim(1) = max(abs(tmp),[],"all");
        ps_e_1 = imagesc(axs_1(1),XData=plot_x,CData=tmp);
        tmp = late_fc(:,acc_late_index)';
        tmp_clim(2) = max(abs(tmp),[],"all");
        ps_l_1 = imagesc(axs_1(2),tmp_clim,XData=plot_x,CData=tmp);
        hold(axs_1,"off")
        % actualy give early late same scale
        tmp_clim([1,2]) = max(tmp_clim([1,2]));
        tmp_clim = cat(2,-tmp_clim,tmp_clim);
        tmp = [axs_1];
        if r>2
            for ax_i = 1:2
                ax = tmp(ax_i);
                colormap(ax,'redblue')
                ax.CLim = tmp_clim(ax_i,:);
            end
        end
        for ax_i = 1:2
            ax = tmp(ax_i);
            ax.YDir = 'reverse';
            cb = colorbar(ax);
            tk_label = linspace(cb.Limits(1),cb.Limits(2),5);
            tk_label_text = common_functions.decimal_to_str(tk_label,3);
            cb.Ticks = tk_label;
            cb.TickLabels = string(tk_label_text);
            xlim(ax,xlim_bound)
            xline(ax,0)
            yline(ax,0)
        end
        ylim(axs_1(1),[0,length(acc_early_index)]+0.5);
        ylim(axs_1(2),[0,length(acc_late_index)]+0.5);
        linkaxes(axs_1,"x")
        title(axs_1(2),title_text(ri));
        delete(axs_1(1)); % remove unwanted axes
        k = k+1;
        if k == nr * nc + 1 || ri == length(r_to_plot)
            k = 0;
            saveas(fig_1,[cf,'sup_fig7I_upper_',num2str(mouse_name_i),'.png'])
            delete([fig_1])
        end
    end
end

% --------trace plot-------- %
phase_names = ["unpred"];
mouse_names = ["G19","G24"];
plot_rois = {[49],[23]};

output_bkup = struct;

cdf_xline_quantile = [0.2,0.8];
vel_thres_quantile = struct;
for phase_name_i = 1:length(phase_names)
    phase_name = phase_names(phase_name_i);
    this_cwa = cwa.(phase_name);
    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end

    decel_struct = struct;
    decel_struct.var_names = ["time_start","time_acc_peak","time_end",...,
            "acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"];
    decel_struct.across = [];
    for mouse_name = ["G17","G19","G21","G22","G23","G24"]
        linv = squeeze(this_cwa.(mouse_name).activity(:,2,:));
        this_sr = main_sr(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr; vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        this_vel = smoothdata(linv,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
        this_vel_data = common_functions.find_acc_start_end(this_vel,vel_window_frame,this_sr);
        decel_struct.(mouse_name) = this_vel_data;
        decel_struct.across = cat(1,decel_struct.across,this_vel_data);
        vel_thres_quantile.(phase_name).(mouse_name) = quantile(this_vel_data(:,4),cdf_xline_quantile);
    end
    vel_thres_quantile.(phase_name).across = quantile(decel_struct.across(:,4),cdf_xline_quantile);

    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end
    high_low_decel_thres = vel_thres_quantile.(phase_name).across;

    % seperate high low decel trials
    seperated_struct = struct;
    seperated_struct.high_decel = struct; seperated_struct.medium_decel = struct; seperated_struct.low_decel = struct;
    for mouse_name = mouse_names
        this_decel = decel_struct.(mouse_name)(:,4);
        low_bit = this_decel >= high_low_decel_thres(2);
        medium_bit = this_decel > high_low_decel_thres(1) & this_decel < high_low_decel_thres(2);
        high_bit = this_decel <= high_low_decel_thres(1);
    
        seperated_struct.low_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,low_bit);
        seperated_struct.medium_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,medium_bit);
        seperated_struct.high_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,high_bit);
    end
    
    % plot seperated trials
    hi_me_lo_color = flipud(parula(3));
    decel_text = string(fields(seperated_struct)');
    for mouse_name_i = 1:length(mouse_names)
        mouse_name = mouse_names(mouse_name_i);
        plot_roi = plot_rois{mouse_name_i};
        this_sr = main_sr(mouse_name);
        plot_x = [1:size(seperated_struct.(decel_text(1)).(mouse_name),1)]/this_sr-1;
        r_to_plot = [2,3,plot_roi+3];
        title_text = ["lick","lin V","acc","ROI "+string(plot_roi)];
        nr = 1; nc = 3;
        k = 0; fg = 0;
        for ri = 1:length(r_to_plot)
            r = r_to_plot(ri);
            if k == 0
                fg = fg + 1;
                k = k + 1;
                fig = figure(Position=[100,100,900,300]);
                tiled = tiledlayout(nr,nc,Parent=fig,TileSpacing="tight");
                sgtitle(tiled,mouse_name+" "+phase_name)
            end
            ax = nexttile(tiled,k);
            hold(ax,"on")
            for i = [1,3]
                if r~=3
                    tmp_act_hi = seperated_struct.(decel_text(i)).(mouse_name)(:,r,:);
                else
                    tmp_act_hi = squeeze(seperated_struct.(decel_text(i)).(mouse_name)(:,2,:));
                    tmp_act_hi = diff(tmp_act_hi,1,1);
                    tmp_act_hi = cat(1,tmp_act_hi(1,:),tmp_act_hi);
                    tmp_act_hi = smoothdata(tmp_act_hi,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
                end
                if ~isempty(tmp_act_hi)
                    p = common_functions.plot_data_single(ax,plot_x,tmp_act_hi,plot_color = hi_me_lo_color(i,:)*0.8,sem_color = hi_me_lo_color(i,:));
                end
            end
            hold(ax,"off")
            xline(ax,0)
            yline(ax,0)
            title(ax,title_text(ri))
            k = k+1;
            if k == nr * nc + 1 || ri == length(r_to_plot)
                k = 0;
                saveas(fig,[cf,'sup_fig7I_lower_',num2str(mouse_name_i),'.png'])
                delete(fig)
            end
        end
    end
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7J
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
CT = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
main_sr = common_functions.get_main_samplerate();

cdf_xline_quantile = [0.2,0.8];
phase_names = ["cue1late","cue2late"];
mouse_names = ["G17","G19","G22","G21","G23","G24"];
vel_smooth_info = {"sgolay",1/4};

for phase_name_i = 1:length(phase_names)
    phase_name = phase_names(phase_name_i);
    this_cwa = cwa.(phase_name);
    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end

    decel_struct = struct;
    decel_struct.var_names = ["time_start","time_acc_peak","time_end",...,
            "acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"];
    decel_struct.across = [];
    for mouse_name = mouse_names
        linv = squeeze(this_cwa.(mouse_name).activity(:,2,:));
        this_sr = main_sr(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr; vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        this_vel = smoothdata(linv,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
        this_vel_data = common_functions.find_acc_start_end(this_vel,vel_window_frame,this_sr);
        decel_struct.(mouse_name) = this_vel_data;
        decel_struct.across = cat(1,decel_struct.across,this_vel_data);
        vel_thres_quantile.(phase_name).(mouse_name) = quantile(this_vel_data(:,4),cdf_xline_quantile);
    end
    vel_thres_quantile.(phase_name).across = quantile(decel_struct.across(:,4),cdf_xline_quantile);
    high_low_decel_thres = vel_thres_quantile.(phase_name).across;

    % seperate high low decel trials
    seperated_struct = struct;
    seperated_struct.high_decel = struct; seperated_struct.medium_decel = struct; seperated_struct.low_decel = struct;
    for mouse_name = mouse_names
        this_decel = decel_struct.(mouse_name)(:,4);
        low_bit = this_decel >= high_low_decel_thres(2);
        medium_bit = this_decel > high_low_decel_thres(1) & this_decel < high_low_decel_thres(2);
        high_bit = this_decel <= high_low_decel_thres(1);
    
        seperated_struct.low_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,low_bit);
        seperated_struct.medium_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,medium_bit);
        seperated_struct.high_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,high_bit);
    end

    fake_components_window = struct;
    for mouse_name = mouse_names
        this_sr = main_sr(mouse_name);
        n_rois = size(this_cwa.(mouse_name).mu,2)-3;
        % get low decel
        tmp_act_lo = seperated_struct.low_decel.(mouse_name);
        if ~isempty(tmp_act_lo)
            tmp_mu_lo = mean(tmp_act_lo,3,"omitmissing");
            tmp_sem_lo = std(tmp_act_lo,[],3,"omitmissing");
        else
            tmp_act_lo = nan(4*this_sr,n_rois+3,1);
            tmp_mu_lo = tmp_act_lo;
            tmp_sem_lo = tmp_act_lo;
        end
        % get high decel
        tmp_act_hi = seperated_struct.high_decel.(mouse_name);
        if ~isempty(tmp_act_hi)
            tmp_mu_hi = mean(tmp_act_hi,3,"omitmissing");
            tmp_sem_hi = std(tmp_act_hi,[],3,"omitmissing");
        else
            tmp_act_hi = nan(4*this_sr,n_rois+3,1);
            tmp_mu_hi = tmp_act_hi;
            tmp_sem_hi = tmp_act_hi;
        end

        if contains(phase_name,"cue1")
            fake_phase_names = "cue1" + ["early","late"];
        elseif contains(phase_name,"cue2")
            fake_phase_names = "cue2" + ["early","late"];
        elseif contains(phase_name,"unpred")
            fake_phase_names = "unpred" + ["early","late"];
        else
            cprintf("r","Input high_low_decel_thres for this phase.\n")
            keyboard
        end
        % also fake a unpred which is only used for formatting
        fake_components_window.unpred = cwa.unpred;
        % for tmp = "unpred" + ["","early","late","LEDomi","Toneomi"]
        %     fake_components_window.(tmp) = cwa.(tmp);
        % end
        

        field_windows=common_functions.get_comp_timewindow();


        fake_components_window.(fake_phase_names(1)).(mouse_name).activity = tmp_act_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).mu = tmp_mu_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).sem = tmp_sem_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).null = this_cwa.(mouse_name).null;
    
        fake_components_window.(fake_phase_names(2)).(mouse_name).activity = tmp_act_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).mu = tmp_mu_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).sem = tmp_sem_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).null = this_cwa.(mouse_name).null;
    end

    [~,tri_avg_single] = generate_tas_from_cwa(fake_components_window,CT,...
        field_windows=field_windows,...
        mu_std_factor=3);

    output_bkup.(phase_name).fake_cwa = fake_components_window;
    output_bkup.(phase_name).fake_tas = tri_avg_single;

    % stats across ROIs
    comp_names = ["pk","dp","re"];
    fig = figure(position=[100,100,1800,600]);
    tiled = tiledlayout(fig,1,3,TileSpacing="tight");
    sgtitle(tiled,phase_name+" component amplitudes low vs high decel across ROIs")
    axs = gobjects(1,3);
    for axi=1:3
        axs(axi) = nexttile(tiled,axi);
    end
    hold(axs,"on")
    [slow_data,this_table] = combine_across_mouse(tri_avg_single.(fake_phase_names(1)),CT_table=CT);
    fast_data = combine_across_mouse(tri_avg_single.(fake_phase_names(2)));
    slow_data = slow_data(:,logical(this_table{:,"significance"}));
    fast_data = fast_data(:,logical(this_table{:,"significance"}));
    for ci=1:3
        ax = axs(ci);
        %%%%%%
        [p,~,stats] = ranksum(slow_data(ci,:)',fast_data(ci,:)');
        asters = random_functions.p_to_asterisk(p);
        %%%%%%
        % [p,~,stats] = anova1([slow_data(ci,:);fast_data(ci,:)]');
        % asters = random_functions.p_to_asterisk(p);
        %%%%%%
        % [p,~,stats] = anova1([slow_data(ci,:);fast_data(ci,:)]');
        % c = multcompare(stats,display="off");
        % asters = random_functions.p_to_asterisk(c(:,6)');
        %%%%%%
        

        this_test_data_tag = categorical([repmat("low vel",size(slow_data,2),1);repmat("high vel",size(fast_data,2),1)],["low vel","high vel"]);
        boxchart(ax,this_test_data_tag,[slow_data(ci,:),fast_data(ci,:)]',MarkerStyle='none');
        swarmchart(ax,repmat([1,2],size(slow_data,2),1),[slow_data(ci,:);fast_data(ci,:)]',10,'k','filled','o',XJitterWidth=0.4);
        random_functions.add_significance_to_ax(ax,{[1.1,1.9]},[1.05],asters);
        title(ax,comp_names(ci)+" p="+p+" n_{low}="+size(slow_data,2)+" n_{high}="+size(fast_data,2));
    end
    hold(axs,"off")
    saveas(fig,[cf,'sup_fig7J_',char(phase_name),'.png'])
    delete(fig)
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7K
% top - reward
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
CT = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
main_sr = common_functions.get_main_samplerate();

cdf_xline_quantile = [0.2,0.8];
phase_names = ["unpred"];
mouse_names = ["G17","G19","G22","G21","G23","G24"];
vel_smooth_info = {"sgolay",1/4};

for phase_name_i = 1:length(phase_names)
    phase_name = phase_names(phase_name_i);
    this_cwa = cwa.(phase_name);
    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end

    decel_struct = struct;
    decel_struct.var_names = ["time_start","time_acc_peak","time_end",...,
            "acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"];
    decel_struct.across = [];
    for mouse_name = mouse_names
        linv = squeeze(this_cwa.(mouse_name).activity(:,2,:));
        this_sr = main_sr(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr; vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        this_vel = smoothdata(linv,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
        this_vel_data = common_functions.find_acc_start_end(this_vel,vel_window_frame,this_sr);
        decel_struct.(mouse_name) = this_vel_data;
        decel_struct.across = cat(1,decel_struct.across,this_vel_data);
        vel_thres_quantile.(phase_name).(mouse_name) = quantile(this_vel_data(:,4),cdf_xline_quantile);
    end
    vel_thres_quantile.(phase_name).across = quantile(decel_struct.across(:,4),cdf_xline_quantile);
    high_low_decel_thres = vel_thres_quantile.(phase_name).across;

    % seperate high low decel trials
    seperated_struct = struct;
    seperated_struct.high_decel = struct; seperated_struct.medium_decel = struct; seperated_struct.low_decel = struct;
    for mouse_name = mouse_names
        this_decel = decel_struct.(mouse_name)(:,4);
        low_bit = this_decel >= high_low_decel_thres(2);
        medium_bit = this_decel > high_low_decel_thres(1) & this_decel < high_low_decel_thres(2);
        high_bit = this_decel <= high_low_decel_thres(1);
    
        seperated_struct.low_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,low_bit);
        seperated_struct.medium_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,medium_bit);
        seperated_struct.high_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,high_bit);
    end

    fake_components_window = struct;
    for mouse_name = mouse_names
        this_sr = main_sr(mouse_name);
        n_rois = size(this_cwa.(mouse_name).mu,2)-3;
        % get low decel
        tmp_act_lo = seperated_struct.low_decel.(mouse_name);
        if ~isempty(tmp_act_lo)
            tmp_mu_lo = mean(tmp_act_lo,3,"omitmissing");
            tmp_sem_lo = std(tmp_act_lo,[],3,"omitmissing");
        else
            tmp_act_lo = nan(4*this_sr,n_rois+3,1);
            tmp_mu_lo = tmp_act_lo;
            tmp_sem_lo = tmp_act_lo;
        end
        % get high decel
        tmp_act_hi = seperated_struct.high_decel.(mouse_name);
        if ~isempty(tmp_act_hi)
            tmp_mu_hi = mean(tmp_act_hi,3,"omitmissing");
            tmp_sem_hi = std(tmp_act_hi,[],3,"omitmissing");
        else
            tmp_act_hi = nan(4*this_sr,n_rois+3,1);
            tmp_mu_hi = tmp_act_hi;
            tmp_sem_hi = tmp_act_hi;
        end

        if contains(phase_name,"cue1")
            fake_phase_names = "cue1" + ["early","late"];
        elseif contains(phase_name,"cue2")
            fake_phase_names = "cue2" + ["early","late"];
        elseif contains(phase_name,"unpred")
            fake_phase_names = "unpred" + ["early","late"];
        else
            cprintf("r","Input high_low_decel_thres for this phase.\n")
            keyboard
        end
        % also fake a unpred which is only used for formatting
        fake_components_window.unpred = cwa.unpred;
        % for tmp = "unpred" + ["","early","late","LEDomi","Toneomi"]
        %     fake_components_window.(tmp) = cwa.(tmp);
        % end
        

        field_windows=common_functions.get_comp_timewindow();


        fake_components_window.(fake_phase_names(1)).(mouse_name).activity = tmp_act_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).mu = tmp_mu_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).sem = tmp_sem_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).null = this_cwa.(mouse_name).null;
    
        fake_components_window.(fake_phase_names(2)).(mouse_name).activity = tmp_act_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).mu = tmp_mu_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).sem = tmp_sem_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).null = this_cwa.(mouse_name).null;
    end

    [~,tri_avg_single] = generate_tas_from_cwa(fake_components_window,CT,...
        field_windows=field_windows,...
        mu_std_factor=3);

    output_bkup.(phase_name).fake_cwa = fake_components_window;
    output_bkup.(phase_name).fake_tas = tri_avg_single;

    % stats across ROIs
    comp_names = ["pk","dp","re"];
    fig = figure(position=[100,100,1800,600]);
    tiled = tiledlayout(fig,1,3,TileSpacing="tight");
    sgtitle(tiled,phase_name+" component amplitudes low vs high decel across ROIs")
    axs = gobjects(1,3);
    for axi=1:3
        axs(axi) = nexttile(tiled,axi);
    end
    hold(axs,"on")
    [slow_data,this_table] = combine_across_mouse(tri_avg_single.(fake_phase_names(1)),CT_table=CT);
    fast_data = combine_across_mouse(tri_avg_single.(fake_phase_names(2)));
    slow_data = slow_data(:,logical(this_table{:,"significance"}));
    fast_data = fast_data(:,logical(this_table{:,"significance"}));
    for ci=1:3
        ax = axs(ci);
        %%%%%%
        [p,~,stats] = ranksum(slow_data(ci,:)',fast_data(ci,:)');
        asters = random_functions.p_to_asterisk(p);
        %%%%%%
        % [p,~,stats] = anova1([slow_data(ci,:);fast_data(ci,:)]');
        % asters = random_functions.p_to_asterisk(p);
        %%%%%%
        % [p,~,stats] = anova1([slow_data(ci,:);fast_data(ci,:)]');
        % c = multcompare(stats,display="off");
        % asters = random_functions.p_to_asterisk(c(:,6)');
        %%%%%%
        

        this_test_data_tag = categorical([repmat("low vel",size(slow_data,2),1);repmat("high vel",size(fast_data,2),1)],["low vel","high vel"]);
        boxchart(ax,this_test_data_tag,[slow_data(ci,:),fast_data(ci,:)]',MarkerStyle='none');
        swarmchart(ax,repmat([1,2],size(slow_data,2),1),[slow_data(ci,:);fast_data(ci,:)]',10,'k','filled','o',XJitterWidth=0.4);
        random_functions.add_significance_to_ax(ax,{[1.1,1.9]},[1.05],asters);
        title(ax,comp_names(ci)+" p="+p+" n_{low}="+size(slow_data,2)+" n_{high}="+size(fast_data,2));
    end
    hold(axs,"off")
    saveas(fig,[cf,'sup_fig7K_top.png'])
    delete(fig)
end

% bottom
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\ITI_licking_components_window_activity_filtered.mat']);
CT = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
main_sr = common_functions.get_main_samplerate();

cdf_xline_quantile = [0.2,0.8];
phase_names = ["early"];
mouse_names = ["G17","G19","G22","G21","G23","G24"];
vel_smooth_info = {"sgolay",1/4};

for phase_name_i = 1:length(phase_names)
    phase_name = phase_names(phase_name_i);
    this_cwa = cwa.(phase_name);
    vel_window = [-0.5,0];

    decel_struct = struct;
    decel_struct.var_names = ["time_start","time_acc_peak","time_end",...,
            "acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"];
    decel_struct.across = [];
    for mouse_name = mouse_names
        linv = squeeze(this_cwa.(mouse_name).activity(:,2,:));
        this_sr = main_sr(mouse_name);
        vel_window_frame = (vel_window+1)*this_sr; vel_window_frame = floor(vel_window_frame(1)):ceil(vel_window_frame(2));
        this_vel = smoothdata(linv,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
        this_vel_data = common_functions.find_acc_start_end(this_vel,vel_window_frame,this_sr);
        decel_struct.(mouse_name) = this_vel_data;
        decel_struct.across = cat(1,decel_struct.across,this_vel_data);
        vel_thres_quantile.(phase_name).(mouse_name) = quantile(this_vel_data(:,4),cdf_xline_quantile);
    end
    vel_thres_quantile.(phase_name).across = quantile(decel_struct.across(:,4),cdf_xline_quantile);
    high_low_decel_thres = vel_thres_quantile.(phase_name).across;

    % seperate high low decel trials
    seperated_struct = struct;
    seperated_struct.high_decel = struct; seperated_struct.medium_decel = struct; seperated_struct.low_decel = struct;
    for mouse_name = mouse_names
        this_decel = decel_struct.(mouse_name)(:,4);
        low_bit = this_decel >= high_low_decel_thres(2);
        medium_bit = this_decel > high_low_decel_thres(1) & this_decel < high_low_decel_thres(2);
        high_bit = this_decel <= high_low_decel_thres(1);
    
        seperated_struct.low_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,low_bit);
        seperated_struct.medium_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,medium_bit);
        seperated_struct.high_decel.(mouse_name) = this_cwa.(mouse_name).activity(:,:,high_bit);
    end

    fake_components_window = struct;
    for mouse_name = mouse_names
        this_sr = main_sr(mouse_name);
        n_rois = size(this_cwa.(mouse_name).mu,2)-3;
        % get low decel
        tmp_act_lo = seperated_struct.low_decel.(mouse_name);
        if ~isempty(tmp_act_lo)
            tmp_mu_lo = mean(tmp_act_lo,3,"omitmissing");
            tmp_sem_lo = std(tmp_act_lo,[],3,"omitmissing");
        else
            tmp_act_lo = nan(4*this_sr,n_rois+3,1);
            tmp_mu_lo = tmp_act_lo;
            tmp_sem_lo = tmp_act_lo;
        end
        % get high decel
        tmp_act_hi = seperated_struct.high_decel.(mouse_name);
        if ~isempty(tmp_act_hi)
            tmp_mu_hi = mean(tmp_act_hi,3,"omitmissing");
            tmp_sem_hi = std(tmp_act_hi,[],3,"omitmissing");
        else
            tmp_act_hi = nan(4*this_sr,n_rois+3,1);
            tmp_mu_hi = tmp_act_hi;
            tmp_sem_hi = tmp_act_hi;
        end

        fake_phase_names = "cue1" + ["early","late"];
        fake_components_window.unpred = cwa.early;
        field_windows=struct;
        field_windows.unpred = common_functions.get_ITI_licking_comp_timewindow()+1;
        field_windows.(fake_phase_names(1)) = [-0.27,0.23;0,0.6;0.2,0.8]+1;
        field_windows.(fake_phase_names(2)) = [-0.27,0.23;0,0.6;0.2,0.8]+1;



        fake_components_window.(fake_phase_names(1)).(mouse_name).activity = tmp_act_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).mu = tmp_mu_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).sem = tmp_sem_lo;
        fake_components_window.(fake_phase_names(1)).(mouse_name).null = this_cwa.(mouse_name).null;
    
        fake_components_window.(fake_phase_names(2)).(mouse_name).activity = tmp_act_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).mu = tmp_mu_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).sem = tmp_sem_hi;
        fake_components_window.(fake_phase_names(2)).(mouse_name).null = this_cwa.(mouse_name).null;
    end

    tri_avg_single = build_tas_ITI_lick_version(fake_components_window,["G12","G15","G20","G25","G26","G27"]);

    output_bkup.(phase_name).fake_cwa = fake_components_window;
    output_bkup.(phase_name).fake_tas = tri_avg_single;

    % stats across ROIs
    comp_names = ["pk","dp","re"];
    fig = figure(position=[100,100,1800,600]);
    tiled = tiledlayout(fig,1,3,TileSpacing="tight");
    sgtitle(tiled,phase_name+" component amplitudes low vs high decel across ROIs")
    axs = gobjects(1,3);
    for axi=1:3
        axs(axi) = nexttile(tiled,axi);
    end
    hold(axs,"on")
    [slow_data,this_table] = combine_across_mouse(tri_avg_single.(fake_phase_names(1)),CT_table=CT);
    fast_data = combine_across_mouse(tri_avg_single.(fake_phase_names(2)));
    slow_data = slow_data(:,logical(this_table{:,"significance"}));
    fast_data = fast_data(:,logical(this_table{:,"significance"}));
    for ci=1:3
        ax = axs(ci);
        %%%%%%
        [p,~,stats] = ranksum(slow_data(ci,:)',fast_data(ci,:)');
        asters = random_functions.p_to_asterisk(p);
        %%%%%%
        % [p,~,stats] = anova1([slow_data(ci,:);fast_data(ci,:)]');
        % asters = random_functions.p_to_asterisk(p);
        %%%%%%
        % [p,~,stats] = anova1([slow_data(ci,:);fast_data(ci,:)]');
        % c = multcompare(stats,display="off");
        % asters = random_functions.p_to_asterisk(c(:,6)');
        %%%%%%
        

        this_test_data_tag = categorical([repmat("low vel",size(slow_data,2),1);repmat("high vel",size(fast_data,2),1)],["low vel","high vel"]);
        boxchart(ax,this_test_data_tag,[slow_data(ci,:),fast_data(ci,:)]',MarkerStyle='none');
        swarmchart(ax,repmat([1,2],size(slow_data,2),1),[slow_data(ci,:);fast_data(ci,:)]',10,'k','filled','o',XJitterWidth=0.4);
        random_functions.add_significance_to_ax(ax,{[1.1,1.9]},[1.05],asters);
        title(ax,comp_names(ci)+" p="+p+" n_{low}="+size(slow_data,2)+" n_{high}="+size(fast_data,2));
    end
    hold(axs,"off")
    saveas(fig,[cf,'sup_fig7K_bottom.png'])
    delete(fig)
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 7LM
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat']);
decel_minmax = load([cf,'processed_and_organized_data\decel_circle_min_max.mat']).circle_value;
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
mouse_names = ["G17","G19","G22","G21","G23","G24"];

plot_phases = ["cue1late","rew1late"];
circle_ct_table = {};
circle_value = struct;
for p=plot_phases
    circle_value.(p) = [];
end

for mi = 1:length(mouse_names)
    mname = mouse_names(mi);
    this_table = ct_table(ct_table.mouse_name==mname & ct_table.significance==1,:);
    circle_ct_table = cat(1,circle_ct_table,this_table);
    for pname=plot_phases
        tmp1 = tas.(pname).(mname).mu_location_value(:,this_table.significance==1,1);
        tmp2 = tas.(pname).(mname).mu_location_value_significance(:,this_table.significance==1);
        tmp1(~tmp2)=nan;
        circle_value.(pname) = cat(2,circle_value.(pname),tmp1); clear("tmp1","tmp2");
    end
end

% --------------------- plot ROI traces ----------------------- %
for session_circle_plot = []
for pname=plot_phases
    fig = figure(Position=[100,100,1600,900]);
    tiled = tiledlayout(fig,2,3,TileSpacing="tight");
    sgtitle(tiled,pname+" "+strjoin(mouse_names," "));
    axs = gobjects(2,3); for i=1:6;axs(i)=nexttile(tiled,i);end;
    hold(axs,"on")
    for comp_i = 1:3
        plot_functions.scatter_3d_with_datatip(axs(comp_i),circle_value.(pname)(comp_i,:),circle_ct_table,...
            cmapBounds=[-1,1]*0.02,viewAngle=[0,90],colormapOption="redblue",sorting_tag="abs");
        plot_functions.scatter_3d_with_datatip(axs(comp_i+3),circle_value.(pname)(comp_i,:),circle_ct_table,...
            cmapBounds=[-1,1]*0.02,viewAngle=[-90,0],colormapOption="redblue",sorting_tag="abs");
    end
    hold(axs,"off")
    for ax=axs
        axis(ax,"vis3d");
    end
    saveas(fig,[save_dir,'circle_',char(pname),'.fig']);
    saveas(fig,[save_dir,'circle_',char(pname),'.png']);
    delete(fig)
end
end
% --------------------- plot ROI traces ----------------------- %

% do other comparison
for pname=plot_phases
    task_signal_original = circle_value.(pname);
    decel_signal = decel_minmax';

    % scatter and fit line
    fig = figure(Position=[100,100,1200,800]);
    tiled = tiledlayout(fig,2,3,TileSpacing="compact");
    sgtitle(tiled,pname+" "+strjoin(mouse_names," "));
    axs = gobjects(1,6); for i=1:6;axs(i)=nexttile(tiled,i);end;
    for compi = 1:3
        scatter(axs(compi),task_signal_original(compi,:),decel_signal(1,:),50,"filled");
        scatter(axs(compi+3),task_signal_original(compi,:),decel_signal(2,:),50,"filled")

        tmp = ~isnan(task_signal_original(compi,:)) & ~isnan(decel_signal(1,:));
        if sum(tmp) >=2 
            [coef,pcoef] = corr(task_signal_original(compi,tmp)',decel_signal(1,tmp)',Type="Pearson");
            title(axs(compi),sprintf("rho=%0.3f, p=%0.3f",coef,pcoef));
        end
        tmp = ~isnan(task_signal_original(compi,:)) & ~isnan(decel_signal(2,:));
        if sum(tmp) >=2 
            [coef,pcoef] = corr(task_signal_original(compi,tmp)',decel_signal(2,tmp)',Type="Pearson");
            title(axs(compi+3),sprintf("rho=%0.3f, p=%0.3f",coef,pcoef));
        end
    end
    axis(axs,"equal");
    for ax=axs
        xline(ax,0,'k-',LineWidth=2);yline(ax,0,'k-',LineWidth=2);
        xlabel(ax,"task Fc");
        ylabel(ax,"decel Fc");
    end
    saveas(fig,[save_dir,'stats_scatter_corr_',char(pname),'.fig']);
    saveas(fig,[save_dir,'stats_scatter_corr_',char(pname),'.png']);
    delete(fig)
end



%% 
% ███████╗██╗░░░██╗███╗░░██╗░█████╗░████████╗██╗░█████╗░███╗░░██╗
% ██╔════╝██║░░░██║████╗░██║██╔══██╗╚══██╔══╝██║██╔══██╗████╗░██║
% █████╗░░██║░░░██║██╔██╗██║██║░░╚═╝░░░██║░░░██║██║░░██║██╔██╗██║
% ██╔══╝░░██║░░░██║██║╚████║██║░░██╗░░░██║░░░██║██║░░██║██║╚████║
% ██║░░░░░╚██████╔╝██║░╚███║╚█████╔╝░░░██║░░░██║╚█████╔╝██║░╚███║
% ╚═╝░░░░░░╚═════╝░╚═╝░░╚══╝░╚════╝░░░░╚═╝░░░╚═╝░╚════╝░╚═╝░░╚══╝

function [valid_bit,valid_bits] = subselect_velocity(accvel,varargin)
    ip = inputParser;
    ip.addParameter("sr",18)
    ip.addParameter("decel_peak",2)
    ip.addParameter("vel_before",0.5)
    ip.addParameter("vel_after",0.3)
    ip.addParameter("before_after_window",[0.3,0.6])
    ip.parse(varargin{:});
    for j = fields(ip.Results)'
        eval([j{1}, '= ip.Results.', j{1}, ';']);
    end
    vel = permute(accvel(:,2,:),[1,3,2]);
    acc = permute(accvel(:,1,:),[1,3,2]);
    % get vel before and after by taking before_after_window
    before_after_frame = round(before_after_window*sr);
    before_frame = fliplr((sr+1)-[before_after_frame(1):before_after_frame(2)]);
    after_frame = (sr+1)+[before_after_frame(1):before_after_frame(2)];
    % filter for before/after vel and peak acc
    before_vel = mean(vel(before_frame,:),1,"omitmissing");
    before_vel_max = max(vel(before_frame,:),[],1,"omitmissing");

    after_vel = mean(vel(after_frame,:),1,"omitmissing");
    after_vel_max = max(vel(after_frame,:),[],1,"omitmissing");

    peak_acc = acc(sr+1,:);
    
    before_cri = before_vel>=vel_before & before_vel_max >= vel(sr+1,:);
    after_cri = after_vel<=vel_after & after_vel_max <= vel(sr+1,:);
    peak_cri = peak_acc<=-decel_peak;
    valid_bits = cat(1,before_cri,after_cri,peak_cri);
    valid_bit = all(valid_bits,1);
end

function valid_bit = get_consecutive_bit(binary_in,n_threshod)
    valid_bit = false(size(binary_in));
    [L,n] = bwlabel(binary_in);
    for i=1:n
        if sum(L==i) >= n_threshod
            valid_bit(L==i) = true;
        end
    end
end

function output = find_acc_start_end(vel,window,sr) % taken from test_group_by_mouse
    % Given a velocity trace([-1,3], has to start with -1), a frame window and SR. Find some velocity related fields.
    % Note vel is a n_frame*n_trial matrix.
    % "time_start","time_acc_peak","time_end","acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"
    delta_vel_threshold = 0.08;
    n_trials = size(vel,2);
    output = nan(n_trials,8);
    % find end first then go back and find peak
    [end_value,end_frame] = min(vel(window,:),[],1); end_frame = end_frame + window(1) - 1;
    for i = 1:n_trials
        this_end_frame = end_frame(i);
        this_vel = vel(:,i);
        [start_value,start_frame] = max(this_vel(window(1):this_end_frame-1)); start_frame = start_frame + window(1) - 1;
        if (start_value - end_value(i)) < delta_vel_threshold
            warning("Diff vel too small.")
            output(i,[1,3,4,5,7,8]) = [start_frame/sr,this_end_frame/sr,0,start_value,end_value(i),mean(this_vel(start_frame:this_end_frame))];
            continue
        end

        % get acc
        acc = diff(this_vel)*sr; acc = [acc(1);acc];
        target_acc = acc(start_frame:this_end_frame);
        try
        [~, negative_parts] = separate_positive_negative(target_acc); negative_parts = negative_parts + start_frame - 1;
        catch err
            keyboard
        end
        n_decel = size(negative_parts,1);
        delta_vel = nan(1,n_decel);
        for j = 1:n_decel
            delta_vel(j) = sum(acc(negative_parts(j,1):negative_parts(j,2)))/sr;
        end
        [delta_vel_value, delta_vel_I] = min(delta_vel);
        if abs(delta_vel_value) < delta_vel_threshold
            % warning("Probably no continouse acceleration. Set deceleration as 0.")
            acc_value = 0; acc_frame = nan;vel_pk = nan;
        else
            this_vel_id = negative_parts(delta_vel_I,:);
            [acc_value,acc_frame] = min(acc(this_vel_id(1):this_vel_id(2))); acc_frame = acc_frame + this_vel_id(1) - 1;
            vel_pk = this_vel(acc_frame);
        end    
        output(i,:) = [start_frame/sr,acc_frame/sr,this_end_frame/sr,acc_value,start_value,vel_pk,end_value(i),mean(this_vel(start_frame:this_end_frame))];
    end
    output = array2table(output,VariableNames=["time_start","time_acc_peak","time_end","acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"]);

    function [positive_parts, negative_parts] = separate_positive_negative(array)
        array = reshape(array,1,[]);
        % Check if the array is empty or has only one element
        if isempty(array) || numel(array) == 1
            return;
        end
    
        % Find the indices where the sign of the elements change
        sign_changes = diff(sign(array));
    
        % Find the indices where the positive parts start and end
        positive_starts = find(sign_changes > 0) + 1;
        positive_ends = find(sign_changes < 0);
    
        % Find the indices where the negative parts start and end
        negative_starts = find(sign_changes < 0) + 1;
        negative_ends = find(sign_changes > 0);
    
        % If the array starts with a positive part, add the first positive part
        if sign(array(1)) == 1
            positive_starts = [1, positive_starts];
        end
    
        % If the array ends with a positive part, add the last positive part
        if sign(array(end)) == 1
            positive_ends = [positive_ends, numel(array)];
        end
    
        % If the array starts with a negative part, add the first negative part
        if sign(array(1)) == -1
            negative_starts = [1, negative_starts];
        end
    
        % If the array ends with a negative part, add the last negative part
        if sign(array(end)) == -1
            negative_ends = [negative_ends, numel(array)];
        end
    
        positive_parts = [positive_starts;positive_ends]';
        negative_parts = [negative_starts;negative_ends]';
    end
end

function [sig_frames,tri_avg_single] = generate_tas_from_cwa(across_struct,CT,varargin)
    field_windows = [];
    mu_std_factor = 3;
    mouse_name_to_plot = [];

    ip = inputParser;
    ip.addParameter('field_windows',[]);
    ip.addParameter('mu_std_factor',3);
    ip.addParameter('mouse_name_to_plot',[]);
    ip.parse(varargin{:})
    for j = fields(ip.Results)'
        eval([j{1},'=ip.Results.',j{1},';']);
    end

    % get max peak,min dip and max rebound from across day tri_avg and
    % auc_loci info
    early_late_cue1_cue2_0 = string(fields(across_struct)');
    early_late_cue1_cue2_reorder = contains(early_late_cue1_cue2_0,"unpred") | contains(early_late_cue1_cue2_0,"cue");
    early_late_cue1_cue2 = [early_late_cue1_cue2_0(early_late_cue1_cue2_reorder),early_late_cue1_cue2_0(~early_late_cue1_cue2_reorder)];

    sig_frames = struct;
    std_profile = struct;
    tri_avg_single = struct;
    
    for elc1c2 = early_late_cue1_cue2
        mouse_names = string(fields(across_struct.(elc1c2))');
        for mouse_name = mouse_names
            if isempty(mouse_name_to_plot)
                mouse_name_to_plot = mouse_names;
            end
            fr = get_frequency_mouse(mouse_name);
            n_frame = size(across_struct.(elc1c2).(mouse_name).activity,1);
            n_rois = size(across_struct.(elc1c2).(mouse_name).activity,2)-3;
            n_trial = size(across_struct.(elc1c2).(mouse_name).activity,3);
            mu = across_struct.(elc1c2).(mouse_name).mu;
            activity = across_struct.(elc1c2).(mouse_name).activity;

            count_time_pos = zeros(n_frame,n_rois);
            count_time_neg = zeros(n_frame,n_rois);
            count_time = zeros(n_frame,n_rois);
            mu_stds = nan(2,n_rois);

            this_CT = CT(CT{:,"mouse_name"}==mouse_name,:);
            sig_rois = find(this_CT{:,"significance"}');

            single_values = nan(3,n_rois,n_trial);
            mu_loc = nan(3,n_rois,2);
            mu_significance = zeros(3,n_rois);
            half_max_dur = nan(3,n_rois,5); % comp X rois X [half_max_duration, start_s, end_s, start_frame, end_frame]
            half_max_dur_single = nan(3,n_rois,3,n_trial); % just keep halfmax duration would be enough. It would be too hard or unnecessary (I guess) to validate for each single trials
            for r = sig_rois
                r3 = r + 3;
                activity_r = activity(:,r3,:);
                activity_r = permute(activity_r,[1,3,2]);
                mu_r = mu(:,r3);
                if contains(elc1c2,"unpred") || contains(elc1c2,"cue")
                    mu_std_0 = std(mu_r(1:round(fr)));
                    mu_std_1 = across_struct.(elc1c2).(mouse_name).null(r3);
                else
                    switch elc1c2
                        case "rew1early"
                            std_component_str = "cue1early";
                        case "rew1late"
                            std_component_str = "cue1late";
                        case "rew2early"
                            std_component_str = "cue2early";
                        case "rew2late"
                            std_component_str = "cue2late";
                    end
                    mu_std_0 = std_profile.(std_component_str).(mouse_name)(1,r);
                    mu_std_1 = std_profile.(std_component_str).(mouse_name)(2,r);
                end
                mu_stds(1,r) = mu_std_0;
                mu_stds(2,r) = mu_std_1;

                mu_std = mu_std_factor*mu_std_1;

                sig_frame_pos = mu_r > mu_std;
                sig_frame_neg = mu_r < -mu_std;

                sig_frame = sig_frame_pos | sig_frame_neg;
                count_time(:,r) = sig_frame;
                count_time_pos(:,r) = sig_frame_pos;
                count_time_neg(:,r) = sig_frame_neg;

                if ~isempty(field_windows)
                    if contains(elc1c2,"unpred")
                        window = round(field_windows.unpred*fr);
                    else
                        window = round(field_windows.(elc1c2)*fr);
                    end
                    
                    x1 = window(1,1):window(1,2);
                    x2 = window(2,1):window(2,2);
                    x3 = window(3,1):window(3,2);

                    [M1,I1] = max(mu_r(x1)); I1 = I1 + x1(1) - 1;
                    [M2,I2] = min(mu_r(x2)); I2 = I2 + x2(1) - 1;
                    [M3,I3] = max(mu_r(x3)); I3 = I3 + x3(1) - 1;

                    dip_recovery_window = fr;
                    dip_recovery_factor = 2/3;
                    after_I2_end = min([(I2 + dip_recovery_window),length(mu_r)]);
                    after_I2 = [(I2 + 1):after_I2_end];
                    [M2_1,I2_1] = max(mu_r(after_I2)); I2_1 = I2_1 + I2;
                    
                    mu_loc(1,r,1) = M1; mu_loc(1,r,2) = I1;
                    single_values(1,r,:) = activity_r(I1,:);
                    if M1 > mu_std && sum(sig_frame_pos(I1-1:I1+1)) >= 2
                        mu_significance(1,r) = 1;
                    end
                    
                    mu_loc(2,r,1) = M2; mu_loc(2,r,2) = I2;
                    single_values(2,r,:) = activity_r(I2,:);
                    if M2 < -mu_std && sum(sig_frame_neg((I2-1:I2+1))) >= 2 && (M2_1 > M2 * dip_recovery_factor || M2_1 > -mu_std)
                        mu_significance(2,r) = 1;
                    end

                    if I3 > I2 || isnan(I2) % rebound after dip
                        mu_loc(3,r,1) = M3; mu_loc(3,r,2) = I3;
                        single_values(3,r,:) = activity_r(I3,:);
                        if M3 > mu_std && sum(sig_frame_pos((I3-1:I3+1))) >= 2
                            mu_significance(3,r) = 1;
                        end
                    end
                    % merge window3 and window1 if there is a slow peak
                    if ~isnan(mu_loc(1,r,1)) && ~isnan(mu_loc(3,r,1)) && all(mu_r(I1:I3) > mu_std)
                        [M4,I4] = max(mu_r(I1:I3)); I4 = I4 + I1 - 1;

                        mu_loc(1,r,1) = M4; mu_loc(1,r,2) = I4;
                        single_values(1,r,:) = activity_r(I4,:);

                        mu_loc(3,r,1) = nan; mu_loc(3,r,2) = nan;
                        single_values(3,r,:) = nan(size(single_values(3,r,:)));
                        mu_significance(3,r) = 0;
                    end


                    % get half_max durations (note the tri-avg duration is not reflection 
                    % single trial level activity because single trials are not perfectly aligned)
                    Ms = reshape(mu_loc(1:3,r,1),1,[]); Is = reshape(mu_loc(1:3,r,2),1,[]);
                    for comp_i = 1:3
                        if ~isnan(Ms(comp_i))
                            this_I = Is(comp_i);
                            this_M = Ms(comp_i) * (-1)^(comp_i-1);
                            this_fc = mu_r * (-1)^(comp_i-1);
                            this_act = activity_r * (-1)^(comp_i-1);

                            % calculate intersect and get half_max duration
                            [x_x_l,x_y_l,x_x_r,x_y_r] = get_half_max_duration(this_fc,this_I,this_M);
                            half_max_dur(comp_i,r,[1,2,3,4,5]) = [x_x_r-x_x_l,x_x_l,x_x_r,x_x_l,x_x_r];
                            %       half_max dur for single trials (use same amplitude and time as in average level)
                            for t_i = 1:n_trial
                                [x_x_l,x_y_l,x_x_r,x_y_r] = get_half_max_duration(this_act(:,t_i),this_I,this_M);
                                half_max_dur_single(comp_i,r,[2,3],t_i) =[x_x_l,x_x_r];
                            end
                            half_max_dur_single(comp_i,r,1,:) = half_max_dur_single(comp_i,r,3,:) - half_max_dur_single(comp_i,r,2,:);
                        end
                    end
                end
            end % loop thru rois
            sig_frames.(elc1c2).(mouse_name).total = count_time;
            sig_frames.(elc1c2).(mouse_name).pos = count_time_pos;
            sig_frames.(elc1c2).(mouse_name).neg = count_time_neg;
            std_profile.(elc1c2).(mouse_name) = mu_stds;
            if ~isempty(field_windows)
                tri_avg_single.(elc1c2).(mouse_name).single_values = single_values;
                tri_avg_single.(elc1c2).(mouse_name).mu_location_value = mu_loc;
                tri_avg_single.(elc1c2).(mouse_name).mu_location_value(:,:,2) = tri_avg_single.(elc1c2).(mouse_name).mu_location_value(:,:,2)/fr;
                tri_avg_single.(elc1c2).(mouse_name).mu_location_value_significance = mu_significance;
                half_max_dur(:,:,[1,2,3]) = half_max_dur(:,:,[1,2,3])/fr;
                half_max_dur_single = half_max_dur_single/fr;
                tri_avg_single.(elc1c2).(mouse_name).half_max_dur = half_max_dur;
                tri_avg_single.(elc1c2).(mouse_name).half_max_dur_single = half_max_dur_single;
            end
        end % loop thru mouse
    end % loop thru phases (early late cue1 cue2)
    tri_avg_single = tri_avg_single_class.add_diff_placeholder(tri_avg_single);

    % helper functions
    function [x_x_l,x_y_l,x_x_r,x_y_r] = get_half_max_duration(this_fc,this_I,this_M)
        % helper function to get halfmax duration of avg and single trials

        % if trace is even below threshold, just return nan
        if this_fc(this_I) <= this_M/2
            [x_x_l,x_y_l,x_x_r,x_y_r] = deal(nan);
            return
        end

        hm_l = find(this_fc(1:this_I) <= this_M/2, 1, "last");
        hm_r = find(this_fc(this_I:end) <= this_M/2, 1, "first"); hm_r = hm_r + this_I - 1;
        if ~isempty(hm_l)
            ls_x = [hm_l,hm_l+1];
            ls_y1 = this_fc(ls_x);
            ls_y2 = [this_M/2,this_M/2];
            [x_x_l,x_y_l,on_both_line_l] = random_functions.line_intersect(ls_x,ls_y1,ls_x,ls_y2);
            if ~on_both_line_l
                [x_x_l,x_y_l] = deal(nan);
            end
        else
            [x_x_l,x_y_l] = deal(nan);
        end

        if ~isempty(hm_r)
            ls_x = [hm_r-1,hm_r];
            ls_y1 = this_fc(ls_x);
            ls_y2 = [this_M/2,this_M/2];
            [x_x_r,x_y_r,on_both_line_r] = random_functions.line_intersect(ls_x,ls_y1,ls_x,ls_y2);
            if ~on_both_line_r
                [x_x_r,x_y_r] = deal(nan);
            end
        else % this is more likely to happen than left
            [x_x_r,x_y_r] = deal(nan);
        end

        % % test purpose
        % if ~isempty(hm_l) && ~isempty(hm_r) && (~on_both_line_l || ~on_both_line_r)
        %     text_fig = figure();
        %     plot(this_fc)
        %     hold on
        %     line([hm_l,hm_r],ls_y2)
        %     xline([hm_l,hm_r])
        %     keyboard
        %     close(text_fig)
        % end
    end

    function tas_out = add_diff_placeholder(tas_in)
        % add to tri_avg_single placeholders of differences to be used in plot app
        tag_1s = ["unpred","LED","Tone","rewLED","rewTone"];
        tag_2s = ["_late_early_diff","_LED_omission_late_diff","_Tone_omission_late_diff"];
        tag_ori_1s = ["unpred","cue1","cue2","rew1","rew2"];
        tag_ori_2s = ["early","late";"late","LEDomi";"late","Toneomi"];
    
        for i1 = 1:length(tag_1s)
            for i2 = 1:length(tag_2s)
                fnames = tag_ori_1s(i1)+tag_ori_2s(i2,:);
                if ~isfield(tas_in,fnames(1)) || ~isfield(tas_in,fnames(2))
                    continue
                end
                included_names = intersect(string(fields(tas_in.(fnames(1)))'),string(fields(tas_in.(fnames(2)))'));
                mouse_name_place = cat(1,included_names,cellstr(repmat("placeholder",size(included_names))));
                mouse_name_place = reshape(mouse_name_place,1,[]);
                tmp = struct(mouse_name_place{:});
                tmp = structfun(@string,tmp,UniformOutput=0);
                tas_out.(tag_1s(i1)+tag_2s(i2)) = tmp;
            end
        end
    
        f_names = string(fields(tas_in)');
        for f_n = f_names
            tas_out.(f_n) = tas_in.(f_n);
        end
    end

    function fr = get_frequency_mouse(mouse_name)
        % get the main sample rates of the input mouse
        Hz11 = ["G09"];
        Hz18 = ["G17","G19","G20","G22","G21","G23","G24","G25","G26","G27","UG20","UG21","UG22"];
        Hz30 = ["G11","G12","G13","G15"];
        if contains(mouse_name,Hz11)
            fr = 11;
        elseif contains(mouse_name,Hz18)
            fr = 18;
        elseif contains(mouse_name,Hz30)
            fr = 30;
        else
            error("IfElseException: mouse name not included.")
        end
    end

end

function [tas_out,ct_out] = combine_across_mouse(tas,varargin)
    % combining single tas of single mouse, setting insig as nans
    set_insig_as_nan = 0;
    CT_table = {};

    ip = inputParser;
    ip.addParameter("set_insig_as_nan",0)
    ip.addParameter("CT_table",{})
    ip.parse(varargin{:})
    for j = fields(ip.Results)'
        eval([j{1},'=ip.Results.',j{1},';']);
    end

    tas_out = [];
    ct_out = {};
    mouse_names = string(fields(tas)');
    for m_name = mouse_names
        tmp1 = tas.(m_name).mu_location_value(:,:,1);
        tmp2 = tas.(m_name).mu_location_value_significance;
        if set_insig_as_nan
            tmp1(~logical(tmp2)) = nan;
        end
        tas_out = cat(2,tas_out,tmp1);

        if ~isempty(CT_table)
            ct_out = cat(1,ct_out,CT_table(CT_table{:,"mouse_name"}==m_name,:));
        end
    end
end
