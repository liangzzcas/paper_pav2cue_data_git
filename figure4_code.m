%% figure 4B
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = "G23";
lick_data = load([cf,'processed_and_organized_data\across_mice_lick_index_data_whole_ITI.mat']);
training_info = common_functions.get_training_info();
plot_color = lines(2);

lgs = gobjects(1,2);
fig = figure(position=[100,100,800,600]);
ax = axes(fig);
hold(ax,"on")
for mouse_name = mouse_names
    last_day = training_info{[training_info{:,1}]==mouse_name,4};

    plot_data = lick_data.(mouse_name).single_trial_struct.across_1s.cue1_index_across(:,1:last_day);
    % datamu = mean(plot_data,1,"omitmissing");
    % datasem = std(plot_data,[],1,"omitmissing")/sqrt(size(plot_data,1));
    lgs(1)=common_functions.plot_data_single(ax,1:last_day,plot_data',plot_color=plot_color(1,:)*0.8,sem_color=plot_color(1,:),Marker='s');
    plot_data = lick_data.(mouse_name).single_trial_struct.across_1s.cue2_index_across(:,1:last_day);
    % datamu = mean(plot_data,1,"omitmissing");
    % datasem = std(plot_data,[],1,"omitmissing")/sqrt(size(plot_data,1));
    lgs(2)=common_functions.plot_data_single(ax,1:last_day,plot_data',plot_color=plot_color(2,:)*0.8,sem_color=plot_color(2,:),Marker='s');
end
hold(ax,"on")
xline(ax,[training_info{[training_info{:,1}]==mouse_name,2:3}]+1)
legend(ax,lgs,["Light","Tone"])
saveas(fig,[cf,'fig4B.png'],'png')
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4C
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
cwa_all_days = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']).cwa_raw;
cwa_all_days_fields = ["cueOn","cue1","activity"];
session_info = common_functions.get_training_info();
plot_info_1 = {"G23",1,18;};
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
    for i=1:3
        ylabel(this_axes(i),plot_info_2(i));
    end
end
hold(axes,"off");
saveas(fig,'fig4C.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4EG
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
cwa_all_days = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']).cwa_raw;
cwa_all_days_fields = ["cueOn","cue1","activity"];
session_info = common_functions.get_training_info();
% plot_info_1 = {"G23",48,18;"G23",36,18;"G21",41,18;"G21",9,18;};
plot_info_1 = {"G23",48,18;"G23",36,18;};
plot_imagesc_clim = {0.03;0.05;};
plot_info_2 = ["cue1late","cue1LEDomi","cue1Toneomi"];
plot_info_3 = {11:15,16:20;11:15,16:20;};
plot_colors = lines(3);

fig = figure(Position=[100,100,1600,900]);
nr=4;
tiled = tiledlayout(fig,nr,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,nr*size(plot_info_1,1));
for i=1:nr*size(plot_info_1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axes([pi,pi+size(plot_info_1,1),pi+2*size(plot_info_1,1),pi+3*size(plot_info_1,1)]);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    
    plot_x = (1:sr*4)/sr-1;
    tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
    
    imagesc(this_axes(2),plot_x,1:nt1,permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
    this_axes(2).YDir = "reverse";
    colormap(this_axes(2),common_functions.redblue());
    colorbar(this_axes(2));
    



    
    tmp = cwa.(plot_info_2(2)).(this_info{1}).activity(:,this_info{2}+3,:);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));
    % get all days for omission
    tmp = arrayfun(@(x) getfield(cwa_all_days.(this_info{1}).(x),cwa_all_days_fields{:}),"file"+plot_info_3{pi,1}, UniformOutput=false);
    tmp = cell2mat(reshape(tmp,1,1,[])); tmp = tmp(:,this_info{2}+3,:);
    nt2 = size(tmp,3);
    imagesc(this_axes(3),plot_x,1:nt2,permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
    this_axes(3).YDir = "reverse";
    colormap(this_axes(3),common_functions.redblue());
    colorbar(this_axes(3));
    




    tmp = cwa.(plot_info_2(3)).(this_info{1}).activity(:,this_info{2}+3,:);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(3,:)*0.8,sem_color=plot_colors(3,:));
    % get all days for omission
    tmp = arrayfun(@(x) getfield(cwa_all_days.(this_info{1}).(x),cwa_all_days_fields{:}),"file"+plot_info_3{pi,2}, UniformOutput=false);
    tmp = cell2mat(reshape(tmp,1,1,[])); tmp = tmp(:,this_info{2}+3,:);
    nt3 = size(tmp,3);
    imagesc(this_axes(4),plot_x,1:nt3,permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
    this_axes(4).YDir = "reverse";
    colormap(this_axes(4),common_functions.redblue());
    colorbar(this_axes(4));



    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-0.5,2]);
    for ax=this_axes([2,3,4])
        xline(ax,0,'-');
        xlim(ax,plot_x([1,end]));
    end
    ylim(this_axes(2),[1,nt1]);
    ylim(this_axes(3),[1,nt2]);
    ylim(this_axes(4),[1,nt3]);
end
hold(axes,"off");
saveas(fig,'fig4EG.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4FH
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue1LEDomi","cue1late"]},mouse_names,'fig4FH',tabular_sheet_name="fig4FH");


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4K
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.gradient_test(cf,"tas_normal", {{"cue1LEDomi_cue1late",1,2;"cue1LEDomi_cue1late",3,1;},{"cue2Toneomi_cue2late",1,2;"cue2Toneomi_cue2late",3,1;},},mouse_names,...
    {{"cue1_pk1";"cue2_pk1";},{"cue1_pk2";"cue2_pk2";},},"cue_extinction",'fig4K');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4L
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.gradient_test(cf,"tas_normal", {{"cue1LEDomi_cue1late",1,2;"cue2Toneomi_cue2late",1,2;},{"cue1LEDomi_cue1late",2,1;"cue2Toneomi_cue2late",2,1;},{"cue1LEDomi_cue1late",3,1;"cue2Toneomi_cue2late",3,1;},},mouse_names,...
   {{"cue1_pk";"cue1_dip";"cue1_re";},{"cue2_pk";"cue2_dip";"cue2_re"},},"cue_extinction",'fig4L');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4M


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4MNOP
close all;clear;clc;
cf = [pwd,'\'];

include_days = common_functions.get_include_days();
comp_windows = common_functions.get_comp_timewindow();
main_srs = common_functions.get_main_samplerate();
tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
cwa_plus = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']);
lick_info = load([cf,'processed_and_organized_data\across_mice_lick_index_data_whole_ITI.mat']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

cwa_raw_day_info = cwa_plus.cwa_raw_day_info;
cwa_plus = cwa_plus.cwa_raw;


e_l_text = ["early","late"];
phases = ["late","LEDomi","Toneomi"]; phase_cues = ["late","LED","Tone"];
cue_id = {"cue1","cue2"}; cue_type = {"LED trial","Tone trial"}; cwa_phase_name_list = {"cue1"+phases,"cue2"+phases};
mouse_names = ["G12","G23"];

single_tral_window_from_ave_s = 0.25; % looking for single trial comp peak during +-0.25s from tri_avg

omitted_phase_cwas = struct;
omitted_phase_tass = struct;

for ci = 1:length(cue_id)
    omitted_phase_cwa = struct;
    omitted_phase_tas = struct;

    cwa_phase_names = cwa_phase_name_list{ci};
    getfield_names = {"cueOn",cue_id{ci},"activity"};
    for m_i = 1:length(mouse_names)
        mouse_name = mouse_names(m_i);
        main_sr = main_srs(mouse_name);
        mouse_day_info = cwa_raw_day_info(string(cwa_raw_day_info(:,1)) == mouse_name,:);
    
        for p_i = 1:length(phases)
            phase = phases(p_i);
            if contains(phase,["early","late"])
                this_day = include_days.(mouse_name+phase);
            else
                tmp_id = find(contains(string(mouse_day_info([5,6])),phase_cues(p_i)));
                tmp_day = cell2mat(mouse_day_info(2:4));
                this_day = (tmp_day(tmp_id)+1):tmp_day(tmp_id+1);
            end
            if isempty(this_day)
                continue
            end
            
            % build cwa
            this_act = [];
            this_act_num_trial = [];
            for d_i = 1:length(this_day)
                tmp_act = getfield(cwa_plus.(mouse_name).("file"+this_day(d_i)),getfield_names{:});
                this_sr = round(size(tmp_act,1)/4); % Assuming 4 second
                if this_sr ~= main_sr
                    tmp_act = common_functions.interp_ta(tmp_act,this_sr,main_sr);
                end
                tmp_n_trial = size(tmp_act,3);
                this_act = cat(3,this_act,tmp_act);
                this_act_num_trial = cat(2,this_act_num_trial,tmp_n_trial);
            end
    
            % group cwa by batch size
            omitted_phase_cwa.(cwa_phase_names(p_i)).(mouse_name).activity = this_act;
            omitted_phase_cwa.(cwa_phase_names(p_i)).(mouse_name).activity_trial_number = this_act_num_trial;
            
            % build tas based on cwa
            lick_comp_window_pk_dp_re_s = comp_windows.(cwa_phase_names(p_i));
            omitted_phase_tas.(cwa_phase_names(p_i)).(mouse_name).single_values = []; % just for pre-arrange fields order
            omitted_phase_tas.(cwa_phase_names(p_i)).(mouse_name).single_values_trial_number = this_act_num_trial;
            omitted_phase_tas.(cwa_phase_names(p_i)).(mouse_name).mu_location_value_significance = tas.(cwa_phase_names(p_i)).(mouse_name).mu_location_value_significance;
            tas_comp_windows_s = tas.(cwa_phase_names(p_i)).(mouse_name).mu_location_value(:,:,2);
            tas_comp_windows = round(tas_comp_windows_s * main_sr);
            single_tral_window = ceil(single_tral_window_from_ave_s * main_sr);
            single_tral_window = [-single_tral_window:single_tral_window];
            
            n_rois = size(this_act,2)-3;
            n_trial = size(this_act,3);
            this_single_values = nan(3,n_rois,n_trial,2);
            for r = 1:n_rois
                for c_i = 1:3
                    if isnan(tas_comp_windows(c_i,r))
                        continue
                    end
                    this_single_tral_window = tas_comp_windows(c_i,r) + single_tral_window;
                    if c_i ~= 2
                        [M,I] = max(squeeze(this_act(this_single_tral_window,r+3,:)),[],1,"omitmissing");
                    else
                        [M,I] = min(squeeze(this_act(this_single_tral_window,r+3,:)),[],1,"omitmissing");
                    end
                    I = I + this_single_tral_window(1)-1;
                    this_single_values(c_i,r,:,1) = M;
                    this_single_values(c_i,r,:,2) = I;
                end
            end
            omitted_phase_tas.(cwa_phase_names(p_i)).(mouse_name).single_values = this_single_values;
    
            % also split omission phase into early late
            if contains(phase,["late","Toneomi"]) || contains(mouse_name,["G25","G26","G27"])
                continue
            end
            e_l_day = get_early_late_days("cue1_pre_lick","LEDomi",mouse_name);
            tmp = cumsum([0,this_act_num_trial]);
            e_l_trial_id = [tmp(e_l_day{1}(1))+1, tmp(e_l_day{1}(end)+1), tmp(e_l_day{2}(1))+1, tmp(e_l_day{2}(end)+1)];
            % Get cwa for early late
            omitted_phase_cwa.(cwa_phase_names(p_i)+"_early").(mouse_name).activity = this_act(:,:,e_l_trial_id(1):e_l_trial_id(2));
            omitted_phase_cwa.(cwa_phase_names(p_i)+"_late").(mouse_name).activity = this_act(:,:,e_l_trial_id(3):e_l_trial_id(4));
            % Get tas for early late. This is (probably) only used for early late diff,
            % so the single trial value is (probably) only used for significance test. 
            % As a result, just use same frame point instead of single_tral_window.
            null_std_factor = 3;
            comp_time_window = get_comp_time_window(cwa_phase_names(p_i));
            pk_dp_re_window = round(comp_time_window * main_sr);
            pk_window = pk_dp_re_window(1,1):pk_dp_re_window(1,2);
            dp_window = pk_dp_re_window(2,1):pk_dp_re_window(2,2);
            re_window = pk_dp_re_window(3,1):pk_dp_re_window(3,2);
            for e_l_i = 1:2
                this_el = e_l_text(e_l_i);
                this_act = omitted_phase_cwa.(cwa_phase_names(p_i)+"_"+this_el).(mouse_name).activity(:,1:end,:);
                null_thres = mean(this_act(1:main_sr,:,:),3,"omitmissing");
                null_thres = (std(null_thres,[],1)*null_std_factor);
                [pk_pks,pk_locs_frame,pk_pks_single,pk_sig] = get_peak(this_act,pk_window,null_thres);
                pk_locs = pk_locs_frame/main_sr-1;
                [dp_pks,dp_locs,dp_pks_single,dp_sig] = get_peak(-this_act,dp_window,null_thres);
                dp_locs = dp_locs/main_sr-1;
                dp_pks = -dp_pks; dp_pks_single = -dp_pks_single;
                % get sig peak to be excluded when finding rebound
                tmp_pk_locs = pk_locs_frame;
                tmp_pk_locs(~logical(pk_sig)) = nan;
                [re_pks,re_locs,re_pks_single,re_sig] = get_peak(this_act,re_window,null_thres,exclude_peaks=tmp_pk_locs);
                re_locs = re_locs/main_sr-1;
                
                single_values = cat(3,pk_pks_single,dp_pks_single,re_pks_single);
                single_values = permute (single_values,[3,2,1]);
                mu_location_value = cat(3,cat(1,pk_pks,dp_pks,re_pks),cat(1,pk_locs,dp_locs,re_locs));
                mu_location_value_significance = cat(1,pk_sig,dp_sig,re_sig); % to be assigned
    
                tmp_tas.single_values = single_values;
                tmp_tas.mu_location_value = mu_location_value;
                tmp_tas.mu_location_value_significance = mu_location_value_significance;
                omitted_phase_tas.(cwa_phase_names(p_i)+"_"+this_el).(mouse_name) = tmp_tas;
            end
        end
    end
    omitted_phase_cwas.(cue_id{ci}) = omitted_phase_cwa;
    omitted_phase_tass.(cue_id{ci}) = omitted_phase_tas;
end

% 
batch_cwa = omitted_phase_cwas;
batch_tas = omitted_phase_tass;
clear("omitted_phase_cwas","omitted_phase_tass");

batch_tas_new = struct;
for ci = 1:2
    fnames = string(fields(batch_tas.("cue"+ci))');
    for fname = fnames
        batch_tas_new.(fname) = batch_tas.("cue"+ci).(fname);
    end
end
batch_tas = batch_tas_new;
clear("batch_tas_new");
session_info = cwa_raw_day_info;
phases = ["late","LED","Tone"]; phases_text = ["late","LEDomi","Toneomi"]; % late, LED omi, Tone omi

cue_name = ["LED","Tone"];

% also group lick_info into batches (dont need lick here)
batch_lick_info = struct;
for ci = 1:2
    for m_i = 1:length(mouse_names)
        this_m = mouse_names(m_i);
        mouse_day_info = session_info(string(session_info(:,1)) == this_m,:);
        for p_i = 1:length(phases)
            phase = phases(p_i);
            if contains(phase,["early","late"])
                this_day = include_days.(this_m+phase);
            else
                tmp_id = find(contains(string(mouse_day_info([5,6])),phase));
                tmp_day = cell2mat(mouse_day_info(2:4));
                this_day = (tmp_day(tmp_id)+1):tmp_day(tmp_id+1);
            end
    
            batch_lick_across_day = [];
            batch_lick_across_day_n = nan(1,length(this_day));
            for d_i = 1:length(this_day)
                d = this_day(d_i);
                this_lick = lick_info.(this_m).single_trial_taking_all_ITI_frame.("day"+d).("cue"+ci+"_lick");
                n_batch = ceil(size(this_lick,1));
                batch_trial_id = [0:n_batch]; 
                batch_trial_id(end) = size(this_lick,1);
                % modify this_lick to get iti_dur etc
                this_lick{:,2} = this_lick{:,3}./this_lick{:,4};
                this_lick{:,2}(isnan(this_lick{:,2})) = 0;
                tmp_lick_batch = nan([n_batch,3]); % iti_rate, pre_rate, index
                for b_i = 1:n_batch
                    tmp_trial_id = (batch_trial_id(b_i)+1):batch_trial_id(b_i+1);
                    tmp = this_lick{tmp_trial_id,[2,3,6]};
                    tmp_lick_batch(b_i,2) = mean(tmp(:,3),"omitmissing");
                    tmp_lick_batch(b_i,1) = sum(tmp(:,2),[],"omitmissing")/sum(tmp(:,1),[],"omitmissing");
                end
                tmp_lick_batch(isnan(tmp_lick_batch(:,1)),1) = 0;
                tmp_lick_batch(:,3) = log(tmp_lick_batch(:,2)).*(tmp_lick_batch(:,2)-tmp_lick_batch(:,1))./(tmp_lick_batch(:,2)+tmp_lick_batch(:,1));
                batch_lick_across_day = cat(1,batch_lick_across_day,tmp_lick_batch);
                batch_lick_across_day_n(d_i) = size(tmp_lick_batch,1);
            end
            tmp = batch_lick_across_day(:,3);
            batch_lick_across_day(isnan(tmp)|isinf(tmp),3) = 0;
            batch_lick_info.("cue"+ci+phases_text(p_i)).(this_m).lick = batch_lick_across_day;
            batch_lick_info.("cue"+ci+phases_text(p_i)).(this_m).lick_n_trial = batch_lick_across_day_n;
        end
    end
end

for ci = 1:2
    cit = "cue"+ci; cit1 = "cue"+ci+cue_name(ci)+"omi";
    % get all ROIs whose rebound increases from late to omission
    across_mice_struct = struct;
    across_mice_struct.mouse_names = []; % hmmmmm i should just use mouse_names here. Anyway, just be more "precise"
    across_mice_struct.table = [];
    across_mice_struct.diff_mu = [];
    across_mice_struct.diff_ranksum_p = [];
    
    for m_i = 1:length(mouse_names)
        this_m = mouse_names(m_i);
        % get this_tas data used for ranksum
        late_vs_omi = struct;
        for p_name = [cit+"late",cit1]
            this_mu = squeeze(tas.(p_name).(this_m).mu_location_value(3,:,1));
            this_single = squeeze(tas.(p_name).(this_m).single_values(3,:,:));
            this_sig = logical(tas.(p_name).(this_m).mu_location_value_significance(3,:));
    
            n_trial = size(this_single,2);
            this_mu(~this_sig) = nan;
            this_single(~this_sig,:) = nan(sum(~this_sig),n_trial);
    
            late_vs_omi.(p_name).mu = this_mu';
            late_vs_omi.(p_name).single = this_single;
        end
    
        % get this_table and in_striatum filter
        this_table = ct_table(ct_table{:,"mouse_name"}==this_m,:);
        in_striatum_bit = logical(this_table{:,"significance"})';
        this_table = this_table(in_striatum_bit,:);
        diff_mu = common_functions.nan_minus(late_vs_omi.(cit1).mu(in_striatum_bit),late_vs_omi.(cit+"late").mu(in_striatum_bit));
        diff_ranksum_p = nan(sum(in_striatum_bit),1);
        tmp = find(in_striatum_bit);
        for r_i = 1:sum(in_striatum_bit)
            diff_ranksum_p(r_i) = common_functions.nan_ranksum(late_vs_omi.(cit1).single(tmp(r_i),:),late_vs_omi.(cit+"late").single(tmp(r_i),:));
        end
    
        % cat to across mouse struct
        across_mice_struct.mouse_names = cat(1,across_mice_struct.mouse_names,this_m);
        across_mice_struct.table = cat(1,across_mice_struct.table,this_table);
        across_mice_struct.diff_mu = cat(1,across_mice_struct.diff_mu,diff_mu);
        across_mice_struct.diff_ranksum_p = cat(1,across_mice_struct.diff_ranksum_p,diff_ranksum_p);
    end

    % filter diff_mu and ranksum_p to get target ROIs
    p_threshold = 0.01;
    diff_mu_threshold = 0;
    target_roi_bit = across_mice_struct.diff_mu > diff_mu_threshold & across_mice_struct.diff_ranksum_p < p_threshold;
    target_rois = across_mice_struct.table(target_roi_bit,["mouse_name","ROI_original"]);
    
    % smoothing trace before using
    plot_color = lines(4);
    lick_index_type = 3;
    smooth_detail = {30,'lowess',0};
    cusum_detail = {3,2};
    diff_chang_point_across = {};
    for m_i = 1:length(mouse_names)
        this_m = mouse_names(m_i);
        this_rois = target_rois{target_rois{:,"mouse_name"}==this_m,"ROI_original"};
        if isempty(this_rois)
            continue
        end

        re_single = permute(batch_tas.(cit1).(this_m).single_values(3,this_rois,:,1),[3,2,1,4]);
        lick_single = batch_lick_info.(cit1).(this_m).lick(:,lick_index_type);

        % smooth data
        lick_smoothed = smooth(lick_single,smooth_detail{:});
        re_smoothed = nan(size(re_single));
        for ri = 1:size(re_single,2)
            re_smoothed(:,ri) = smooth(re_single(:,ri),smooth_detail{:});
        end
        
        % change detection
        real_control = 15;
        cusum_lick_mean = mean(lick_smoothed(1:real_control),"omitmissing");
        cusum_lick_std = std(lick_smoothed(1:real_control),[],"omitmissing");
        cusum_re_mean = mean(re_smoothed(1:real_control,:),1,"omitmissing");
        cusum_re_std = std(re_smoothed(1:real_control,:),1,"omitmissing");
        % if no previous sessions are used, no need for rigid threshold
        real_control = 0;
        
        % find changing point of licking and ROI Fc via cusum
        [~,~,lick_uppersum,lick_lowersum] = cusum(lick_smoothed,cusum_detail{:},cusum_lick_mean,cusum_lick_std);
        cusum_lick_control = min([lick_lowersum(1:real_control);-cusum_detail{1}*cusum_lick_std],[],"all");
        lick_smooth_ilower = find(lick_lowersum<cusum_lick_control,1,"first");
        n_rois = length(this_rois);
        title_text = this_m + " " + string(this_rois);
        
        change_point = cell(1,length(this_rois)+1); % {lick_change,roi_changes}
        change_point{1} = lick_smooth_ilower;
        for ri = 1:length(this_rois)
            r = ri;
            [~,~,re_uppersum,re_lowersum] = cusum(re_smoothed(:,ri),3,1,cusum_re_mean(ri),cusum_re_std(ri));
            re_cusum_control = max([re_uppersum(1:real_control);cusum_detail{1}*cusum_re_std(ri)],[],"all");
            re_smooth_iupper = find(re_uppersum>re_cusum_control,1,"first");
            change_point{ri+1} = re_smooth_iupper;
        end

        diff_chang_point_across = cat(2,diff_chang_point_across,cellfun(@(x) x-change_point{1},change_point(2:end),UniformOutput=false));

        % plot smoothed trace only
        if ci==1 && this_m=="G23"
            plot_rois = 36;
            save_tag = {'fig4M_LED','fig4N_LED'};
        elseif ci==2 && this_m=="G12"
            plot_rois = 5;
            save_tag = {'fig4O_Tone','fig4P_Tone'};
        else
            continue
        end

        nr = 1; nc = 1;
        k = 0; fg = 0;
        for ri = 1:length(this_rois)
            r = ri;
            if k == 0
                fg = fg + 1;
                k = k + 1;
                fig = figure(Position=[100,100,1600,900]);
                tiled = tiledlayout(nr,nc,Parent=fig,TileSpacing="none");
                sgtitle(tiled,this_m+" "+cit)
            end
            axs = gobjects(1,1);
            for ax_i = 1:1
                axs(ax_i) = nexttile(tiled,k+ax_i-1);
            end
            hold(axs,"on")
            plot(axs(1),re_smoothed(:,ri),color=plot_color(2,:),LineWidth=2)
            scatter(axs(1),change_point{ri+1},re_smoothed(change_point{ri+1},ri),250,"x",MarkerEdgeColor=plot_color(3,:),LineWidth=2)
            ylabel(axs(1),"pk2 DF/F")
            yyaxis(axs(1),"right")
            plot(axs(1),lick_smoothed,color=plot_color(1,:),LineStyle='-',LineWidth=2)
            scatter(axs(1),change_point{1},lick_smoothed(change_point{1}),250,"x",MarkerEdgeColor=plot_color(4,:),LineWidth=2)
            ylabel(axs(1),"lick rate")
            yyaxis(axs(1),"left")
            hold(axs,"off")
            title(axs(1),title_text(ri))
            k = k+1;
            if k == nr * nc + 1 || ri == n_rois
                k = 0;
                saveas(fig,[cf,save_tag{1},'.png']);
                delete(fig)
            end
        end
    end

    % circle plot
    fig = figure(Position=[100,100,800,600]);
    tiled=tiledlayout(1,1);
    sgtitle(tiled,[cit,"trial # when pk2 start to increase minus trial # when lick id start to decrease"])
    ax = nexttile(tiled,1);

    % get trial # diff to plot
    tmp_diff_chang_point_across = diff_chang_point_across;
    tmp_diff_chang_point_across(cellfun(@isempty, tmp_diff_chang_point_across)) = {nan};
    circle_value = nan(size(across_mice_struct.diff_mu));
    circle_value(target_roi_bit) = cell2mat(tmp_diff_chang_point_across);
    circle_value = circle_value(~skipBubble);
    histo_frac = [sum(circle_value<0),sum(circle_value>0)];
    histo_frac = num2cell([histo_frac,histo_frac/sum(~isnan(circle_value))*100]);
    histo_frac = histo_frac([1,3,2,4]);
    % only plot ~skipBubble rois
    hold(ax,"on")
    histogram(ax,circle_value(~isnan(circle_value)),binwidth=10,Normalization="count");
    xline(ax,0,LineWidth=2)
    hold(ax,"off")
    new_xlim = max(abs(ax.XLim));
    xlim(ax,[-new_xlim,new_xlim])
    xlabel(ax,["diff trial #",sprintf("<0:%d(%2.2f%%) >0:%d(%2.2f%%)",histo_frac{:})])
    saveas(fig,[cf,save_tag{2},'.png']);
    delete(fig)
end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% functions
function out = get_early_late_days(index_measure,phase_name,mouse_name)
    % get early late day id based for phase+mouse based on index_measure (and set test alpha=0.05)
    LEDomi_mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
    
    cue1lickindex_LEDomi = array2table(...
        {[1],[2,3,4,5];[1,2,3,4],[5,6];[1,2,3],[4,5,6];[1,2],[3,4];...
        [1,2],[3,4,5];[1],[2,3,4,5];[1,2,3],[4,5];[1,2],[3,4,5]},...
        RowNames=LEDomi_mouse_names,VariableNames=["early","late"]);

    cue1prelick_LEDomi = array2table(...
        {[1],[2,3,4,5];[1,2,3],[4,5,6];[1],[2,3,4,5,6];[1,2],[3,4];...
        [1,2,3],[4,5];[1],[2,3,4,5];[1],[2,3,4,5];[1,2],[3,4,5]},...
        RowNames=LEDomi_mouse_names,VariableNames=["early","late"]);

    early_late_day_id = struct;
    early_late_day_id.cue1_lick_index.LEDomi = cue1lickindex_LEDomi;
    early_late_day_id.cue1_pre_lick.LEDomi = cue1prelick_LEDomi;
    

    switch index_measure
        case "cue1_pre_lick"
            out = early_late_day_id.cue1_pre_lick.(phase_name){mouse_name,:};
        case "cue1_index"
            out = early_late_day_id.cue1_lick_index.(phase_name){mouse_name,:};
        case "manual"
            error("TBI.")
    end
end

function [pk_pks,pk_locs,pk_pks_single,pk_significance] = get_peak(roi_traces,pk_window,min_amp_threshold,varargin)
    exclude_peaks = [];

    ip = inputParser;
    ip.addParameter('exclude_peaks',[]); % if this is not empty, it means to find rebounds
    ip.parse(varargin{:})
    for j = fields(ip.Results)'
        eval([j{1},'=ip.Results.',j{1},';']);
    end

    n_trial = size(roi_traces,3);
    n_rois = size(roi_traces,2)-3;
    roi_traces_mu = mean(roi_traces,3,"omitmissing");
    if isempty(exclude_peaks) % if is to find peak, find it directly
        pk_cell = num2cell(roi_traces_mu(pk_window,:),1);
        [pk_pks,pk_locs] = cellfun(@(x) findpeaks(x,pk_window,NPeaks=1,SortStr="descend"),pk_cell,UniformOutput=0);
    else % if is to find rebound, start from the first frame when trace goes below control since peak
        pk_cell = cell([1,n_rois+3]);
        pk_window_cell = cell([1,n_rois+3]);
        exclude_peaks = [nan,nan,nan,exclude_peaks];
        for ri = 1:n_rois+3
            if isnan(exclude_peaks(ri))
                pk_cell{ri} = roi_traces_mu(pk_window,ri);
                pk_window_cell{ri} = pk_window;
            else
                % find the first frame where
                tmp_pk_window = (exclude_peaks(ri)+1):pk_window(end);
                reset_frame = find(roi_traces_mu(tmp_pk_window,ri) <= min_amp_threshold,1,"first");
                tmp_pk_window = (reset_frame+tmp_pk_window(1)-1):tmp_pk_window(end);
                if length(tmp_pk_window) >= 3 % if enough frames exist, use them; otherwise fake it
                    pk_cell{ri} = roi_traces_mu(tmp_pk_window,ri);
                    pk_window_cell{ri} = tmp_pk_window;
                else
                    pk_cell{ri} = [nan,nan,nan];
                    pk_window_cell{ri} = [1,2,3];
                end
            end
        end
        [pk_pks,pk_locs] = cellfun(@(x,y) findpeaks(x,y,NPeaks=1,SortStr="descend"),pk_cell,pk_window_cell,UniformOutput=0);
    end
    pk_locs_cell = pk_locs;
    empty_id = cellfun(@(x) isempty(x),pk_pks);
    if any(empty_id)
        [pk_pks{empty_id}] = deal(nan);
        [pk_locs{empty_id}] = deal(nan);
    end
    pk_pks = cell2mat(pk_pks);
    pk_locs = cell2mat(pk_locs);
    % apply amplitude threshold
    pk_significance = zeros(size(pk_pks));
    pk_significance(pk_pks>=min_amp_threshold) = 1;    

    pk_pks = pk_pks(4:end);
    pk_locs = pk_locs(4:end);
    pk_significance = pk_significance(4:end);
    pk_locs_cell = pk_locs_cell(4:end);
    
    tmp = num2cell(roi_traces(:,4:end,:),[1,3]);
    tmp1 = cellfun(@(x,y) squeeze(x(y,:,:)),tmp,pk_locs_cell,UniformOutput=0);
    empty_id = cellfun(@(x) isempty(x),tmp1);
    pk_pks_single = tmp1;
    [pk_pks_single{empty_id}] = deal(nan(n_trial,1));
    pk_pks_single = cell2mat(pk_pks_single);
end
