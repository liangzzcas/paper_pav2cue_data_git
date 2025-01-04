%% loading raw data and training GLM to regress out movement related noise
close all;clear;clc;
cf = [pwd,'\'];
raw_path = [cf,'raw_data\'];
CT = readtable([raw_path,'CT_across_GXX_mice.xlsx']);
session_info = common_functions.get_training_info();
mouse_names_ori = string(unique(CT{:,"mouse_name"})');

% hyper parameters
Fourier_hz = [nan,0.3];
output_data = struct;
for mouse_name = mouse_names_ori
    mouse_path = [raw_path,char(mouse_name),'\'];
    tmp = dir(mouse_path);
    tmp = string({tmp.name});
    dates = tmp(~contains(tmp,[".",".."])); clear("tmp");
    mouse_struct = struct;
    data_cated_variablenames = ["sample_rate","wavelength","ITI_bit","vel","acc","470","405","vel_original"];
    data_cated = cell(length(dates),length(data_cated_variablenames));
    load_ids = 1:length(dates);
    for i = 1:length(dates)
        di = dates(i);
        [ROI405,ROI470,behav5,behav7,~,~,tag_405470570,~,~] = ROI_loader([mouse_path,char(di),'\'],Fourier = Fourier_hz,Fourier_1=Fourier_hz);
        sr = common_functions.get_sr(behav7);
        loco = common_functions.get_overall_vel_acc(behav7,...                
                'sampling_freq',sr,...
                'lowpass_freq',1.5,...
                'smooth_window',0.5,...
                'filter_flag',1);

        tmp_ta_input_470 = ROI470.Fc;
        tmp_ta_input_405 = [];
        ITI_bit = behav7.stimulus_led | behav7.stimulus_sound2;
        data_cated{i,1} = sr;
        data_cated{i,2} = tag_405470570;
        data_cated{i,3} = ITI_bit;
        data_cated{i,4} = loco.vel;
        data_cated{i,5} = loco.acc;
        data_cated{i,6} = tmp_ta_input_470;
        data_cated{i,7} = tmp_ta_input_405;
        data_cated{i,8} = loco.vel_original;
    end
    mouse_struct.data = cell2table(data_cated,VariableNames=data_cated_variablenames);
    output_data.(mouse_name) = mouse_struct;
end

data_ori = output_data;
mouse_names = string(fields(data_ori)');
additional_info = struct;
for mouse_name = mouse_names
    % sr info
    if contains(mouse_name,["G12","G15"])
        additional_info.(mouse_name).sr = 30;
    else
        additional_info.(mouse_name).sr = 18;
    end
    additional_info.(mouse_name).include_405 = 0;
end

flag_recons_plot_during_training = 1;

% set hyperparam
xcorr_maxlag_s = [0.5,2]; %[window used, window plotted]
glm_mdls_mice = struct;
for mouse_name = mouse_names
    this_data = data_ori.(mouse_name).data;
    this_additional_info = additional_info.(mouse_name);
    n_rois = size(this_data{1,"470"}{1},2);

    glm_data_id = logical(this_data{:,1} == this_additional_info.sr);
    glm_data_cated = cat_data_across_day(this_data,glm_data_id);
    glm_data_sr = glm_data_cated{1,1};
    % get glm behav and fc data
    xcorr_maxlag = round(xcorr_maxlag_s * glm_data_sr);
    glm_vel = glm_data_cated{1,4}{1}(glm_data_cated{1,3}{1});
    glm_acc = glm_data_cated{1,5}{1}(glm_data_cated{1,3}{1});
    glm_fc470s = glm_data_cated{1,6}{1}(glm_data_cated{1,3}{1},:);
    glm_fc405s = [];
    % loop thru roi to get glm
    glm_mdls_mouse = struct;
    glm_mdls_mouse.data = glm_data_cated;
    glm_mdls_xcorr_tags = ["lags_frame_plotted","lags_frame_considered","vel_auto_corr_coef","acc_auto_corr_coef",...
        "vel_lag_pos_frame","vel_lag_neg_frame","acc_lag_pos_frame","acc_lag_neg_frame"];
    glm_mdls_xcorr = cell(1,length(glm_mdls_xcorr_tags));
    glm_mdls_xcorr_arr = nan((xcorr_maxlag(2)*2+1)*2+4,n_rois);
    glm_mdls_mdls = cell(1,n_rois);

    glm_input_table_tags = ["fc470","fc405","vel_zero","acc_zero","vel_pos","vel_neg","acc_pos","acc_neg"];
    for ri = 1:n_rois
        r = ri;
        glm_input_table_bit = false(size(glm_input_table_tags));
        glm_input_table_bit([3,4]) = true;
        % xcorr to get delay between vel,acc vs signal
        [vel_I_pos,vel_I_neg,vel_r] = find_lags_using_xcorr(glm_fc470s(:,r),glm_vel,xcorr_maxlag);
        [acc_I_pos,acc_I_neg,acc_r] = find_lags_using_xcorr(glm_fc470s(:,r),glm_acc,xcorr_maxlag);
        glm_mdls_xcorr_arr(:,ri) = [vel_r;acc_r;vel_I_pos;vel_I_neg;acc_I_pos;acc_I_neg];
        % set and train glm
        behav_lags = [vel_I_pos,vel_I_neg,acc_I_pos,acc_I_neg];
        glm_input_table_bit(5:8) = ~isnan(behav_lags);
        truncate = {min(behav_lags(behav_lags<0),[],"omitmissing"),max(behav_lags(behav_lags>0),[],"omitmissing")};
        fc_mask = [1,length(glm_fc470s(:,r))];
        if ~isempty(truncate{1})
            fc_mask(1) = 1-truncate{1};
        end
        if ~isempty(truncate{2})
            fc_mask(2) = length(glm_fc470s(:,r)) - truncate{2};
        end
        n_frame_truncated = fc_mask(2) - fc_mask(1) + 1;
        glm_input_table = nan(n_frame_truncated,length(glm_input_table_tags));
        glm_input_table(:,1) = glm_fc470s(fc_mask(1):fc_mask(2),r);
        glm_input_table(:,2) = nan(n_frame_truncated,1);
        glm_input_table(:,3) = glm_vel(fc_mask(1):fc_mask(2));
        glm_input_table(:,4) = glm_acc(fc_mask(1):fc_mask(2));

        for behav_lags_i = 1:length(behav_lags)
            if ~isnan(behav_lags(behav_lags_i))
                if behav_lags_i<=2
                    this_behav = glm_vel;
                else
                    this_behav = glm_acc;
                end
                glm_input_table(:,behav_lags_i+4) = this_behav(fc_mask(1)+behav_lags(behav_lags_i):fc_mask(1)+behav_lags(behav_lags_i)+n_frame_truncated-1);
            end
        end
        glm_input_table = array2table(glm_input_table,VariableNames=glm_input_table_tags);
        mdlSpec = ['fc470 ~ ',char(strjoin(glm_input_table_tags(glm_input_table_bit)," + ")),' + 1'];
        mdl = fitglm(glm_input_table,mdlSpec);
        glm_mdls_mdls{ri} = extract_model(mdl);
    end
    glm_mdls_xcorr(1,[1,2]) = {xcorr_maxlag(2),xcorr_maxlag(1)};
    glm_mdls_xcorr(1,[3:end]) = mat2cell(glm_mdls_xcorr_arr,[xcorr_maxlag(2)*2+1,xcorr_maxlag(2)*2+1,1,1,1,1]);
    glm_mdls_xcorr = cell2table(glm_mdls_xcorr,VariableNames=glm_mdls_xcorr_tags);
    glm_mdls_mouse.xcorr = glm_mdls_xcorr;
    glm_mdls_mouse.mdls = glm_mdls_mdls;

    glm_mdls_mice.(mouse_name) = glm_mdls_mouse;
end

% regress out movement related noise from data
data_path = [cf,'processed_and_organized_data\'];
mouse_names = ["G12","G15","G17","G19","G22","G21","G23","G24","G25","G26","G27"];

ta_highpassed = struct;
ta_highpassed.cwa_raw_day_info = common_functions.get_training_info();
ta_highpassed_ITI_licking = struct;
for mouse_name = mouse_names
    mouse_path = [raw_path,char(mouse_name),'\'];
    tmp = dir(mouse_path);
    tmp = string({tmp.name});
    dates = tmp(~contains(tmp,[".",".."])); clear("tmp");
    load_ids = 1:length(dates);
    ta_highpassed.cwa_raw.(mouse_name) = get_triavg(mouse_path, mouse_name, dates, load_ids, load_ft_fr=[0.3,nan], GLM=glm_mdls_mice.(mouse_name));
    % get tri ave of ITI licking onset events
    tmp = ta_highpassed.cwa_raw_day_info;
    this_session_info = tmp([tmp{:,1}]==mouse_name,2:end);
    ta_highpassed_ITI_licking.(mouse_name) = get_triavg_ITI_lick(mouse_path, mouse_name, dates, load_ids, this_session_info, load_ft_fr=[0.3,nan], GLM=glm_mdls_mice.(mouse_name));
end
save([data_path,'event_aligned_highpass03.mat'],"ta_highpassed",'-v7.3');
save([data_path,'ITI_licking_aligned_highpass03.mat'],"ta_highpassed_ITI_licking",'-v7.3');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% getting licking index for each day from behav file
close all;clear;clc;
cf = [pwd,'\'];
raw_path = [cf,'raw_data\'];
% parameters
index_type = 2;
lick_window = [0,3]; % x second before cue and rew
second_exclude_since_reward_onset = [0,3];
mouse_names = ["G12","G15","G17","G19","G22","G21","G23","G24","G25","G26","G27","glu924","glu926"];
across_mouse_struct = struct;
for m_i = 1:length(mouse_names)
    mouse_name = mouse_names(m_i);
    mouse_path = [raw_path,char(mouse_name),'\'];
    tmp = dir(mouse_path);
    tmp = string({tmp.name});
    dates = tmp(arrayfun(@(x) any(regexp(x,'^\d{6}$')),tmp)); clear("tmp");
    ids = 1:length(dates);
    % preallocate behav cell
    behav_cell = cell(1,length(ids));
    % assign behav path
    for d_i = 1:length(ids)
        [ROI405i,ROI470i,behav5i,behav7i,date_load,~,image_channel] = ROI_loader([mouse_path,char(dates(di)),'\']);
        behav_cell{d_i} = behav7i;
    end
    output = get_lick_index(behav_cell, mouse_name=mouse_name, lick_window=lick_window, second_exclude_since_reward_onset=second_exclude_since_reward_onset, lick_index_type = index_type,lick_onset_gap=[]);
    across_mouse_struct.lick.(mouse_name) = output;
end
save([cf,'processed_and_organized_data\across_mice_lick_index_data_whole_ITI.mat'],'-struct',"across_mouse_struct");

% repeat for DA mice
mouse_names = ["DL18","DL20","DL21","DL23"];
across_mouse_struct = struct;
for m_i = 1:length(mouse_names)
    mouse_name = mouse_names(m_i);
    mouse_path = [raw_path,char(mouse_name),'\'];
    tmp = dir(mouse_path);
    tmp = string({tmp.name});
    dates = tmp(arrayfun(@(x) any(regexp(x,'^\d{6}$')),tmp)); clear("tmp");
    ids = 1:length(dates);
    % preallocate behav cell
    behav_cell = cell(1,length(ids));
    % assign behav path
    for d_i = 1:length(ids)
        [ROI405i,ROI470i,behav5i,behav7i,date_load,~,image_channel] = ROI_loader([mouse_path,char(dates(di)),'\']);
        behav_cell{d_i} = behav7i;
    end
    output = get_lick_index(behav_cell, mouse_name=mouse_name, lick_window=lick_window, second_exclude_since_reward_onset=second_exclude_since_reward_onset, lick_index_type = index_type,lick_onset_gap=[]);
    across_mouse_struct.lick.(mouse_name) = output;
end
save([cf,'processed_and_organized_data\DA_across_mice_lick_index_data_whole_ITI.mat'],'-struct',"across_mouse_struct");

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% reorganize data into different learning phases, determined via the licking of each mouse
close all;clear;clc;
cf = [pwd,'\'];
ta_data_all  = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']).cwa_raw;
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
training_info = common_functions.get_training_info();
include_days = common_functions.get_include_days();
main_sr = common_functions.get_main_samplerate();
plot_components_combinations = ["unpred","unpred";"unpred","early";"unpred","late";"unpred","LEDomi";"unpred","Toneomi";...
    "cue1","early";"cue1","late";"cue2","early";"cue2","late";"rew1","early";"rew1","late";"rew2","early";"rew2","late";...
    "cue1","LEDomi";"cue1","Toneomi";"cue2","LEDomi";"cue2","Toneomi";"rew1","LEDomi";"rew1","Toneomi";"rew2","LEDomi";"rew2","Toneomi"];
mouse_names = ["G12","G15","G17","G19","G22","G21","G23","G24","G25","G26","G27"];

for plot_components_combination_i = 1:size(plot_components_combinations,1) % usually size(plot_components_combinations,1)==1, only exception is deprecated here
    plot_components_combination = plot_components_combinations(plot_components_combination_i,:);
    plot_components = plot_components_combination(1); plot_components_early_late = plot_components_combination(2);

    if contains(plot_components,"unpred")
        cue_rew_unpred = 'unpred';
    elseif contains(plot_components,"rew")
        cue_rew_unpred = 'rew';
    elseif contains(plot_components,"cue")
        cue_rew_unpred = 'cue';
    end

    for mouse_name_i = 1:length(mouse_names)
        mouse_name = mouse_names(mouse_name_i);
        % skip some cases
        if contains(mouse_name,["G25","G26","G27"]) && plot_components_early_late == "Toneomi"
            continue
        end

        % get ta
        switch plot_components_early_late
            case "unpred"
                ids = 1:training_info{[training_info{:,1}]==mouse_name,4};
                load_days = ids;
            case "late"
                ids = include_days.(mouse_name+"late");
                load_days = ids;
            case "early"
                ids = include_days.(mouse_name+"early");
                load_days = ids;
            case "omission"
                ids = training_info{[training_info{:,1}]==mouse_name,[2,4]};
                ids = [ids(1)+1:ids(2)];
                load_days = ids;
            case "LEDomi"
                ids = include_days.(mouse_name+"LEDomi");
                load_days = ids;
            case "Toneomi"
                ids = include_days.(mouse_name+"Toneomi");
                load_days = ids;
        end

        for day_group_i = 1:size(load_days,1)
            load_day = load_days(day_group_i,:);
            load_ids = 1:length(load_day);

            ta = struct;
            for i = 1:length(load_day)
                ta.("file"+i) = ta_data_all.(mouse_name).("file"+load_day(i));
            end

            % get and cat unpred rew/pred rew/cue
            ta_unfilter = get_ta_total(ta,load_ids,plot_components,sample_rate=main_sr(mouse_name));
            if strcmp(plot_components,plot_components_early_late)
                ta_across.(plot_components).(mouse_name) = ta_unfilter;
            else
                ta_across.(plot_components+plot_components_early_late).(mouse_name) = ta_unfilter;
            end
        end
    end
end

rebaseline_flag = 1;
plot_components = ["unpred","cue","rew"];

rebaseline_info = struct;
ta_fields = string(fields(ta_across)');

components_window_activity = struct;
for ta_i = 1:length(ta_fields)
    field_name = ta_fields(ta_i);
    tmp0 = ["cue","rew","unpred"];
    tmp1 = arrayfun(@(x) contains(field_name,x),tmp0);
    plot_component_str = tmp0(tmp1); clear("tmp0","tmp1");
    ta_mouse = ta_across.(field_name);
    ta_mouse_names = string(fields(ta_mouse)');

    components_window_activity.(field_name) = ta_mouse;
    mice_name = ta_mouse_names;
    
    for mouse_name_i = mice_name
        fr = main_sr(string(mouse_name_i));
        % get mu and std
        mu_m = ta_mouse.(mouse_name_i).mu;
        std_m = std(mu_m(1:fr,4:end),[],1,"omitnan");
        std_m = cat(2,[nan,nan,nan],std_m);
        components_window_activity.(field_name).(mouse_name_i).null = std_m;
        % get mu and std for rew consumption (if consumption exist)
        if isfield(ta_mouse.(mouse_name_i),"mu_consump") && ~isempty(ta_mouse.(mouse_name_i).mu_consump)
            rewconsump_tag = 1;
            mu_m = ta_mouse.(mouse_name_i).mu_consump;
            std_m = std(mu_m(1:fr,4:end),[],1,"omitnan");
            std_m = cat(2,[nan,nan,nan],std_m);
            components_window_activity.(field_name).(mouse_name_i).null_consump = std_m;
        end
        % get mu and std for omission (if omission exist)
        if isfield(ta_mouse.(mouse_name_i),"mu_omi") && ~isempty(ta_mouse.(mouse_name_i).mu_omi)
            mu_m = ta_mouse.(mouse_name_i).mu_omi;
            std_m = std(mu_m(1:fr,4:end),[],1,"omitnan");
            std_m = cat(2,[nan,nan,nan],std_m);
            components_window_activity.(field_name).(mouse_name_i).null_omi = std_m;
        end

        % get rebaseline (based on before cue onset only) if rebaseline flag is on
        if rebaseline_flag == 1
            activity = ta_mouse.(mouse_name_i).activity;
            if plot_component_str == "cue"
                activity_shift = mean(activity(1:fr,4:end,:),[1,3],"omitnan");
                rebaseline_info.(field_name).(mouse_name_i) = activity_shift;
            elseif plot_component_str == "rew"
                rebaseline_load_name = char(field_name);
                rebaseline_load_name(1:3) = 'cue';
                rebaseline_load_name = string(rebaseline_load_name);
                assert(isfield(rebaseline_info,rebaseline_load_name),"Process cue before reward to get rebaseline info for reward.")
                activity_shift = rebaseline_info.(rebaseline_load_name).(mouse_name_i);
            elseif plot_component_str == "unpred"
                activity_shift = mean(activity(1:fr,4:end,:),[1,3],"omitnan");
            end
            activity_shift = repmat(activity_shift,[size(activity,1),1,size(activity,3)]);

            activity(:,4:end,:) = activity(:,4:end,:) - activity_shift;
            components_window_activity.(field_name).(mouse_name_i).activity = activity;
            components_window_activity.(field_name).(mouse_name_i).mu = mean(activity,3);

            % shift baseline for omission
            if isfield(ta_mouse.(mouse_name_i),"mu_omi") && ~isempty(ta_mouse.(mouse_name_i).mu_omi)
                activity_omi = ta_mouse.(mouse_name_i).activity_omi;
                activity_shift = repmat(activity_shift(1,:,1),[size(activity_omi,1),1,size(activity_omi,3)]);
                activity_omi(:,4:end,:) = activity_omi(:,4:end,:) - activity_shift;
                components_window_activity.(field_name).(mouse_name_i).activity_omi = activity_omi;
                components_window_activity.(field_name).(mouse_name_i).mu_omi = mean(activity_omi,3);
            end

            % shift baseline for consumption
            if isfield(ta_mouse.(mouse_name_i),"mu_consump") && ~isempty(ta_mouse.(mouse_name_i).mu_consump)
                if plot_component_str == "unpred"
                    null_end = round(0.6*fr); % for consumption, signal might come before onset
                    this_activity_shift = mean(ta_mouse.(mouse_name_i).activity_consump(1:null_end,4:end,:),[1,3],"omitnan");
                else
                    this_activity_shift = activity_shift(1,:,1);
                end
                activity_consump = ta_mouse.(mouse_name_i).activity_consump;
                this_activity_shift = repmat(this_activity_shift,[size(activity_consump,1),1,size(activity_consump,3)]);
                activity_consump(:,4:end,:) = activity_consump(:,4:end,:) - this_activity_shift;
                components_window_activity.(field_name).(mouse_name_i).activity_consump = activity_consump;
                components_window_activity.(field_name).(mouse_name_i).mu_consump = mean(activity_consump,3);
            end
        end
    end
end

[cwa_delivery,cwa_consump] = separate_delivery_consumption(components_window_activity);
save([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat'],"-struct","cwa_consump")

% identify transient components from triggered average
mu_std_factor = 3;
components_window = cwa_consump;
field_windows = common_functions.get_comp_timewindow();

[~,tri_avg_single] = generate_tas_from_cwa(components_window,ct_table,field_windows=field_windows,mu_std_factor=mu_std_factor);
save([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat'],"-struct","tri_avg_single")


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% merge rew1 and rew2
close all;clear;clc;
cf = [pwd,'\'];
seperated_cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
mouse_names = string(fields(seperated_cwa.rew1late)');
p_names = ["rew1early","rew2early","rew1late","rew2late"];

main_srs = common_functions.get_main_samplerate();
field_windows = common_functions.get_comp_timewindow();
mu_std_factor = 3;

e_l_text = ["early","late"];
% build cwa
merged_cwa = seperated_cwa;
for mouse_name = mouse_names
    this_sr = main_srs(mouse_name);
    for e_l_i = 1:2
        new_rew_struct = struct;
        new_rew_struct.activity = cat(3,...
            seperated_cwa.("rew1"+e_l_text(e_l_i)).(mouse_name).activity,...
            seperated_cwa.("rew2"+e_l_text(e_l_i)).(mouse_name).activity);
        new_rew_struct.mu = mean(new_rew_struct.activity,3,"omitmissing");
        new_rew_struct.sem = std(new_rew_struct.activity,[],3,"omitmissing")/sqrt(size(new_rew_struct.activity,3));
        new_rew_struct.activity_omi = []; new_rew_struct.mu_omi = []; new_rew_struct.sem_omi = [];
        % combine std, just use math (or we can recalculate this from cue activity, will be the same)
        null_1 = seperated_cwa.("cue1"+e_l_text(e_l_i)).(mouse_name).null;
        null_2 = seperated_cwa.("cue2"+e_l_text(e_l_i)).(mouse_name).null;
        null_merged = sqrt(       (this_sr-1)*(null_1.^2+null_2.^2)/(2*this_sr-1)   +  this_sr^2*(null_1-null_1).^2/((2*this_sr)*(2*this_sr-1))     );
        new_rew_struct.null = null_merged;

        % assign back and replace rew1|rew2
        merged_cwa.("cue1"+e_l_text(e_l_i)).(mouse_name).null = null_merged;
        merged_cwa.("cue2"+e_l_text(e_l_i)).(mouse_name).null = null_merged;
        
        merged_cwa.("rew1"+e_l_text(e_l_i)).(mouse_name) = new_rew_struct;
        merged_cwa.("rew2"+e_l_text(e_l_i)).(mouse_name) = new_rew_struct;
    end
end

% build tas
[~,tri_avg_single] = generate_tas_from_cwa(merged_cwa,ct_table,field_windows=field_windows,mu_std_factor=mu_std_factor);

save([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump_rew_merged.mat'],'-struct','merged_cwa');
save([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat'],'-struct','tri_avg_single');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% reorganize ITI licking data into different learning phases, determined via the licking of each mouse
close all;clear;clc;
cf = [pwd,'\'];
lick_data = load([cf,'ITI_licking_aligned_highpass03.mat']);
CT = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
include_days = common_functions.get_include_day_ids(); % use version=3 here to get all pav days merged (named as late)
null_window_s = 0.7;
ITI_licking_comp_window = common_functions.get_ITI_licking_comp_timewindow();
null_std_factor = 3;

plot_data = build_plotdata_by_session(lick_data,include_days,null_window_s);
cwa_struct = cwa_table_to_struct(plot_data);
mouse_names = string(fields(plot_data)');

axs_title_texts = ["early","late","LEDomi","Toneomi"];
lick_comp_window_pk_dp_re_s = ITI_licking_comp_window+1;

tri_avg_single = specialized_add_to_tri_avg_place_holder();

for mi = 1:length(mouse_names)
    mouse_name = mouse_names(mi);
    this_tbl = plot_data.(mouse_name);
    session_names = string(this_tbl.Properties.RowNames');
    for i = 1:4
        s_i = find(contains(session_names,axs_title_texts(i)));
        if isempty(s_i)
            continue
        end
        roi_sr = this_tbl{s_i,1};
        roi_traces = cell2mat(this_tbl{s_i,2});
        roi_traces_mu = mean(roi_traces,3,"omitmissing");
        
        pk_dp_re_window = round(lick_comp_window_pk_dp_re_s * roi_sr);
        pk_window = pk_dp_re_window(1,1):pk_dp_re_window(1,2);
        dp_window = pk_dp_re_window(2,1):pk_dp_re_window(2,2);
        re_window = pk_dp_re_window(3,1):pk_dp_re_window(3,2);
        
        null_thres = cell2mat(this_tbl{s_i,3}) * null_std_factor;

        [pk_pks,pk_locs_frame,pk_pks_single,pk_sig] = get_peak(roi_traces,roi_traces_mu,pk_window,null_thres);
        pk_locs = pk_locs_frame/roi_sr-1;
        [dp_pks,dp_locs,dp_pks_single,dp_sig] = get_peak(-roi_traces,-roi_traces_mu,dp_window,null_thres);
        dp_locs = dp_locs/roi_sr-1;
        dp_pks = -dp_pks; dp_pks_single = -dp_pks_single;
        % get sig peak to be used in finding rebound
        tmp_pk_locs = pk_locs_frame;
        tmp_pk_locs(~logical(pk_sig)) = nan;
        [re_pks,re_locs,re_pks_single,re_sig] = get_peak(roi_traces,roi_traces_mu,re_window,null_thres,exclude_peaks=tmp_pk_locs);
        re_locs = re_locs/roi_sr-1;
        
        single_values = cat(3,pk_pks_single,dp_pks_single,re_pks_single);
        single_values = permute (single_values,[3,2,1]);
        mu_location_value = cat(3,cat(1,pk_pks,dp_pks,re_pks),cat(1,pk_locs,dp_locs,re_locs));
        mu_location_value_significance = cat(1,pk_sig,dp_sig,re_sig); % to be assigned

        tri_avg_single.(axs_title_texts(i)).(mouse_name).single_values = single_values;
        tri_avg_single.(axs_title_texts(i)).(mouse_name).mu_location_value = mu_location_value;
        tri_avg_single.(axs_title_texts(i)).(mouse_name).mu_location_value_significance = mu_location_value_significance;
    end
end
save([cf,'processed_and_organized_data\ITI_licking_components_window_activity_filtered_table.mat'],'-struct','plot_data');
save([cf,'processed_and_organized_data\ITI_licking_components_window_activity_filtered.mat'],'-struct','cwa_struct');
save([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat'],'-struct','tri_avg_single');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% Dopamine data: load and then reorganize data into different learning phases, determined via the licking of each mouse
close all;clear;clc;
cf = [pwd,'\'];
CT = readtable([cf,'raw_data\CT_across_DLXX_mice.xlsx']);

mouse_names_ori = ["DL18","DL20","DL21","DL23"];
fc_loading_fft_1 = [0.3,nan];
vel_type = "lin_vel";

% build data
output_filtered = struct;
for mouse_name = mouse_names_ori
    main_sr = 30;
    mouse_struct = struct;

    this_info = common_functions.get_training_info_DA(mouse_name);
    this_days = 1:max([cell2mat(this_info("LED_omi")),cell2mat(this_info("Tone_omi"))]);
    mouse_struct.include_days = this_days;
    
    tas_filtered = struct;
    for i = 1:length(this_days)
        di = this_days(i);
        [ROIBehav,~] = ROI_loader_DA([cf,'raw_data/',char(mouse_name),'/'],di);
        ROI470 = ROIBehav.roi; behav7 = ROIBehav.behavior; this_sr = common_functions.get_sr(behav7);
        ROI470_filtered = frequency_pass(ROI470,fc_loading_fft_1,this_sr);

        loco = common_functions.get_overall_vel_acc(behav7,...
            lowpass_freq=1.5, smooth_window=0.5, filter_flag=3, vel_type=vel_type);
        behavfc = cat(2,behav7.lick_count,loco.vel,loco.acc,ROI470.Fc);
        behavfc_filtered = cat(2,behav7.lick_count,loco.vel,loco.acc,ROI470_filtered.Fc);

        event_time = extract_task_info(ROIBehav.task_info);
        event_text = string(fields(event_time))';
        this_ta = struct;
        this_ta_filtered = struct;
        for ei = 1:length(event_text)
            this_e_text = event_text(ei);
            this_e = event_time.(this_e_text);
            this_ta.(this_e_text) = eventTriggeredAverage(behavfc,this_e(~isnan(this_e)),-(this_sr-1),3*this_sr);
            this_ta_filtered.(this_e_text) = eventTriggeredAverage(behavfc_filtered,this_e(~isnan(this_e)),-(this_sr-1),3*this_sr);
        end
        tas_filtered.("day"+di) = this_ta_filtered;
        % merge across days
        if i == 1
            tas_filtered.across = structfun(@(x) x.activity,this_ta_filtered,UniformOutput=false);
        else
            for ei = 1:length(event_text)
                this_e_text = event_text(ei);
                tas_filtered.across.(this_e_text) = cat(3,tas_filtered.across.(this_e_text),this_ta_filtered.(this_e_text).activity);
            end
        end
        keyboard
    end
    output_filtered.(mouse_name) = tas_filtered;
end
save([cf,'processed_and_organized_data\DA_event_aligned_highpass03.mat'],'-struct',"output_filtered")

% build cwa by dividing tri_ave into phases based on learning
cwa = struct;
for mouse_name = mouse_names_ori
    main_sr = 30;
    phase_info = common_functions.get_include_days_DA(mouse_name);
    p_names = phase_info.keys;

    mouse_struct = struct;
    for pi = 1:length(p_names)
        p_name = p_names(pi);
        ids = phase_info{p_name};
        merged_datas = struct;
        merged_datas_unfilter = struct;
        for di = 1:length(ids)
            id = ids(di);
            this_data = output_filtered.(mouse_name).("day"+id);
            merged_datas = merge_pav2cue_data(this_data,merged_datas,main_sr);
        end
        cwa.(mouse_name).(p_name) = merged_datas;
    end
end
save([cf,'processed_and_organized_data\DA_components_window_activity_filtered.mat'],'-struct',"cwa");

% identify transient components from triggered average
include_mask = CT{:,"fiber_bottom_AP"}>0 & CT{:,"fiber_bottom_ML"}<2 & CT{:,"fiber_bottom_DV"}<4;
CT(:,"fiber_include_bitmask") = num2cell(include_mask);
include_info = [sum(include_mask)/size(CT,1),sum(include_mask),size(CT,1)];

pnames_to_plot = ["pre","post","LED_omi","Tone_omi"];
main_sr = 30;

% put all ROIs togather
pnames = pnames_to_plot;
mouse_names = string(fields(cwa))';

data_struct = struct;
data_struct.plot_ct = {};
for mouse_name = mouse_names
    this_CT = CT(CT{:,"mouse_name"}==mouse_name,:);
    this_CT_bitmask = this_CT{:,"fiber_include_bitmask"};
    this_CT_bitmask_3 = [false;false;false;this_CT_bitmask];
    plot_CT = this_CT(this_CT_bitmask,:);
    data_struct.plot_ct = cat(1,data_struct.plot_ct,plot_CT);

    for pi = 1:length(pnames)
        p_name = pnames(pi);
        if ~isfield(data_struct,p_name)
            data_struct.(p_name).cue1 = {};
            data_struct.(p_name).cue2 = {};
        end
        for ci = 1:2
            c = "cue"+ci;
            this_act = cwa.(mouse_name).(p_name).cueOn.(c).activity(:,this_CT_bitmask_3,:);
            % rebaseline Fc
            tmp = mean(this_act(1:main_sr,:,:),[1,3],"omitmissing");
            tmp = repmat(tmp,[size(this_act,1),1,size(this_act,3)]);
            this_act = this_act-tmp;
            this_acts.(p_name).(c) = this_act;
            data_struct.(p_name).("cue"+ci) = cat(2,data_struct.(p_name).("cue"+ci),...
                num2cell(this_act,3));
        end
    end
end

n_rois = size(data_struct.plot_ct,1);
plot_x = (1:4*main_sr)/main_sr-1;

% get component
pk_dp_window_s = common_functions.DA_get_comp_timewindow();
pk_dp_window_frame = (pk_dp_window_s+1)*main_sr;
control_frame = 1:round(0.5*main_sr);
control_factor = 2; % 2 std

tas = struct;
for pi = 1:length(pnames)
    pname = pnames(pi);
    for ci=1:2
        c="cue"+ci;
        tas.(pname).(c).pk_loc = nan(2,n_rois,3);
        tas.(pname).(c).single = cell(2,n_rois);

        this_roi = data_struct.(pname).("cue"+ci);
        this_mus = cellfun(@mean,this_roi);
        for r=1:n_rois
            this_mu = this_mus(:,r);
            this_act = permute(cell2mat(this_roi(:,r)),[1,3,2]);
            control = control_factor*std(this_mu(control_frame));

            % peak
            tmp = pk_dp_window_frame(1,1):pk_dp_window_frame(1,2);
            [pk,loc] = findpeaks(this_mu(tmp),tmp,NPeaks=1,SortStr='descend');
            if isempty(pk)
                pk=nan;loc=nan;
                single = [];
            else
                single = this_act(loc,:);
            end
            sig = pk>control;
            tas.(pname).(c).pk_loc(1,r,:) = [pk;loc;sig];
            tas.(pname).(c).single{1,r} = single;
            % dip
            tmp = pk_dp_window_frame(2,1):pk_dp_window_frame(2,2);
            [pk,loc] = findpeaks(-this_mu(tmp),tmp,NPeaks=1,SortStr='descend');
            pk = -pk;
            if isempty(pk)
                pk=nan;loc=nan;
                single = [];
            else
                single = this_act(loc,:);
            end
            sig = -pk>control;
            tas.(pname).(c).pk_loc(2,r,:) = [pk;loc;sig];
            tas.(pname).(c).single{2,r} = single;
        end
    end
end
save([cf,'processed_and_organized_data\DA_tri_avg_single_filtered_aDMs_only.mat'],'-struct',"tas");


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% Dopamine data: building ITI licking aligned data for DA mice
close all;clear;clc;
cf = [pwd,'\'];
CT = readtable([cf,'raw_data\CT_across_DLXX_mice.xlsx']);
mouse_names_ori = ["DL18","DL20","DL21","DL23"];

% hyper parameters
fc_loading_fft_1 = [0.3,nan];
lick_onset_gap_second = 2; % second
after_lick_window_s = 0.5;
after_lick_ratio_threshold = 0;

iti_remove_cue_s = [-0.5,1];
iti_remove_rew_s = [-0.5,6];

vel_type = "lin_vel";

el_text = ["pre","post"];
% build data
output_data = struct;
for mouse_name = mouse_names_ori
    this_info = common_functions.get_training_info_DA(mouse_name);
    early_late = common_functions.get_include_days_DA(mouse_name);
    
    for el = 1:2 % loop early late
        this_days = cell2mat(early_late(el_text(el)));
        out_data = cell(length(this_days),4);
        for i = 1:length(this_days)
            di = this_days(i);
            ROIBehav = ROI_loader_DA([cf,'raw_data/',char(mouse_name),'/'],di);
            ROI470 = ROIBehav.roi; behav7 = ROIBehav.behavior; this_sr = behav_functions.get_sr(behav7);
            % maybe filter ROI 470
            ROI470_filtered = frequency_pass(ROI470,fc_loading_fft_1,this_sr);

            loco = common_functions.get_overall_vel_acc(behav7,...
                lowpass_freq=1.5, smooth_window=0.5, filter_flag=3,vel_type=vel_type);
            tmp_ta_input = cat(2,behav7.lick_count,loco.vel,loco.acc,ROI470.Fc);
            tmp_ta_input_filtered = cat(2,behav7.lick_count,loco.vel,loco.acc,ROI470_filtered.Fc);

            ITI_bit = get_ITI_pav(behav7,remove_lick=[],remove_cue=iti_remove_cue_s,remove_reward=iti_remove_rew_s);
            % get lick onset
            lick_onsets = get_lick_onset(behav7.lick,lick_onset_gap_second,this_sr,ITI_bit,include_non_ITI=0,...
                after_lick_window_s=after_lick_window_s,after_lick_ratio_threshold=after_lick_ratio_threshold);
            lick_ta = eventTriggeredAverage(tmp_ta_input,lick_onsets,-this_sr+1,this_sr);
            lick_ta_filtered = eventTriggeredAverage(tmp_ta_input_filtered,lick_onsets,-this_sr+1,this_sr);
            out_data{i,1} = this_sr;
            out_data{i,2} = lick_onsets;
            out_data{i,3} = lick_ta.activity;
            out_data{i,4} = lick_ta_filtered.activity;
        end
        out_data = cell2table(out_data,VariableNames=["sr","lick_frame","data","data_filtered"],RowNames="day"+string(this_days));
        mouse_struct.(el_text(el)) = out_data;
    end
    output_data.(mouse_name) = mouse_struct;
end

cwa = struct;
% merge_table into activity struct
for mouse_name = mouse_names_ori
    for el = 1:2 % loop early late
        tbl = output_data.(mouse_name).(el_text(el));
        n_day = size(tbl,1);
        act = []; act_filtered = [];
        for ni = 1:n_day
            act = cat(3,act,cell2mat(tbl{ni,3}));
            act_filtered = cat(3,act_filtered,cell2mat(tbl{ni,4}));
        end
        cwa.non_filtered.(mouse_name).(el_text(el)) = act;
        cwa.filtered.(mouse_name).(el_text(el)) = act_filtered;
    end
end
save([cf,'processed_and_organized_data\DA_ITI_lick_ta_table.mat'],'-struct','output_data','-v7.3');
save([cf,'processed_and_organized_data\DA_ITI_lick_ta.mat'],'-struct','cwa','-v7.3');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% Dopamine data: get component amplitude of licking aligned data for DA mice
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\DA_ITI_lick_ta.mat']).filtered;
mouse_names_ori = string(fields(cwa)');
el_tag = ["pre","post"];
null_std_factor = 2.5;

tas = struct;
for mi = 1:length(mouse_names_ori)
    m_name = mouse_names_ori(mi);
    for el = 1:2
        roi_sr = 30;
        pk_dp_window = round((common_functions.DA_get_comp_timewindow_ITI()+1) * roi_sr);
        pk_window = pk_dp_window(1,1):pk_dp_window(1,2);
        dp_window = pk_dp_window(2,1):pk_dp_window(2,2);
        null_window = 1:round(0.7*roi_sr);
        roi_traces = cwa.(m_name).(el_tag(el));
        roi_traces_mu = mean(roi_traces,3,"omitmissing");
        null_thres = std(roi_traces_mu(null_window,:),[],1,"omitmissing") * null_std_factor;

        [pk_pks,pk_locs_frame,pk_pks_single,pk_sig] = get_peak(roi_traces,roi_traces_mu,pk_window,null_thres);
        pk_locs = pk_locs_frame/roi_sr-1;
        [dp_pks,dp_locs,dp_pks_single,dp_sig] = get_peak(-roi_traces,-roi_traces_mu,dp_window,null_thres);
        dp_locs = dp_locs/roi_sr-1;
        dp_pks = -dp_pks; dp_pks_single = -dp_pks_single;

        single_values = cat(3,pk_pks_single,dp_pks_single);
        single_values = permute(single_values,[3,2,1]);
        mu_loc = cat(1,cat(3,pk_pks,pk_locs),cat(3,dp_pks,dp_locs));
        significance = cat(1,pk_sig,dp_sig);
        tas.(m_name).(el_tag(el)).single_values = single_values;
        tas.(m_name).(el_tag(el)).mu_location_value = mu_loc;
        tas.(m_name).(el_tag(el)).significance = significance;
    end
end
save([cf,'processed_and_organized_data\DA_ITI_lick_tas.mat'],'-struct',"tas");


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% Glu data
close all;clear;clc;
cf = [pwd,'\'];
raw_path = [cf,'raw_data\'];
mouse_names_ori = ["glu924","glu926"];

% hyper parameters
fc_loading_fft = [0.3,nan];
vel_type = "lin_vel";

% build data
output_data = struct;
for mouse_name = mouse_names_ori
    main_sr = 18;
    mouse_path = [raw_path,char(mouse_name),'\'];
    tmp = dir(mouse_path);
    tmp = string({tmp.name});
    dates = tmp(~contains(tmp,[".",".."])); clear("tmp");
    load_ids = 1:length(dates);    
    mouse_struct = struct;
    mouse_struct.include_days = load_ids;
    mouse_struct.data = get_triavg(mouse_path, mouse_name, dates, load_ids, load_ft_fr=[0.3,nan],vel_type=vel_type);
    output_data.(mouse_name) = mouse_struct;
end

cwa = struct;
% merge_table into activity struct
for mouse_name = mouse_names_ori
    main_sr = 18;
    early_late = common_functions.get_training_info_Glu(mouse_name);
    phase_names = early_late.keys;
    all_days = output_data.(mouse_name).include_days;
    all_datas = output_data.(mouse_name).data;
    for pi = 1:length(phase_names)
        p_name = phase_names(pi);
        ids = early_late{p_name};
        p_days = find(ismember(all_days,ids));
        merged_datas = struct;
        for di = 1:length(p_days)
            p_day = p_days(di);
            this_data = all_datas.("file"+p_day);
            merged_datas = merge_pav2cue_data(this_data,merged_datas,main_sr);
        end
        cwa.(mouse_name).(p_name) = merged_datas;
    end
end
save([cf,'processed_and_organized_data\Glu_task_ta_table.mat'],'-struct','output_data','-v7.3');
save([cf,'processed_and_organized_data\Glu_task_ta.mat'],'-struct','cwa','-v7.3');





% ███████╗██╗░░░██╗███╗░░██╗░█████╗░████████╗██╗░█████╗░███╗░░██╗░██████╗
% ██╔════╝██║░░░██║████╗░██║██╔══██╗╚══██╔══╝██║██╔══██╗████╗░██║██╔════╝
% █████╗░░██║░░░██║██╔██╗██║██║░░╚═╝░░░██║░░░██║██║░░██║██╔██╗██║╚█████╗░
% ██╔══╝░░██║░░░██║██║╚████║██║░░██╗░░░██║░░░██║██║░░██║██║╚████║░╚═══██╗
% ██║░░░░░╚██████╔╝██║░╚███║╚█████╔╝░░░██║░░░██║╚█████╔╝██║░╚███║██████╔╝
% ╚═╝░░░░░░╚═════╝░╚═╝░░╚══╝░╚════╝░░░░╚═╝░░░╚═╝░╚════╝░╚═╝░░╚══╝╚═════╝░
function ETAstruct = eventTriggeredAverage(data,events,window1,window2,varargin)
    % ETAstruct = eventTriggeredAverage(data,events,window1,window2,varargin);
    %
    % This function calculates the event-triggered activity and averages, as 
    % well as a bootstrapped null distribution.
    % The bootstrap null distribution is calculated by taking repeated samples
    % of the same size as the events, and checking for significance against a
    % percentile interval of the sample means.
    %
    % takes as input 
    %   - data (n x p matrix): n = datapoints; p = diff timeseries (eg ROIs)
    %   - events: a vector the indices of event occurrences
    %   - window1: index of beginning of averaging window 
    %       (<0 means before the event, >0 means after the event)
    %   - window2: index of end of averaging window
    %       (<0 means before the event, >0 means after the event)
    %
    % optional input
    %   - otherExcl: other events to exclude from the bootstrap
    %   - nullDistr: whether you want to estimate a bootstrapped null distrib 
    %                   (0 = no, 1  = yes)
    %   - bootstrapN: how many bootstrapped samples (default 5000)
    %   - bootstrapSig: what %interval for significance (default .95)
    %   - plotResults: whether you want to plot the results (default: 0 = no)
    % 
    %
    % returns a struct (ETAstruct), with fields
    %   - activity: event triggered activity: r x c x n
    %       - r = another pt in the averaging window
    %       - c = separate timeseries, eg, ROIs
    %       - n = event occurrance 
    %       - so if we gave it 2 ROIs, and want a window of 5 before and 
    %         5 after, with 20 events happening, this matrix would be 11x2x20
    %   - mean: the event triggered average: r x c
    %   - std: event-triggered standard deviation: r x c
    %   - events: the indices of the events included in the triggered average
    %   - windowIdx: the relative indices of window
    %   - bootstrap.means: the means from each of the bootstrapped samples 
    %   - bootstrap.meansMu: the means of the sample means
    %   - bootstrap.meansU: the upper of the 95th percentile (or whatever
    %           percentile interval)
    %   - bootstrap.meansL: the lower of the 95th percentile (or whatever
    %           percentile interval)
    %   - bootstrap.meanSig: a boolean indicating where the event-triggered average is
    %              significant,i.e., mean falls outside the significance interval
    %               you set of the values of the null distribution. For
    %               example, if you set this interval to be 95%, then a value
    %               of -1 means the triggered average mean falls below the
    %               2.5th percentile and a value of 1 means that your triggered
    %               average mean falls above the 97.5th percentile.
    %   - bootstrap.meanSemSig: a boolean indicating where the event-triggered average is
    %               significant,i.e., mean +/- SEM falls outside the significance interval
    %               you set of the values of the null distribution. For
    %               example, if you set this interval to be 95%, then a value
    %               of -1 means the triggered average mean+sem falls below the
    %               2.5th percentile and a value of 1 means that your triggered
    %               average mean-sem falls above the 97.5th percentile.
    %   - bootstrap.sampleable: the indices of parts of the timeseries that can
    %               be sampled for the null distribution (i.e., not the
    %               excluded timepoints)
    %
    % Mai-Anh Vu, 4/30/20
    % updated 5/15/20
    % updated 8/19/20 to check data size for single array data
    % updated 1/21/21 by Mai-Anh to do better bootstrap
    % updated 3/23/21 by Mai-Anh to allow option to plot output
    % updated 3/25/21 by Mai-Anh: fixed commenting above
    % updated 7/11/21 by Mai-Anh: option of saving out figures
    % updated 9/27/22 by Mai-Anh: saves out the indices that are sampleable

    % modified version based on Mai-Anh's original code by Zack


    %%%  parse optional inputs %%%
    ip = inputParser;
    ip.addParameter('otherExcl',[]);
    ip.addParameter('nullDistr',0);
    ip.addParameter('bootstrapN',5000);
    ip.addParameter('bootstrapSig',.95);
    ip.addParameter('saveDir',[]);
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % some setup
    if size(events,2)>size(events,1)
        events = events';
    end
    if size(data,1)==1 
        if size(data,2)>size(data,1)
            data = transpose(data);
        else
            disp('error: data is only 1 number.')
            return
        end
    end
    % output
    ETAstruct = struct;
    % directly return if isempty(event)
    if isempty(events)
        ETAstruct.activity = [];
        ETAstruct.mean = [];
        ETAstruct.std = [];
        warning("Input event is empty.")
    elseif all(isnan(events))
        ETAstruct.activity = nan(window2-window1+1,size(data,2),1);
        ETAstruct.mean = nan(window2-window1+1,size(data,2));
        ETAstruct.std = nan(window2-window1+1,size(data,2));
        warning("Input event is all nans.")
        return
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get event-related activity: index matrix strategy to avoid looping 
    % the row indices (timepoints) of our data matrix
    windowmat = repmat((window1:window2)',1,size(data,2),numel(events));
    
    eventmat = permute(repmat(events,1,size(data,2),window2-window1+1),[3 2 1]);
    
    % disp(size(windowmat))
    % disp(size(eventmat))
    
    if isempty(windowmat) && isempty(eventmat)
        eventmat = windowmat;
    end
    
    rmat = windowmat+eventmat; 
    rnan = rmat<1 | rmat>size(data,1); % indices falling outside our data
    rmat(rnan) = 1;
    
    % the column indices
    cmat = repmat(1:size(data,2),window2-window1+1,1,numel(events));
    % get linear indices
    linmat = sub2ind(size(data),rmat,cmat);
    % get the activity
    try
        ETAstruct.activity = data(linmat);
        ETAstruct.activity(rnan) = nan;
    catch err
        keyboard()
    end
    
    % average
    ETAstruct.mean = nanmean(ETAstruct.activity,3);
    % standard deviation
    ETAstruct.std = nanstd(ETAstruct.activity,[],3);
    
    % the window indices
    ETAstruct.events = events;
    ETAstruct.windowIdx = transpose(window1:window2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bootstrapped null distribution
    
    if nullDistr==1
        % some setup
        ETAstruct.bootstrap.means = nan(window2-window1+1,size(data,2),bootstrapN);
        % excluded events
        if ~isempty(otherExcl) && size(otherExcl,2)>size(otherExcl,1)
            otherExcl = otherExcl';
        end
        excludedEvents = unique([events; otherExcl]);
        includedEvents = (max([0 -window1])+1):(size(data,1)-window2);
        includedEvents = setdiff(includedEvents,excludedEvents);
        includedEvents = includedEvents(includedEvents>0);
        ETAstruct.bootstrap.sampleable = includedEvents;
        ETAstruct.bootstrap.sampleable = ETAstruct.bootstrap.sampleable(...
            ETAstruct.bootstrap.sampleable>(window2-window1+1));
        ETAstruct.bootstrap.sampleable = ETAstruct.bootstrap.sampleable(...
            ETAstruct.bootstrap.sampleable<(size(data,1)-window2+window1));
        
    
        for i = 1:bootstrapN
            % resample from our includedEvents with replacement
            bsevents = includedEvents(randi(numel(includedEvents),numel(events),1))';
            % the row indices (timepoints) of our data matrix
            windowmat = repmat((window1:window2)',1,size(data,2),numel(bsevents));
            eventmat = permute(repmat(bsevents,1,size(data,2),window2-window1+1),[3 2 1]);
            rmat = windowmat+eventmat; 
            % the column indices
            cmat = repmat(1:size(data,2),window2-window1+1,1,numel(bsevents));
            % get linear indices
            linmat = sub2ind(size(data),rmat,cmat);
            % activity
            temp = data(linmat);
            ETAstruct.bootstrap.means(:,:,i) = mean(temp,3);
        end
        % bootstrap average, , and upper and lower 95% percentiles
        ETAstruct.bootstrap.meansMu = nanmean(ETAstruct.bootstrap.means,3);
        ETAstruct.bootstrap.meansU = quantile(ETAstruct.bootstrap.means,1-(1-bootstrapSig)/2,3);
        ETAstruct.bootstrap.meansL = quantile(ETAstruct.bootstrap.means,(1-bootstrapSig)/2,3);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % significant difference: -1 for sig neg, +1 for sig pos
        ETAstruct.bootstrap.meanSig = zeros(size(ETAstruct.mean));
        ETAstruct.bootstrap.meanSig(ETAstruct.mean<ETAstruct.bootstrap.meansL)=-1;
        ETAstruct.bootstrap.meanSig(ETAstruct.mean>ETAstruct.bootstrap.meansU)=1;
        % significant difference: -1 for sig neg, +1 for sig pos
        ETAstruct.bootstrap.meanSemSig = zeros(size(ETAstruct.mean));
        ETAstruct.bootstrap.meanSemSig(ETAstruct.mean+ETAstruct.std/sqrt(numel(events))<ETAstruct.bootstrap.meansL)=-1;
        ETAstruct.bootstrap.meanSemSig(ETAstruct.mean-ETAstruct.std/sqrt(numel(events))>ETAstruct.bootstrap.meansU)=1;
    end
end

function [ROI405,ROI470,behav5,behav7,sample_rate,tag_405470570,tag_470570_aligned] = ROI_loader(path,varargin)
    % load ROI and behav, same as ROI_loader() below but automatically scan for recording type info
    %   INPUT: id is index for Simul, which matches with training day
    Fc_algo = ["exp","exp"];
    Fc_exp_centering = 1;
    Fourier = [nan,nan];
    Fourier_1 = [nan,nan];
    flag_completing_reward = 1;

    ip = inputParser();
    ip.addParameter('Fc_algo',["exp","exp"]);
    ip.addParameter('Fc_exp_centering',1);
    ip.addParameter('Fourier',[nan,nan]); % ft for ROI470
    ip.addParameter('Fourier_1',[nan,nan]); % ft for ROI405 (usually not used because no ft for 405 0r 570)
    ip.addParameter('flag_completing_reward',1);
    parse(ip,varargin{:})
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    dirs = dir(path);
    [behav_405_470,ROI_405_470,tag_405470570,tag_sync,tag_470570_aligned] = find_behav_tif_smart(path);

    behav5 = path +"\"+ dirs(behav_405_470(1)).name;                        
    behav7 = path +"\"+ dirs(behav_405_470(2)).name;
    ROI405 = path +"\"+ dirs(ROI_405_470(1)).name;
    ROI470 = path +"\"+ dirs(ROI_405_470(2)).name;
    ROI405 = load(ROI405);
    ROI470 = load(ROI470);
    behav5 = load(behav5);
    behav7 = load(behav7);
    sample_rate = round(1/(mean(diff(behav7.timestamp),'omitnan')));
    % preprocess                 
    [ROI405,ROI470,behav5,behav7] = preprocess(ROI405,ROI470,behav5,behav7,tag_405470570,tag_sync,tag_470570_aligned,flag_completing_reward);

    fr = sample_rate;
    % if Fc algorithm is exponential, use FtoFc_exp
    if ip.Results.Fc_algo(1) == "exp"
        [Fc1, scale1] = common_functions.FtoFc_exp(ROI405.F,centering=Fc_exp_centering);
        ROI405.Fc = Fc1;
        ROI405.F_baseline = scale1;
    end
    if ip.Results.Fc_algo(2) == "exp"
        [Fc2, scale2] = common_functions.FtoFc_exp(ROI470.F,centering=Fc_exp_centering);
        ROI470.Fc = Fc2;
        ROI470.F_baseline = scale2;
    end

    % if specified then highlowbandpass
    Fourier_flag = any([Fourier(1),Fourier(2)]);
    if Fourier_flag
        padsize7 = [round(size(ROI470.Fc,1)/2),0]; % pad data before any filtering with half_length symmetric
        Fc_padded7 = padarray(ROI470.Fc,padsize7,"symmetric","both");
        if ~isnan(Fourier(1)) && ~isnan(Fourier(2)) % bandpass
            Fc_filtered7 = bandpass(Fc_padded7,ip.Results.Fourier,fr);
        elseif ~isnan(Fourier(1)) && isnan(Fourier(2)) % highpass
            Fc_filtered7 = highpass(Fc_padded7,ip.Results.Fourier(1),fr);
        elseif isnan(Fourier(1)) && ~isnan(Fourier(2)) % lowpass
            Fc_filtered7 = lowpass(Fc_padded7,ip.Results.Fourier(2),fr);
        else
            Fc_filtered7 = Fc_padded7;
        end
        ROI470.Fc = Fc_filtered7(padsize7(1)+1:end-padsize7(1),:);
    end

    % helper functions
    function [behav_1_2,ROI_1_2,tag_405470570,tag_sync,tag_470570_aligned] = find_behav_tif_smart(PATH)
        % find behav and tif at given path
        datainfo_dict = [];
    
        tag_405470570 = nan;
        tag_sync = nan;
        tag_470570_aligned = "NA";
    
        dirs = dir(PATH);
        behav_1_2 = [0,0];
        ROI_1_2 = [0,0];
        
        possible_470_behav = [];
        possible_570_behav = [];
        possible_470_behav_aligned = [];
        possible_570_behav_aligned = [];
        possible_405470_behav_405 = [];
        possible_405470_behav_470 = [];
    
        possible_470_roi = [];
        possible_570_roi = [];
        possible_470_roi_aligned = [];
        possible_570_roi_aligned = [];
    
        virmen_date_pattern = "\d{4}.\d{2}.\d{2}_\d{2}.\d{2}.\d{2}_ttlIn\d{1}";
    
        regexp470_beh = "\d+_ttlIn1_movie1.mat";
        regexp570_beh = "\d+_ttlIn2_movie1.mat";
        regexp470_beh_aligned = "(aligned_ttlIn1_movie1|ttlIn1_movie1_aligned).mat";
        regexp570_beh_aligned = "(aligned_ttlIn2_movie1|ttlIn2_movie1_aligned).mat";
        regexp405470_beh_405 = "\d+_ttlIn1_movie1_405.mat";
        regexp405470_beh_470 = "\d+_ttlIn1_movie1_470.mat";
        
        regexp470_roi = "^(?!.*aligned)(D|d)ata[\w.]+_crop_MC_ROIs.mat"; % note this excludes aligned 470/570
        regexp570_roi = "^(?!.*aligned)R[\w.]+_crop_MC_ROIs.mat"; % note this excludes aligned 470/570
        regexp470_roi_aligned = "^(D|d)ata[\w.]+(aligned_crop_MC_ROIs|crop_MC_ROIs_aligned).mat";
        regexp570_roi_aligned = "^R[\w.]+(aligned_crop_MC_ROIs|crop_MC_ROIs_aligned).mat";
        regexp405470_roi_405 = "^(D|d)ata[\w.]+405[\w.]*_crop_MC_ROIs.mat";
        regexp405470_roi_470 = "^(D|d)ata[\w.]+470[\w.]*_crop_MC_ROIs.mat";
    
        % scan for all possible files
        for i = 1:length(dirs)
            % get behav
            if regexp(string(dirs(i).name),regexp570_beh,"once")
                possible_570_behav = cat(2,possible_570_behav,i);
            elseif regexp(string(dirs(i).name),regexp470_beh,"once")
                possible_470_behav = cat(2,possible_470_behav,i);
            elseif regexp(string(dirs(i).name),regexp570_beh_aligned,"once")
                possible_570_behav_aligned = cat(2,possible_570_behav_aligned,i);
            elseif regexp(string(dirs(i).name),regexp470_beh_aligned,"once")
                possible_470_behav_aligned = cat(2,possible_470_behav_aligned,i);
            elseif regexp(string(dirs(i).name),regexp405470_beh_405,"once")
                possible_405470_behav_405 = cat(2,possible_405470_behav_405,i);
            elseif regexp(string(dirs(i).name),regexp405470_beh_470,"once")
                possible_405470_behav_470 = cat(2,possible_405470_behav_470,i);
    
            % get tif
            elseif regexp(string(dirs(i).name),regexp470_roi,"once")
                possible_470_roi = cat(2,possible_470_roi,i);
            elseif regexp(string(dirs(i).name),regexp570_roi,"once")
                possible_570_roi = cat(2,possible_570_roi,i);
            elseif regexp(string(dirs(i).name),regexp470_roi_aligned,"once")
                possible_470_roi_aligned = cat(2,possible_470_roi_aligned,i);
            elseif regexp(string(dirs(i).name),regexp570_roi_aligned,"once")
                possible_570_roi_aligned = cat(2,possible_570_roi_aligned,i);
            end
        end
    
        % get recording info from scan result
        if length(possible_405470_behav_405) == 1 && length(possible_405470_behav_470) == 1 % (405+470)/(sync)
            tag_405470570 = "405470";
            tag_sync = "sync";
            behav_1_2(1) = possible_405470_behav_405;
            behav_1_2(2) = possible_405470_behav_470;
            if regexp(string(dirs(possible_470_roi(1)).name),regexp405470_roi_405)
                ROI_1_2(1) = possible_470_roi(1);
                ROI_1_2(2) = possible_470_roi(2);
            elseif regexp(string(dirs(possible_470_roi(1)).name),regexp405470_roi_470)
                ROI_1_2(1) = possible_470_roi(2);
                ROI_1_2(2) = possible_470_roi(1);
            else
                error("IfElseError.")
            end
            
        elseif length(possible_470_behav_aligned) == 1 && length(possible_570_behav_aligned) == 1 % (405|470+570)/(sync+aligned)
            % actually use data_info (if provided) to tell if it is 405 or 470
            tag_405470570 = "470570";
            
            if ~isempty(datainfo_dict) 
                this_wavelength = str2double(strsplit(datainfo_dict("wavelength"),"-"));
                if length(this_wavelength) == 1
                    tag_405470570 = "405570";
                end
            end
            tag_sync = "sync";
            tag_470570_aligned = "aligned";
            behav_1_2(2) = possible_470_behav_aligned;
            behav_1_2(1) = possible_570_behav_aligned;
            ROI_1_2(2) = possible_470_roi_aligned;
            ROI_1_2(1) = possible_570_roi_aligned;
    
        elseif length(possible_470_behav) == 1 && length(possible_570_behav) == 1 % (405|470+570)/(sync|async)
            % actually use data_info (if provided) to tell if it is 405 or 470
            assert(isempty(possible_405470_behav_405) && isempty(possible_405470_behav_470),"More than one recording detected.")
            file_name_470 = dirs(possible_470_behav).name;
            file_name_570 = dirs(possible_570_behav).name;
            timestamp_470 = regexp(string(file_name_470),virmen_date_pattern,"once");
            timestamp_570 = regexp(string(file_name_570),virmen_date_pattern,"once");
            timestamp_470 = file_name_470(timestamp_470:timestamp_470+18);
            timestamp_570 = file_name_570(timestamp_570:timestamp_570+18);
            tag_405470570 = "470570";
            
            if ~isempty(datainfo_dict)
                this_wavelength = str2double(strsplit(datainfo_dict("wavelength"),"-"));
                if length(this_wavelength) == 1
                    tag_405470570 = "405570";
                end
            end
            if strcmp(timestamp_470,timestamp_570)
                tag_sync = "sync";
                tag_470570_aligned = "unaligned";
            else
                tag_sync = "async"; % this is going to be changed, put each session in a single folder
            end
            behav_1_2(2) = possible_470_behav;
            behav_1_2(1) = possible_570_behav;
            ROI_1_2(2) = possible_470_roi;
            ROI_1_2(1) = possible_570_roi;
    
        elseif length(possible_470_behav) == 1
            assert(isempty(possible_405470_behav_405) && isempty(possible_405470_behav_470)&&isempty(possible_570_behav),"More than one recording detected.")
            tag_405470570 = "470";
            tag_sync = "NA";
            behav_1_2(2) = possible_470_behav;
            ROI_1_2(2) = possible_470_roi;
        elseif length(possible_570_behav) == 1
            assert(isempty(possible_405470_behav_405) && isempty(possible_405470_behav_470)&&isempty(possible_470_behav),"More than one recording detected.")
            tag_405470570 = "570";
            tag_sync = "NA";
            behav_1_2(1) = possible_570_behav;
            ROI_1_2(1) = possible_570_roi;
        else
            error("This folder: "+string(PATH)+" doesn't seem to have standard organization.")
        end
    end

    function [ROI405,ROI470,behav5,behav7] = preprocess(ROI405,ROI470,behav5,behav7,tag_405470570,tag_sync,tag_470570_aligned,completing_reward)
    % add missing field, truncate frame, remove nan
        if ~isfield(ROI405,'F_baseline')
            ROI405.F_baseline = zeros(size(ROI405.Fc));
        end
        if ~isfield(ROI405,'F_baseline')
            ROI470.F_baseline = zeros(size(ROI470.Fc));
        end
        fields_5 = fieldnames(behav5);
        fields_7 = fieldnames(behav7);
        behav5_missing = reshape(string(setdiff(fields_7,fields_5)),1,[]);
        behav7_missing = reshape(string(setdiff(fields_5,fields_7)),1,[]);
        for field_name = behav5_missing
            behav5.(field_name) = "missing";
        end
        for field_name = behav7_missing
            behav7.(field_name) = "missing";
        end
        behav5_missing_joined = strjoin(behav5_missing," ");
        behav7_missing_joined = strjoin(behav7_missing," ");
        if ~isempty(behav5_missing)
            fprintf("Adding "+behav5_missing_joined+" to behav1.\n")
        end
        if ~isempty(behav7_missing)
            fprintf("Adding "+behav7_missing_joined+" to behav2.\n")
        end
        % get truncate frame number
        fields = fieldnames(behav5);
        mask = ones(1,length(fields));
        for i = 1:length(fields)
            if contains(string(fields{i}),"ignore") || contains(string(fields{i}),"LED")
                mask(i) = 0;
            end
        end
        fields = fields(logical(mask));

        frame_numbers = inf(1,2*length(fields)+6);
        nan_detector = 0;
        
        frame_numbers(2*length(fields)+1) = size(ROI405.F,1);
        frame_numbers(2*length(fields)+2) = size(ROI470.F,1);
        frame_numbers(2*length(fields)+3) = size(ROI405.F_baseline,1);
        frame_numbers(2*length(fields)+4) = size(ROI470.F_baseline,1);
        frame_numbers(2*length(fields)+5) = size(ROI405.Fc,1);
        frame_numbers(2*length(fields)+6) = size(ROI470.Fc,1);
        for i = 1:length(fields)
            if isnumeric(behav7.(fields{i})) && isnumeric(behav5.(fields{i}))
                if size(behav5.(fields{i}),2) > size(behav5.(fields{i}),1)
                    behav5.(fields{i}) = behav5.(fields{i})';
                end
                if size(behav7.(fields{i}),2) > size(behav7.(fields{i}),1)
                    behav7.(fields{i}) = behav7.(fields{i})';
                end
                frame_numbers(2*i-1) = size(behav5.(fields{i}),1);
                frame_numbers(2*i) = size(behav7.(fields{i}),1);
            end
        end
        
        % truncate and remove nan
        if ~(tag_sync == "async" || tag_470570_aligned == "unaligned")
            truncate_frame = min(frame_numbers);
            for i = 1:length(fields)
                if isnumeric(behav7.(fields{i})) && isnumeric(behav5.(fields{i})) % do not remove nan from behav
                    behav5.(fields{i}) = behav5.(fields{i})(1:truncate_frame,:);
                    behav7.(fields{i}) = behav7.(fields{i})(1:truncate_frame,:);
                end
            end
            ROI405.F = ROI405.F(1:truncate_frame,:);
            ROI470.F = ROI470.F(1:truncate_frame,:);
            ROI405.F_baseline = ROI405.F_baseline(1:truncate_frame,:);
            ROI470.F_baseline = ROI470.F_baseline(1:truncate_frame,:);
            ROI405.Fc = ROI405.Fc(1:truncate_frame,:);
            ROI470.Fc = ROI470.Fc(1:truncate_frame,:);
        end
        
        if any(isnan(ROI405.F),'all') || any(isnan(ROI470.F),'all') || any(isnan(ROI405.F_baseline),'all') ...
                || any(isnan(ROI470.F_baseline),'all') || any(isnan(ROI405.Fc),'all') || any(isnan(ROI470.Fc),'all')
            nan_detector = 1;
            ROI405.F(isnan(ROI405.F)) = 0;
            ROI470.F(isnan(ROI470.F)) = 0;
            ROI405.F_baseline(isnan(ROI405.F_baseline)) = 0;
            ROI470.F_baseline(isnan(ROI470.F_baseline)) = 0;
            ROI405.Fc(isnan(ROI405.Fc)) = 0;
            ROI470.Fc(isnan(ROI470.Fc)) = 0;
        end

        if completing_reward && (tag_405470570 == "405470" || (tag_405470570 == "470570" && tag_470570_aligned == "unaligned"))
            [behav5,behav7] = ultimate_mice_class.completing_reward(behav5,behav7);
            fprintf("Completing rewards.\n")
        end
        
        if nan_detector == 1
            fprintf("\n**************************************\n* Nan in behav or roi, replaced by 0 *\n**************************************\n\n")
        end
    end

end

function ta = get_triavg(path, mouse, dates, ids, varargin)
    cue2Tone_shift_forward_flag=1;
    rewTrgConsumption=1;
    get_405=0;
    load_ft_fr=[nan,nan];
    load_ft_fr_1=[nan,nan];
    GLM = [];
    vel_type = "total_vel_1";

    ip = inputParser;
    ip.addParameter('cue2Tone_shift_forward_flag',1);
    ip.addParameter('rewTrgConsumption',1);
    ip.addParameter('get_405',0)
    ip.addParameter('load_ft_fr',[nan,nan])
    ip.addParameter('load_ft_fr_1',[nan,nan])
    ip.addParameter('GLM',[])
    ip.addParameter('vel_type',"total_vel_1")
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    %%% input modification %%%
    n_days = length(ids);
    if contains(mouse,["G12","G15"])
        main_sr = 30;
    else
        main_sr = 18;
    end
    ROI_struct = struct;
    behav_struct = struct;
    for i = 1:n_days
        id = ids(i);
        di = "day"+string(id);
        this_path = [path,char(dates(id)),'\'];
        [ROI405i,ROI470i,behav5i,behav7i] = ROI_loader(this_path, Fourier = load_ft_fr, Fourier_1 = load_ft_fr_1);

        % if GML is provided, regress things out
        if ~isempty(GLM)
            reduced_fc = ROI470i.Fc;
            n_rois = size(reduced_fc,2);
            n_frames = size(reduced_fc,1);
            sr = common_functions.get_sr(behav7i);
            flag_session_with_405=0;
            loco = common_functions.get_overall_vel_acc(behav7i,...
                    'lowpass_freq',1.5,...
                    'smooth_window',0.5,...
                    'filter_flag',3,...
                    'vel_type',vel_type);
            behav_delay_value_ori = cell2mat(table2cell(GLM.xcorr(1,5:8))')';
            behav_delay_value = round(behav_delay_value_ori/main_sr*sr);
            behav_delay_name = GLM.xcorr(1,5:8).Properties.VariableNames;
            behav_delay_name = cellfun(@(x) string(x([1:4,9:11])),behav_delay_name);
            behav_delay_table = array2table(behav_delay_value,VariableNames=behav_delay_name);
            for r = 1:n_rois
                this_model = GLM.mdls{r};
                this_glm_input_value = [loco.vel,loco.acc];
                this_glm_input_name = ["vel_zero","acc_zero"];
                this_delay = behav_delay_table{r,:}; this_delay(this_delay==0) = nan;
                shift = {min(this_delay(this_delay<0),[],"omitmissing"),max(this_delay(this_delay>0),[],"omitmissing")};
                fc_mask = [1,n_frames];
                if ~isempty(shift{1})
                    fc_mask(1) = 1-shift{1};
                end
                if ~isempty(shift{2})
                    fc_mask(2) = n_frames - shift{2};
                end
                n_frame_truncated = fc_mask(2) - fc_mask(1) + 1;
                this_glm_input_value = this_glm_input_value(fc_mask(1):fc_mask(2),:);

                for behav_lags_i = 1:length(this_delay)
                    if ~isnan(this_delay(behav_lags_i))
                        if behav_lags_i<=2
                            this_behav = loco.vel;
                        else
                            this_behav = loco.acc;
                        end
                        this_glm_input_value = cat(2,this_glm_input_value,...
                            this_behav(fc_mask(1)+this_delay(behav_lags_i):fc_mask(1)+this_delay(behav_lags_i)+n_frame_truncated-1));
                        this_glm_input_name = cat(2,this_glm_input_name,behav_delay_name(behav_lags_i));
                    end
                end
                coef_values = this_model.Coefficients{this_glm_input_name,"Estimate"};
                to_be_regressedout = this_glm_input_value*coef_values + this_model.Coefficients{"(Intercept)","Estimate"};
                reduced_fc(fc_mask(1):fc_mask(2),r) = reduced_fc(fc_mask(1):fc_mask(2),r) - to_be_regressedout;
                ROI470i.Fc = reduced_fc;
            end
        end
        ROI_struct.(di).ROI1 = ROI405i; % well I guess this one is not needed anyway
        ROI_struct.(di).ROI2 = ROI470i;
        behav_struct.(di).behav1 = behav5i;
        behav_struct.(di).behav2 = behav7i;
    end

    %%% COMPILE CUE INFO FIRST %%%
    cueInfo = struct;
    cueInfo405 = struct;
    for f = 1:n_days
        di = "day"+string(f);
        behav1 = behav_struct.(di).behav1;
        behav2 = behav_struct.(di).behav2;

        cuetmp = common_functions.pav2cue_getTrialTimes(behav1,cue2Tone_shift_forward_flag=cue2Tone_shift_forward_flag);
        cueInfo405.(['file' num2str(f)]) = cuetmp;
    
        cuetmp = common_functions.pav2cue_getTrialTimes(behav2,cue2Tone_shift_forward_flag=cue2Tone_shift_forward_flag);
        cueInfo.(['file' num2str(f)]) = cuetmp;
    end

    for f = 1:n_days
        cueInfo.(['file' num2str(f)]).trialTimes = cueInfo.(['file' num2str(f)]).cueEnd-cueInfo.(['file' num2str(f)]).cueStart;
    end
    
    
    %%% TRG AVG: EPOCHS %%%
    ta = struct;
    time_before_event = 1; %get also 1 second before cue onset
    % loop over file list
    for f = 1:n_days
        di = "day"+string(f);
        % load roi
        if ~get_405
            behav = behav_struct.(di).behav2;
            roi = ROI_struct.(di).ROI2;
        elseif get_405
            behav = behav_struct.(di).behav1;
            roi = ROI_struct.(di).ROI1;
        end
        loco = common_functions.ball2xy(behav);
        roi.Fc = [behav.lick_count(1:size(roi.Fc,1)),...
            loco.linear_velocity(1:size(roi.Fc,1)),...
            loco.angular_velocity(1:size(roi.Fc,1)),...
            roi.Fc];

        ta.(['file' num2str(f)]).sample_rate = common_functions.get_sr(behav);
        % get triggered averages
        for cue = 1:2
            % cue On
            tmp = eventTriggeredAverage(roi.Fc,cueInfo.(['file' num2str(f)]).cueStart(cueInfo.(['file' num2str(f)]).cue==cue & ~isnan(cueInfo.(['file' num2str(f)]).cueStart)),-time_before_event*cueInfo.(['file' num2str(f)]).fr,cueInfo.(['file' num2str(f)]).fr*3-1,'nullDistr',0);
            ta.(['file' num2str(f)]).cueOn.(['cue' num2str(cue)]).activity = tmp.activity;
            ta.(['file' num2str(f)]).cueOn.(['cue' num2str(cue)]).mu = tmp.mean;
            ta.(['file' num2str(f)]).cueOn.(['cue' num2str(cue)]).sem = tmp.std/sqrt(size(tmp.activity,3));
            % rew onset                    
            tmp = eventTriggeredAverage(roi.Fc,cueInfo.(['file' num2str(f)]).rewOn(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==1 & ~isnan(cueInfo.(['file' num2str(f)]).rewOn)),-time_before_event*cueInfo.(['file' num2str(f)]).fr,cueInfo.(['file' num2str(f)]).fr*3-1,'nullDistr',0);
            ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).rew.activity = tmp.activity;
            ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).rew.mu = tmp.mean;
            ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).rew.sem = tmp.std/sqrt(size(tmp.activity,3));
            if rewTrgConsumption == 1
                tmp = eventTriggeredAverage(roi.Fc,cueInfo.(['file' num2str(f)]).rewLick(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==1 & ~isnan(cueInfo.(['file' num2str(f)]).rewLick)),-time_before_event*cueInfo.(['file' num2str(f)]).fr,cueInfo.(['file' num2str(f)]).fr*3-1,'nullDistr',0);
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).rew_consump.activity = tmp.activity;
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).rew_consump.mu = tmp.mean;
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).rew_consump.sem = tmp.std/sqrt(size(tmp.activity,3));
            end
            % omission onset
            if sum(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==0) > 0
                tmp = eventTriggeredAverage(roi.Fc,cueInfo.(['file' num2str(f)]).cueStart(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==0 & ~isnan(cueInfo.(['file' num2str(f)]).cueStart))+ round(nanmean(cueInfo.(['file' num2str(f)]).rewOn(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==1)-cueInfo.(['file' num2str(f)]).cueStart(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==1))),-time_before_event*cueInfo.(['file' num2str(f)]).fr,cueInfo.(['file' num2str(f)]).fr*3-1,'nullDistr',0);
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit.activity = tmp.activity;
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit.mu = tmp.mean;
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit.sem = tmp.std/sqrt(size(tmp.activity,3));
                if rewTrgConsumption == 1
                    tmp = eventTriggeredAverage(roi.Fc,cueInfo.(['file' num2str(f)]).omitLick(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==0 & ~isnan(cueInfo.(['file' num2str(f)]).omitLick)),-time_before_event*cueInfo.(['file' num2str(f)]).fr,cueInfo.(['file' num2str(f)]).fr*3-1,'nullDistr',0);
                    ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit_consump.activity = tmp.activity;
                    ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit_consump.mu = tmp.mean;
                    ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit_consump.sem = tmp.std/sqrt(size(tmp.activity,3));
                end
            else
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit.activity = [];
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit.mu = [];
                ta.(['file' num2str(f)]).rewOn.(['cue' num2str(cue)]).omit.sem = [];
            end
            % cue off
            rewStr = {'omit','rew'};
            for rew = 0:1
                if sum(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==rew) > 0
                    tmp = eventTriggeredAverage(roi.Fc,cueInfo.(['file' num2str(f)]).cueEnd(cueInfo.(['file' num2str(f)]).cue==cue & cueInfo.(['file' num2str(f)]).rew==rew & ~isnan(cueInfo.(['file' num2str(f)]).cueEnd)),-time_before_event*cueInfo.(['file' num2str(f)]).fr,cueInfo.(['file' num2str(f)]).fr*3-1,'nullDistr',0);
                    ta.(['file' num2str(f)]).cueOff.(['cue' num2str(cue)]).(rewStr{rew+1}).activity = tmp.activity;
                    ta.(['file' num2str(f)]).cueOff.(['cue' num2str(cue)]).(rewStr{rew+1}).mu = tmp.mean;
                    ta.(['file' num2str(f)]).cueOff.(['cue' num2str(cue)]).(rewStr{rew+1}).sem = tmp.std/sqrt(size(tmp.activity,3));
                else
                    ta.(['file' num2str(f)]).cueOff.(['cue' num2str(cue)]).(rewStr{rew+1}).activity = [];
                    ta.(['file' num2str(f)]).cueOff.(['cue' num2str(cue)]).(rewStr{rew+1}).mu = [];
                    ta.(['file' num2str(f)]).cueOff.(['cue' num2str(cue)]).(rewStr{rew+1}).sem = [];
                end
            end
    
        end
        % unpredicted reward
        tmp = eventTriggeredAverage(roi.Fc,cueInfo.(['file' num2str(f)]).unpredRewOn(~isnan(cueInfo.(['file' num2str(f)]).unpredRewOn)),-time_before_event*cueInfo.(['file' num2str(f)]).fr,cueInfo.(['file' num2str(f)]).fr*3-1,'nullDistr',0);
        ta.(['file' num2str(f)]).rewOn.unpred.rew.activity = tmp.activity;
        ta.(['file' num2str(f)]).rewOn.unpred.rew.mu = tmp.mean;
        ta.(['file' num2str(f)]).rewOn.unpred.rew.sem = tmp.std/sqrt(size(tmp.activity,3));
        if rewTrgConsumption == 1
            tmp = eventTriggeredAverage(roi.Fc,cueInfo.(['file' num2str(f)]).unpredRewLick(~isnan(cueInfo.(['file' num2str(f)]).unpredRewLick)),-time_before_event*cueInfo.(['file' num2str(f)]).fr,cueInfo.(['file' num2str(f)]).fr*3-1,'nullDistr',0);
            ta.(['file' num2str(f)]).rewOn.unpred.rew_consump.activity = tmp.activity;
            ta.(['file' num2str(f)]).rewOn.unpred.rew_consump.mu = tmp.mean;
            ta.(['file' num2str(f)]).rewOn.unpred.rew_consump.sem = tmp.std/sqrt(size(tmp.activity,3));
        end
    end
end % get_triavg end

function ITI_bit = get_ITI_pav(behav,sr,varargin)
    % get ITI bit filter for single pav behav file

    % remove whole cue session and 3 second afer unpred reward
    ip = inputParser();
    ip.addParameter('remove_reward',[0,3]);
    ip.addParameter('remove_lick',[-0.5,0.5]);
    ip.addParameter('remove_cue',[])
    ip.parse(varargin{:})
    ip = ip.Results;

    remove_reward = round(ip.remove_reward * sr);
    remove_lick = round(ip.remove_lick * sr);
    remove_cue = round(ip.remove_cue * sr);

    rew_window = remove_reward(2)-remove_reward(1)+1; % 3 seconds
    conv_filter = ones(rew_window,1); 
    reward_bit = behav.reward;
    reward_conv = conv(reward_bit,conv_filter,"full");
    reward_conv = reward_conv(-remove_reward(1)+1:end-remove_reward(2));
    reward_conv = reward_conv~=0;

    if ~isempty(remove_lick)
        lick_window = remove_lick; % 1 second before and 1s after
        conv_filter = ones(lick_window(2)-lick_window(1)+1,1);
        lick_bit = behav.lick;
        lick_conv = conv(lick_bit,conv_filter,"full");
        lick_conv = lick_conv(-lick_window(1)+1:end-lick_window(2));

        lick_conv = lick_conv~=0;
    else
        lick_conv = zeros(size(behav.timestamp));
    end
    
    if ~isempty(remove_cue)
        cue_bit = behav.stimulus_led | behav.stimulus_sound2;
        cue_window = remove_cue;
        conv_filter = ones(cue_window(2)-cue_window(1)+1,1);
        cue_conv = conv(cue_bit,conv_filter,"full");
        cue_conv = cue_conv(-cue_window(1)+1:end-cue_window(2));
        cue_conv = cue_conv~=0;
    else
        cue_conv = zeros(size(behav.timestamp));
    end

    to_remove = reward_conv | cue_conv | lick_conv;

    ITI_bit = ones(size(behav.timestamp));
    ITI_bit(to_remove) = 0;
    ITI_bit = logical(ITI_bit);
end

function lick_onsets = get_lick_onset(lick_frame,gap_threshold_s,sr,valid_frame_mask,varargin)
    % get lick onsets from all licking frames
    ip = inputParser;
    ip.addParameter("include_non_ITI",false)
    ip.addParameter("after_lick_window_s",0.5)
    ip.addParameter("after_lick_ratio_threshold",0.5)
    ip.parse(varargin{:})
    for j = fields(ip.Results)'
        eval([j{1}, '= ip.Results.', j{1}, ';']);
    end
    
    gap_threshold = round(gap_threshold_s * sr);
    after_lick_window = round(after_lick_window_s * sr);

    conv_filter_1 = ones(gap_threshold,1);
    conv_filter_2 = ones(after_lick_window,1);
    conved_1 = conv(lick_frame,conv_filter_1,"full");
    conved_2 = conv(lick_frame,conv_filter_2,"valid");
    conved_2 = conved_2/after_lick_window >= after_lick_ratio_threshold;
    conved_1 = conved_1(1:length(conved_2));

    lick_onsets = find(conved_1==1 & conved_2);
    lick_onsets = lick_onsets(logical(lick_frame(lick_onsets)));
    if ~include_non_ITI
        lick_onsets = lick_onsets(valid_frame_mask(lick_onsets));
    end
end

function output_data = get_triavg_ITI_lick(path, mouse, dates, ids, this_session_info, varargin)
    load_ft_fr=[nan,nan];
    load_ft_fr_1=[nan,nan];
    GLM = [];
    vel_type = "total_vel_1";

    ip = inputParser;
    ip.addParameter('GLM',[])
    ip.addParameter('vel_type',"total_vel_1")
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    lick_onset_gap_second = 2; % second
    after_lick_window_s = 0.5;
    after_lick_ratio_threshold = 0;
    
    iti_remove_cue_s = [-0.5,1];
    iti_remove_rew_s = [-0.5,6];

    output_data = struct;
    output_data.include_days = ids;
    output_data.info_days = this_session_info;
    output_data.data_tags = ["sr","lick_frames","activity"];
    out_data = cell(length(ids),3);

    n_days = length(ids);
    if contains(mouse,["G12","G15"])
        main_sr = 30;
    else
        main_sr = 18;
    end

    for i = 1:n_days
        id = ids(i);
        this_path = [path,char(dates(id)),'\'];
        [ROI405i,ROI470i,behav5i,behav7i] = ROI_loader(this_path, Fourier = load_ft_fr, Fourier_1 = load_ft_fr_1);

        % if GML is provided, regress things out
        if ~isempty(GLM)
            reduced_fc = ROI470i.Fc;
            n_rois = size(reduced_fc,2);
            n_frames = size(reduced_fc,1);
            sr = common_functions.get_sr(behav7i);
            flag_session_with_405=0;
            loco = common_functions.get_overall_vel_acc(behav7i,...
                    'lowpass_freq',1.5,...
                    'smooth_window',0.5,...
                    'filter_flag',3,...
                    'vel_type',vel_type);
            behav_delay_value_ori = cell2mat(table2cell(GLM.xcorr(1,5:8))')';
            behav_delay_value = round(behav_delay_value_ori/main_sr*sr);
            behav_delay_name = GLM.xcorr(1,5:8).Properties.VariableNames;
            behav_delay_name = cellfun(@(x) string(x([1:4,9:11])),behav_delay_name);
            behav_delay_table = array2table(behav_delay_value,VariableNames=behav_delay_name);
            for r = 1:n_rois
                this_model = GLM.mdls{r};
                this_glm_input_value = [loco.vel,loco.acc];
                this_glm_input_name = ["vel_zero","acc_zero"];
                this_delay = behav_delay_table{r,:}; this_delay(this_delay==0) = nan;
                shift = {min(this_delay(this_delay<0),[],"omitmissing"),max(this_delay(this_delay>0),[],"omitmissing")};
                fc_mask = [1,n_frames];
                if ~isempty(shift{1})
                    fc_mask(1) = 1-shift{1};
                end
                if ~isempty(shift{2})
                    fc_mask(2) = n_frames - shift{2};
                end
                n_frame_truncated = fc_mask(2) - fc_mask(1) + 1;
                this_glm_input_value = this_glm_input_value(fc_mask(1):fc_mask(2),:);

                for behav_lags_i = 1:length(this_delay)
                    if ~isnan(this_delay(behav_lags_i))
                        if behav_lags_i<=2
                            this_behav = loco.vel;
                        else
                            this_behav = loco.acc;
                        end
                        this_glm_input_value = cat(2,this_glm_input_value,...
                            this_behav(fc_mask(1)+this_delay(behav_lags_i):fc_mask(1)+this_delay(behav_lags_i)+n_frame_truncated-1));
                        this_glm_input_name = cat(2,this_glm_input_name,behav_delay_name(behav_lags_i));
                    end
                end
                coef_values = this_model.Coefficients{this_glm_input_name,"Estimate"};
                to_be_regressedout = this_glm_input_value*coef_values + this_model.Coefficients{"(Intercept)","Estimate"};
                reduced_fc(fc_mask(1):fc_mask(2),r) = reduced_fc(fc_mask(1):fc_mask(2),r) - to_be_regressedout;
                ROI470i.Fc = reduced_fc;
            end
        end
        sr = behav_functions.get_sr(behav7i);
        mov_data = common_functions.ball2xy(behav7i);
        tmp_ta_input = cat(2,behav7i.lick_count,mov_data.linear_velocity,mov_data.angular_velocity,ROI470i.Fc);
        ITI_bit = get_ITI_pav(behav7,sr,remove_lick=[],remove_cue=iti_remove_cue_s,remove_reward=iti_remove_rew_s);

        % get lick onset during ITI
        lick_onsets = get_lick_onset(behav7i.lick,lick_onset_gap_second,sr,ITI_bit,include_non_ITI=0,...
            after_lick_window_s=after_lick_window_s,after_lick_ratio_threshold=after_lick_ratio_threshold);
        lick_ta = eventTriggeredAverage(tmp_ta_input,lick_onsets,-sr+1,sr);
        out_data{i,1} = sr;
        out_data{i,2} = lick_onsets;
        out_data{i,3} = lick_ta.activity;
    end
    output_data.data = out_data;
end % get_triavg_ITI_lick end

function data_cated = cat_data_across_day(data_in,data_in_id)
    % concat data
    data_cated = data_in(1,:);
    data_cated{1,1} = unique(data_in{data_in_id,1});
    data_cated{1,2} = strjoin(unique(data_in{data_in_id,2}),"_");

    numerical_field_names = ["ITI_bit","vel","acc","470","405","vel_original"];
    tmp = data_in{data_in_id,numerical_field_names};
    for f_name_i = 1:length(numerical_field_names)
        f_name = numerical_field_names(f_name_i);
        data_cated{1,f_name} = {cat(1,tmp{:,f_name_i})};
    end
end

function [I_pos,I_neg,r,lags] = find_lags_using_xcorr(fc,vel,xcorr_maxlag)
    % find max pos lag and min neg lag between fc and vel using xcorr
    [r,lags] = xcorr(fc,vel,xcorr_maxlag(2));
    r_included = r(round(length(lags)/2)-xcorr_maxlag(1):round(length(lags)/2)+xcorr_maxlag(1));
    lags_included = lags(round(length(lags)/2)-xcorr_maxlag(1):round(length(lags)/2)+xcorr_maxlag(1));
    tmp_r = r_included; tmp_r_pos=tmp_r; tmp_r_neg=tmp_r;
    tmp_r_pos(tmp_r<0)=nan;
    tmp_r_neg(tmp_r>0)=nan;
    [~,I_pos] = max(tmp_r_pos,[],"omitmissing");
    [~,I_neg] = min(tmp_r_neg,[],"omitmissing");
    if all(isnan(tmp_r_pos)) || any(I_pos == [1,length(tmp_r_pos)])
        I_pos = nan;
    else
        I_pos = lags_included(I_pos);
    end
    if all(isnan(tmp_r_neg)) || any(I_neg == [1,length(tmp_r_pos)])
        I_neg = nan;
    else
        I_neg = lags_included(I_neg);
    end
end

function ta_total = get_ta_total(ta,ids,mode,varargin)
    % For input mode (actaully just different components cue/rew/unpre etc.), 
    % combine each trials of ta specified by ids into one ta representing a phase including all specified ids
    rescale_factors = '';

    ip = inputParser;
    ip.addParameter('rescale_factors','');
    ip.addParameter('sample_rate',nan); % target sample rate for interpolation if sr doesn't match
    ip.parse(varargin{:})
    for j = fields(ip.Results)'
        eval([j{1}, '= ip.Results.',j{1},';']);
    end
    interpolation_flag = 0;

    % do a sample rate check
    n_days = size(fields(ta),1);
    srs = nan(1,n_days);
    fi = 0;
    for f = fields(ta)'
        fi = fi+1;
        srs(fi) = ta.(f{1}).sample_rate;
    end
    if length(unique(srs))~=1 && isnan(sample_rate)
        error("Tri Avg contains data of different SR, input sample_rate needed for interpolation.")
    end
    % parse mode
    if contains(mode,["rew1","rew2"])
        mode_parsed = "rew12";
    elseif contains(mode,["cue1","cue2"])
        mode_parsed = "cue12";
    else
        mode_parsed = mode;
    end
    ta_total = struct;
    switch mode_parsed
        case "unpred"
            ta_value = [];
            ta_value_consump = [];
            for i = 1:length(ids)
                this_sample_rate = ta.("file"+string(i)).sample_rate;
                % process activity
                activity_to_combine = ta.("file"+string(i)).rewOn.unpred.rew.activity;
                if ~isempty(rescale_factors)
                    for r = 4:size(activity_to_combine,2)
                        temp1 = permute(activity_to_combine(:,r,:),[1,3,2]);
                        temp2 = temp1 * rescale_factors(ids(i),r-3);
                        activity_to_combine(:,r,:) = temp2;
                    end
                end
                if ~isnan(sample_rate) && this_sample_rate ~= sample_rate
                    activity_to_combine = common_functions.interp_ta(activity_to_combine,this_sample_rate,sample_rate);
                    interpolation_flag = interpolation_flag + 1;
                end
                ta_value = cat(3,ta_value,activity_to_combine);

                % process activity for unpred (consump trigger if any)
                if isfield(ta.("file"+string(i)).rewOn.unpred,"rew_consump")
                    unpred_consump_flag = 1;
                    activity_to_combine = ta.("file"+string(i)).rewOn.unpred.rew_consump.activity;
                    if ~isempty(rescale_factors)
                        for r = 4:size(activity_to_combine,2)
                            temp1 = permute(activity_to_combine(:,r,:),[1,3,2]);
                            temp2 = temp1 * rescale_factors(ids(i),r-3);
                            activity_to_combine(:,r,:) = temp2;
                        end
                    end
                    if ~isnan(sample_rate) && this_sample_rate ~= sample_rate
                        activity_to_combine = common_functions.interp_ta(activity_to_combine,this_sample_rate,sample_rate);
                        interpolation_flag = interpolation_flag + 1;
                    end
                    ta_value_consump = cat(3,ta_value_consump,activity_to_combine);
                end
            end
            ta_total.activity = ta_value;
            ta_total.mu = mean(ta_value,3,'omitnan');
            ta_total.sem = std(ta_value,[],3,'omitnan')/sqrt(size(ta_value,3)-1);
            if unpred_consump_flag
                ta_total.activity_consump = ta_value_consump;
                ta_total.mu_consump = mean(ta_value_consump,3,'omitnan');
                ta_total.sem_consump = std(ta_value_consump,[],3,'omitnan')/sqrt(size(ta_value_consump,3)-1);
            end
        case "rew12"
            if mode == "rew1"
                mode_ci = 1;
            elseif mode == "rew2"
                mode_ci = 2;
            else
                error("mode_ci should only be 1 or 2.")
            end
            rew_value = [];
            rew_value_consump = [];
            rew_value_omi = [];
            rew_value_omi_consump = [];
            for i = 1:length(ids)
                this_sample_rate = ta.("file"+string(i)).sample_rate;
                % process activity for rew (onset trigger)
                activity_to_combine = ta.("file"+string(i)).rewOn.("cue"+mode_ci).rew.activity;
                if  ~isempty(rescale_factors)
                    for r = 4:size(activity_to_combine,2)
                        temp1 = permute(activity_to_combine(:,r,:),[1,3,2]);
                        temp2 = temp1 * rescale_factors(ids(i),r-3);
                        activity_to_combine(:,r,:) = temp2;
                    end
                end
                if ~isnan(sample_rate) && this_sample_rate ~= sample_rate
                    activity_to_combine = common_functions.interp_ta(activity_to_combine,this_sample_rate,sample_rate);
                    interpolation_flag = interpolation_flag + 1;
                end
                rew_value = cat(3,rew_value,activity_to_combine);
                
                % process activity for rew (consump trigger if any)
                if isfield(ta.("file"+string(i)).rewOn.("cue"+mode_ci),"rew_consump")
                    rew_consump_flag = 1;
                    activity_to_combine = ta.("file"+string(i)).rewOn.("cue"+mode_ci).rew_consump.activity;
                    if  ~isempty(rescale_factors)
                        for r = 4:size(activity_to_combine,2)
                            temp1 = permute(activity_to_combine(:,r,:),[1,3,2]);
                            temp2 = temp1 * rescale_factors(ids(i),r-3);
                            activity_to_combine(:,r,:) = temp2;
                        end
                    end
                    if ~isnan(sample_rate) && this_sample_rate ~= sample_rate
                        activity_to_combine = common_functions.interp_ta(activity_to_combine,this_sample_rate,sample_rate);
                        interpolation_flag = interpolation_flag + 1;
                    end
                    rew_value_consump = cat(3,rew_value_consump,activity_to_combine);
                end

                if isfield(ta.("file"+string(i)).rewOn.("cue"+mode_ci),"omit") % if has omission
                    omission_flag = 1;
                    activity_to_combine_omi = ta.("file"+string(i)).rewOn.("cue"+mode_ci).omit.activity;
                    if isfield(ta.("file"+string(i)).rewOn.("cue"+mode_ci),"omit_consump")
                        activity_to_combine_omi_consump = ta.("file"+string(i)).rewOn.("cue"+mode_ci).omit_consump.activity;
                    else
                        activity_to_combine_omi_consump = [];
                    end
                    if  ~isempty(rescale_factors)
                        for r = 4:size(activity_to_combine_omi,2)
                            temp1 = permute(activity_to_combine_omi(:,r,:),[1,3,2]);
                            temp2 = temp1 * rescale_factors(ids(i),r-3);
                            activity_to_combine_omi(:,r,:) = temp2;

                            temp1 = permute(activity_to_combine_omi_consump(:,r,:),[1,3,2]);
                            temp2 = temp1 * rescale_factors(ids(i),r-3);
                            activity_to_combine_omi_consump(:,r,:) = temp2;
                        end
                    end
                    if ~isnan(sample_rate) && this_sample_rate ~= sample_rate
                        activity_to_combine_omi = common_functions.interp_ta(activity_to_combine_omi,this_sample_rate,sample_rate);
                        activity_to_combine_omi_consump = common_functions.interp_ta(activity_to_combine_omi_consump,this_sample_rate,sample_rate);
                    end
                    if isempty(rew_value_omi) && ~isempty(activity_to_combine_omi)
                        rew_value_omi = activity_to_combine_omi;
                    else
                        rew_value_omi = cat(3,rew_value_omi,activity_to_combine_omi);
                    end

                    if isempty(rew_value_omi_consump) && ~isempty(activity_to_combine_omi_consump)
                        rew_value_omi_consump = activity_to_combine_omi_consump;
                    else
                        rew_value_omi_consump = cat(3,rew_value_omi_consump,activity_to_combine_omi_consump);
                    end
                end
            end
            ta_total.activity = rew_value;
            ta_total.mu = mean(rew_value,3,'omitnan');
            ta_total.sem = std(rew_value,[],3,'omitnan')/sqrt(size(rew_value,3)-1);
            if rew_consump_flag
                ta_total.activity_consump = rew_value_consump;
                ta_total.mu_consump = mean(rew_value_consump,3,'omitnan');
                ta_total.sem_consump = std(rew_value_consump,[],3,'omitnan')/sqrt(size(rew_value_consump,3)-1);
            end
            if omission_flag
                ta_total.activity_omi = rew_value_omi;
                ta_total.mu_omi = mean(rew_value_omi,3,'omitnan');
                ta_total.sem_omi = std(rew_value_omi,[],3,'omitnan')/sqrt(size(rew_value_omi,3)-1);
            end
            if omission_flag && rew_consump_flag
                ta_total.activity_omi_consump = rew_value_omi_consump;
                ta_total.mu_omi_consump = mean(rew_value_omi_consump,3,'omitnan');
                ta_total.sem_omi_consump = std(rew_value_omi_consump,[],3,'omitnan')/sqrt(size(rew_value_omi_consump,3)-1);
            end

        case "cue12"
            if mode == "cue1"
                mode_ci = 1;
            elseif mode == "cue2"
                mode_ci = 2;
            else
                error("mode_ci should only be 1 or 2.")
            end
            ta_value = [];
            for i = 1:length(ids)
                this_sample_rate = ta.("file"+string(i)).sample_rate;
                % process activity
                activity_to_combine = ta.("file"+string(i)).cueOn.("cue"+mode_ci).activity;
                if  ~isempty(rescale_factors)
                    for r = 4:size(activity_to_combine,2)
                        temp1 = permute(activity_to_combine(:,r,:),[1,3,2]);
                        temp2 = temp1 * rescale_factors(ids(i),r-3);
                        activity_to_combine(:,r,:) = temp2;
                    end
                end
                if ~isnan(sample_rate) && this_sample_rate ~= sample_rate
                    activity_to_combine = common_functions.interp_ta(activity_to_combine,this_sample_rate,sample_rate);
                    interpolation_flag = interpolation_flag + 1;
                end
                ta_value = cat(3,ta_value,activity_to_combine);
            end
            ta_total.activity = ta_value;
            ta_total.mu = mean(ta_value,3,'omitnan');
            ta_total.sem = std(ta_value,[],3,'omitnan')/sqrt(size(ta_value,3)-1);

        case "cue1"
            ta_value = [];
            for i = 1:length(ids)
                this_sample_rate = ta.("file"+string(i)).sample_rate;
                % process activity
                activity_to_combine = ta.("file"+string(i)).cueOn.cue1.activity;
                if  ~isempty(rescale_factors)
                    for r = 4:size(activity_to_combine,2)
                        temp1 = permute(activity_to_combine(:,r,:),[1,3,2]);
                        temp2 = temp1 * rescale_factors(ids(i),r-3);
                        activity_to_combine(:,r,:) = temp2;
                    end
                end
                if ~isnan(sample_rate) && this_sample_rate ~= sample_rate
                    activity_to_combine = common_functions.interp_ta(activity_to_combine,this_sample_rate,sample_rate);
                    interpolation_flag = interpolation_flag + 1;
                end
                ta_value = cat(3,ta_value,activity_to_combine);
            end
            ta_total.activity = ta_value;
            ta_total.mu = mean(ta_value,3,'omitnan');
            ta_total.sem = std(ta_value,[],3,'omitnan')/sqrt(size(ta_value,3)-1);

        case "cue2"
            ta_value = [];
            for i = 1:length(ids)
                this_sample_rate = ta.("file"+string(i)).sample_rate;
                % process activity
                activity_to_combine = ta.("file"+string(i)).cueOn.cue2.activity;
                if  ~isempty(rescale_factors)
                    for r = 4:size(activity_to_combine,2)
                        temp1 = permute(activity_to_combine(:,r,:),[1,3,2]);
                        temp2 = temp1 * rescale_factors(ids(i),r-3);
                        activity_to_combine(:,r,:) = temp2;
                    end
                end
                if ~isnan(sample_rate) && this_sample_rate ~= sample_rate
                    activity_to_combine = common_functions.interp_ta(activity_to_combine,this_sample_rate,sample_rate);
                    interpolation_flag = interpolation_flag + 1;
                end
                ta_value = cat(3,ta_value,activity_to_combine);
            end
            ta_total.activity = ta_value;
            ta_total.mu = mean(ta_value,3,'omitnan');
            ta_total.sem = std(ta_value,[],3,'omitnan')/sqrt(size(ta_value,3)-1);
        otherwise
            error("TBI")
    end
    % if interpolation_flag~= 0
    %     cprintf('magenta',"Different sample rates detected, spline interpolation applied to %d days.\n",interpolation_flag)
    % end
end

function [data_deliv,data_consump] = separate_delivery_consumption(data_in)
% for data_in containing both delivery and consumption, seperate them into two copies
    data_deliv = data_in;
    data_consump = data_in;

    f_names = string(fields(data_in)');
    % f_names_to_modify = f_names(contains(f_names,["unpred","rew"]));
    f_names_to_modify = f_names(contains(f_names,["rew","unpred"]));
    for f_name = f_names_to_modify
        string(fields(data_deliv.(f_name))')
        data_deliv.(f_name) = structfun(@remove_consump,data_deliv.(f_name),UniformOutput=false);
        data_consump.(f_name) = structfun(@modify_consump,data_consump.(f_name),UniformOutput=false);
    end
    
    function ta_mouse_out = remove_consump(ta_mouse_in)
        fnames = string(fields(ta_mouse_in)');
        tmp = fnames(contains(fnames,'consump'));
        ta_mouse_out = rmfield(ta_mouse_in,tmp);
    end

    function ta_mouse_out = modify_consump(ta_mouse_in)
        ta_mouse_out = struct;
        fnames = fields(ta_mouse_in)';
        tmp = cellfun(@(x) x(1:(find(x=='_')-1)),fnames,UniformOutput=false);
        old_fnames = string(fnames(~cellfun(@isempty,tmp)));
        new_fnames = string(tmp(~cellfun(@isempty,tmp)));
        for i=1:length(old_fnames)
            ta_mouse_out.(new_fnames(i)) = ta_mouse_in.(old_fnames(i));
        end
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

function output = get_lick_index(behav_cell, varargin)
    lick_window = [0,3];
    second_exclude_since_reward_onset = [0,3];
    lick_index_type = 2;
    lick_onset_gap = [];
    mouse_name = '';
    
    ip = inputParser;
    ip.addParameter('lick_window',[0,3])
    ip.addParameter('second_exclude_since_reward_onset',[0,3])
    ip.addParameter('lick_index_type',2)
    ip.addParameter('lick_onset_gap',[])
    ip.addParameter('mouse_name','')
    ip.parse(varargin{:})
    for j = fields(ip.Results)'
        eval([j{1},'=ip.Results.',j{1},';'])
    end
    lick_index_fhandle = get_lick_index_fhandle(lick_index_type);
    lick_window_1s = [2,3]; % the initial lick window i.e. 1s before cueon and 1s before rewon
    lick_window_text = (lick_window(1)-3)+" to "+(lick_window(2)-3)+"s";
    lick_window_tag = (lick_window(1)-3)+"_"+(lick_window(2)-3)+"s";
    lick_window_1s_text = (lick_window_1s(1)-3)+" to "+(lick_window_1s(2)-3)+"s";
    lick_window_1s_tag = (lick_window_1s(1)-3)+"_"+(lick_window_1s(2)-3)+"s";
    n_days = length(behav_cell);
    ids = 1:n_days;
    % preallocate
    cue1_pre_across = {}; % merge across days
    cue1_iti_across = {}; % merge across days
    cue2_pre_across = {}; % merge across days
    cue2_iti_across = {}; % merge across days
    cue1_index_across = {}; % merge across days
    cue2_index_across = {}; % merge across days

    single_trial_struct = struct;
    single_trial_taking_all_ITI_frame = struct;
    single_trial_rew = struct;

    % get data
    for i = ids
        this_cell = behav_cell{i};
        if isstring(this_cell)
            behav7 = load(this_cell);
        elseif isstruct(this_cell)
            behav7 = this_cell;
        else
            error("IfElseError.")
        end
        % werid things to do
        fnames = string(fields(behav7)');
        for fname = fnames
            if isnumeric(behav7.(fname)) && size(behav7.(fname),1) >= (size(behav7.timestamp,1)-10)
                if any(isnan(behav7.(fname)))
                    warning("Nans are removed, might cause potentially error.")
                end
                behav7.(fname) = behav7.(fname)(~isnan(behav7.(fname)));
            end
        end
        % get trial info
        frequency = common_functions.get_sr(behav7);
        trials = behav7.experimentSetup.exp.trials;
        n_trials = trials.trialNum(end);
        save_suffix = '';
        % use lick bout onset or not
        if isempty(lick_onset_gap)
            % skip
        else
            save_suffix = [save_suffix,'_bout_onset'];
            gap_threshold = round(frequency * lick_onset_gap);
            lick_onsets = get_lick_onset(find(diff(behav7.lick)==1)+1,gap_threshold);
            lickon_frames = zeros(size(behav7.lick));
            lickon_frames(lick_onsets) = 1;
            behav7.lick_count = lickon_frames;
        end

        % find corresponding cue frame for each trial
        cues = behav7.stimulus_sound2 | behav7.stimulus_led;
        cues_frames = find(diff(cues)==1);

        if length(cues_frames) < n_trials
            warning("Actual cue count doesn't match with trials assigned in behav. Modify trials by removing last trial(s)."...
                + "or maybe the trials have two columns? IDK.")
            trials = trials(1:length(cues_frames),:);
        end

        rewarded_trialnums = trials{trials.rew==1,"trialNum"}';
        cue1_rew = trials(trials.cueRL==1&trials.rew==1,:);
        cue2_rew = trials(trials.cueRL==2&trials.rew==1,:);
        cue1_omi = trials(trials.cueRL==1&trials.rew==0,:);
        cue2_omi = trials(trials.cueRL==2&trials.rew==0,:);
        rewards = behav7.reward & cues;
        rewards_frames = find(diff(rewards)==1);
    
        % get licks
        cue1_rew_frames_cue = [cues_frames(cue1_rew.trialNum)-3*frequency,cues_frames(cue1_rew.trialNum)+3*frequency-1];
        [~,cue1_rew_frames_rew_index] = intersect(rewarded_trialnums,cue1_rew.trialNum,'stable');
        cue1_rew_frames_rew = [rewards_frames(cue1_rew_frames_rew_index),rewards_frames(cue1_rew_frames_rew_index)+3*frequency-1];

        cue2_rew_frames_cue = [cues_frames(cue2_rew.trialNum)-3*frequency,cues_frames(cue2_rew.trialNum)+3*frequency-1];
        [~,cue2_rew_frames_rew_index] = intersect(rewarded_trialnums,cue2_rew.trialNum,'stable');
        cue2_rew_frames_rew = [rewards_frames(cue2_rew_frames_rew_index),rewards_frames(cue2_rew_frames_rew_index)+3*frequency-1];
    
        cue1_rew_lick = nan(size(cue1_rew_frames_cue,1),frequency*9);
        cue2_rew_lick = nan(size(cue2_rew_frames_cue,1),frequency*9);

        for i1 = 1:size(cue1_rew_frames_cue,1)
            if cue1_rew_frames_cue(i1,1) <= 0 % if trial not complete remove it
                continue
            end
            cue1_rew_lick(i1,1:6*frequency) = behav7.lick_count(cue1_rew_frames_cue(i1,1):cue1_rew_frames_cue(i1,2));
            cue1_rew_lick(i1,6*frequency+1:9*frequency) = behav7.lick_count(cue1_rew_frames_rew(i1,1):cue1_rew_frames_rew(i1,2));
        end
        for i1 = 1:size(cue2_rew_frames_cue,1)
            if cue2_rew_frames_cue(i1,1) <= 0
                continue
            end
            cue2_rew_lick(i1,1:6*frequency) = behav7.lick_count(cue2_rew_frames_cue(i1,1):cue2_rew_frames_cue(i1,2));
            cue2_rew_lick(i1,6*frequency+1:9*frequency) = behav7.lick_count(cue2_rew_frames_rew(i1,1):cue2_rew_frames_rew(i1,2));
        end
    
        cue1_omi_lick = nan(size(cue1_omi.trialNum,1),frequency*9);
        cue2_omi_lick = nan(size(cue2_omi.trialNum,1),frequency*9);
        % if omission
        if ~isempty(cue1_omi)
            cue1_omi_frames_cue = [cues_frames(cue1_omi.trialNum)-3*frequency,cues_frames(cue1_omi.trialNum)+3*frequency-1];
            cue1_omi_frames_rew = [cues_frames(cue1_omi.trialNum)+3*frequency, cues_frames(cue1_omi.trialNum) + 6*frequency-1];
    
            cue2_omi_frames_cue = [cues_frames(cue2_omi.trialNum)-3*frequency,cues_frames(cue2_omi.trialNum)+3*frequency-1];
            cue2_omi_frames_rew = [cues_frames(cue2_omi.trialNum) + 3*frequency,cues_frames(cue2_omi.trialNum) + 6*frequency-1];
    
            for i1 = 1:size(cue1_omi_frames_cue,1)
                if cue1_omi_frames_cue(i1,1) <= 0
                    continue
                end
                cue1_omi_lick(i1,1:6*frequency) = behav7.lick_count(cue1_omi_frames_cue(i1,1):cue1_omi_frames_cue(i1,2));
                cue1_omi_lick(i1,6*frequency+1:9*frequency) = behav7.lick_count(cue1_omi_frames_rew(i1,1):cue1_omi_frames_rew(i1,2));
            end
            for i1 = 1:size(cue2_omi_frames_cue,1)
                if cue2_omi_frames_cue(i1,1) <= 0
                    continue
                end
                cue2_omi_lick(i1,1:6*frequency) = behav7.lick_count(cue2_omi_frames_cue(i1,1):cue2_omi_frames_cue(i1,2));
                cue2_omi_lick(i1,6*frequency+1:9*frequency) = behav7.lick_count(cue2_omi_frames_rew(i1,1):cue2_omi_frames_rew(i1,2));
            end
        end
    
        cue1_lick = cue1_rew_lick;
        cue2_lick = cue2_rew_lick;
        if ~isempty(cue1_omi)
            cue1_lick = [cue1_lick;cue1_omi_lick];
            cue2_lick = [cue2_lick;cue2_omi_lick];
        end
        
        cue1_iti_lick = nan([size(cue1_lick,1),2]);
        cue1_pre_lick = nan([size(cue1_lick,1),2]);
        cue1_lick_index = nan([size(cue1_lick,1),2]);
        cue2_iti_lick = nan([size(cue2_lick,1),2]);
        cue2_pre_lick = nan([size(cue2_lick,1),2]);
        cue2_lick_index = nan([size(cue2_lick,1),2]);
        % get total lick of single trials during given lick_window or the 1s before (the initial way)
        try
        cue1_iti_lick(:,1) = sum(cue1_lick(:,lick_window(1)*frequency+1:lick_window(2)*frequency),2,"omitnan");
        cue1_pre_lick(:,1) = sum(cue1_lick(:,(lick_window(1)+3)*frequency+1:(lick_window(2)+3)*frequency),2,"omitnan");
        cue2_iti_lick(:,1) = sum(cue2_lick(:,lick_window(1)*frequency+1:lick_window(2)*frequency),2,"omitnan");
        cue2_pre_lick(:,1) = sum(cue2_lick(:,(lick_window(1)+3)*frequency+1:(lick_window(2)+3)*frequency),2,"omitnan");
        catch err
            keyboard
        end
        %       the initial way took only 1s, meaning lick_window = [2,3]
        cue1_iti_lick(:,2) = sum(cue1_lick(:,lick_window_1s(1)*frequency+1:lick_window_1s(2)*frequency),2,"omitnan");
        cue1_pre_lick(:,2) = sum(cue1_lick(:,(lick_window_1s(1)+3)*frequency+1:(lick_window_1s(2)+3)*frequency),2,"omitnan");
        cue2_iti_lick(:,2) = sum(cue2_lick(:,lick_window_1s(1)*frequency+1:lick_window_1s(2)*frequency),2,"omitnan");
        cue2_pre_lick(:,2) = sum(cue2_lick(:,(lick_window_1s(1)+3)*frequency+1:(lick_window_1s(2)+3)*frequency),2,"omitnan");
        % compute each lick index for the two licks above
        for tmp_i = 1:2
            cue1_lick_index(:,tmp_i) = lick_index_fhandle(cue1_pre_lick(:,tmp_i),cue1_iti_lick(:,tmp_i));
            cue2_lick_index(:,tmp_i) = lick_index_fhandle(cue2_pre_lick(:,tmp_i),cue2_iti_lick(:,tmp_i));
        end

        cue1_pre_across = cat(2,cue1_pre_across,{cue1_pre_lick});
        cue1_iti_across = cat(2,cue1_iti_across,{cue1_iti_lick});
        cue2_pre_across = cat(2,cue2_pre_across,{cue2_pre_lick});
        cue2_iti_across = cat(2,cue2_iti_across,{cue2_iti_lick});
        cue1_index_across = cat(2,cue1_index_across,{cue1_lick_index});
        cue2_index_across = cat(2,cue2_index_across,{cue2_lick_index});

        single_trial_struct.("day"+string(i)).cue1_iti_lick = array2table(cue1_iti_lick,VariableNames=[lick_window_text,lick_window_1s_text]);
        single_trial_struct.("day"+string(i)).cue1_pre_lick = array2table(cue1_pre_lick,VariableNames=[lick_window_text,lick_window_1s_text]);
        single_trial_struct.("day"+string(i)).cue2_iti_lick = array2table(cue2_iti_lick,VariableNames=[lick_window_text,lick_window_1s_text]);
        single_trial_struct.("day"+string(i)).cue2_pre_lick = array2table(cue2_pre_lick,VariableNames=[lick_window_text,lick_window_1s_text]);
        single_trial_struct.("day"+string(i)).cue1_lick_index = array2table(cue1_lick_index,VariableNames=[lick_window_text,lick_window_1s_text]);
        single_trial_struct.("day"+string(i)).cue2_lick_index = array2table(cue2_lick_index,VariableNames=[lick_window_text,lick_window_1s_text]);

        % now get ITI licking which is not restricted within 3 seconds
        ITI_bit = get_ITI_pav(behav7,remove_reward=second_exclude_since_reward_onset,remove_lick=[],remove_cue=1);
        q1_frames = [cues_frames(cue1_rew.trialNum);cues_frames(cue1_omi.trialNum)];
        q2_frames = [cues_frames(cue2_rew.trialNum);cues_frames(cue2_omi.trialNum)];
        unpred_frames = find(diff(behav7.reward & ~behav7.stimulus_led & ~behav7.stimulus_sound2)==1);

        q1_frames_leading = arrayfun(@(x) find_previous_zero(ITI_bit,x),q1_frames);
        q2_frames_leading = arrayfun(@(x) find_previous_zero(ITI_bit,x),q2_frames);
        unpred_frames_leading = arrayfun(@(x) find_previous_zero(ITI_bit,x),unpred_frames);
        
        % cue
        q1_lick = cell2mat(arrayfun(@(x1,x2) lick_count_rate(behav7.lick_count,x1,x2,frequency),q1_frames_leading,q1_frames,UniformOutput=false));
        q2_lick = cell2mat(arrayfun(@(x1,x2) lick_count_rate(behav7.lick_count,x1,x2,frequency),q2_frames_leading,q2_frames,UniformOutput=false));
        unpred_lick = cell2mat(arrayfun(@(x1,x2) lick_count_rate(behav7.lick_count,x1,x2,frequency),unpred_frames_leading,unpred_frames,UniformOutput=false));
        
        q1_out = [q1_frames_leading,q1_frames,q1_lick,cue1_pre_lick(:,1),cue1_pre_lick(:,1)/(lick_window(2)-lick_window(1))];
        q2_out = [q2_frames_leading,q2_frames,q2_lick,cue2_pre_lick(:,1),cue2_pre_lick(:,1)/(lick_window(2)-lick_window(1))];
        unpred_out = [unpred_frames_leading,unpred_frames,unpred_lick];
        q1_index = lick_index_fhandle(q1_out(:,6),q1_out(:,4));
        q2_index = lick_index_fhandle(q2_out(:,6),q2_out(:,4));
        single_trial_taking_all_ITI_frame.("day"+string(i)).cue1_lick = array2table([q1_out,q1_index],VariableNames=["iti_start","iti_end","iti_lick_count","iti_lick_rate","pre_lick_count","pre_lick_rate","lick_index"]);
        single_trial_taking_all_ITI_frame.("day"+string(i)).cue2_lick = array2table([q2_out,q2_index],VariableNames=["iti_start","iti_end","iti_lick_count","iti_lick_rate","pre_lick_count","pre_lick_rate","lick_index"]);
        if ~isempty(unpred_out)
            single_trial_taking_all_ITI_frame.("day"+string(i)).unpred_lick = array2table(unpred_out,VariableNames=["iti_start","iti_end","iti_lick_count","iti_lick_rate"]);
        else
            single_trial_taking_all_ITI_frame.("day"+string(i)).unpred_lick = table;
        end

        % rew
        r1_lick = cell2mat(arrayfun(@(x1,x2) lick_count_rate(behav7.lick_count,x1,x2,frequency),cue1_rew_frames_rew(:,1)+1,cue1_rew_frames_rew(:,2)+1,UniformOutput=false));
        r2_lick = cell2mat(arrayfun(@(x1,x2) lick_count_rate(behav7.lick_count,x1,x2,frequency),cue2_rew_frames_rew(:,1)+1,cue2_rew_frames_rew(:,2)+1,UniformOutput=false));

        if ~isempty(cue1_omi)
            omi1_lick = cell2mat(arrayfun(@(x1,x2) lick_count_rate(behav7.lick_count,x1,x2,frequency),cue1_omi_frames_rew(:,1)+1,cue1_omi_frames_rew(:,2)+1,UniformOutput=false));
            omi2_lick = cell2mat(arrayfun(@(x1,x2) lick_count_rate(behav7.lick_count,x1,x2,frequency),cue2_omi_frames_rew(:,1)+1,cue2_omi_frames_rew(:,2)+1,UniformOutput=false));
        end

        rew_unpred_lick = cell2mat(arrayfun(@(x1,x2) lick_count_rate(behav7.lick_count,x1,x2,frequency),unpred_frames+1,unpred_frames+3*frequency+1,UniformOutput=false));
        
        single_trial_rew.("day"+string(i)).rew1_lick = array2table(r1_lick,VariableNames=["rew_lick_count","rew_lick_rate"]);
        single_trial_rew.("day"+string(i)).rew2_lick = array2table(r2_lick,VariableNames=["rew_lick_count","rew_lick_rate"]);

        if ~isempty(cue1_omi)
            single_trial_rew.("day"+string(i)).omi1_lick = array2table(omi1_lick,VariableNames=["omi_lick_count","omi_lick_rate"]);
            single_trial_rew.("day"+string(i)).omi2_lick = array2table(omi2_lick,VariableNames=["omi_lick_count","omi_lick_rate"]);
        else
            single_trial_rew.("day"+string(i)).omi1_lick = table;
            single_trial_rew.("day"+string(i)).omi2_lick = table;
        end

        if ~isempty(rew_unpred_lick)
            single_trial_rew.("day"+string(i)).unpred_lick = array2table(rew_unpred_lick,VariableNames=["rew_lick_count","rew_lick_rate"]);
        else
            single_trial_rew.("day"+string(i)).unpred_lick = table;
        end
    end

    % classical licking across day
    tmp_i_tag = ["input_lick_window","1s"];
    sts_a = struct;
    [sts_a.cue1_iti_lick,sts_a.cue1_pre_lick,sts_a.cue2_iti_lick,sts_a.cue2_pre_lick,sts_a.cue1_index_across,sts_a.cue2_index_across] = ...
        deal(cue1_iti_across,cue1_pre_across,cue2_iti_across,cue2_pre_across,cue1_index_across,cue2_index_across);
    for tmp_i = 1:2
        single_trial_struct.("across_"+tmp_i_tag(tmp_i)) = ...
            structfun(@(x1) specific_merge_across_day_1(cellfun(@(x) x(:,tmp_i),x1,UniformOutput=false)),sts_a,UniformOutput=false);
    end

    % whole ITI licking across day
    sts_a = struct;
    sts_a.cue1 = cell(2,n_days); sts_a.cue2 = cell(2,n_days); sts_a.unpred = cell(1,n_days);
    for i = 1:n_days
        sts_a.cue1{1,i} = single_trial_taking_all_ITI_frame.("day"+string(i)).cue1_lick{:,"iti_lick_rate"};
        sts_a.cue1{2,i} = single_trial_taking_all_ITI_frame.("day"+string(i)).cue1_lick{:,"pre_lick_rate"};
        sts_a.cue2{1,i} = single_trial_taking_all_ITI_frame.("day"+string(i)).cue2_lick{:,"iti_lick_rate"};
        sts_a.cue2{2,i} = single_trial_taking_all_ITI_frame.("day"+string(i)).cue2_lick{:,"pre_lick_rate"};
        if ~isempty(single_trial_taking_all_ITI_frame.("day"+string(i)).unpred_lick)
            sts_a.unpred{1,i} = single_trial_taking_all_ITI_frame.("day"+string(i)).unpred_lick{:,"iti_lick_rate"};
        else
            sts_a.unpred{1,i} = [];
        end
    end
    cue1cue2unpred_ntrial = [max(cellfun(@(x) size(x,1),sts_a.cue1),[],'all'),max(cellfun(@(x) size(x,1),sts_a.cue2),[],'all'),max(cellfun(@(x) size(x,1),sts_a.unpred),[],'all')];
    
    sts_a.cue1 = cellfun(@(x) addnan(x,cue1cue2unpred_ntrial(1)), sts_a.cue1,UniformOutput=0);
    sts_a.cue2 = cellfun(@(x) addnan(x,cue1cue2unpred_ntrial(2)), sts_a.cue2,UniformOutput=0);
    sts_a.unpred = cellfun(@(x) addnan(x,cue1cue2unpred_ntrial(3)), sts_a.unpred,UniformOutput=0);

    sts_across.cue1_iti = cell2mat(sts_a.cue1(1,:));
    sts_across.cue1_pre = cell2mat(sts_a.cue1(2,:));
    sts_across.cue2_iti = cell2mat(sts_a.cue2(1,:));
    sts_across.cue2_pre = cell2mat(sts_a.cue2(2,:));
    sts_across.unpred_iti = cell2mat(sts_a.unpred(1,:));

    sts_across.cue1_index = lick_index_fhandle(sts_across.cue1_pre,sts_across.cue1_iti);
    sts_across.cue2_index = lick_index_fhandle(sts_across.cue2_pre,sts_across.cue2_iti);
    single_trial_taking_all_ITI_frame.across = sts_across;
    clear sts_across;

    % rew licking across day
    sts_a = struct;
    sts_a.rew1 = cell(1,n_days); sts_a.rew2 = cell(1,n_days);
    sts_a.omi1 = cell(1,n_days); sts_a.omi2 = cell(1,n_days);
    sts_a.unpred = cell(1,n_days);
    for i = 1:n_days
        sts_a.rew1{1,i} = single_trial_rew.("day"+string(i)).rew1_lick{:,"rew_lick_rate"};
        sts_a.rew2{1,i} = single_trial_rew.("day"+string(i)).rew2_lick{:,"rew_lick_rate"};
        if ~isempty(single_trial_rew.("day"+string(i)).omi1_lick)
            sts_a.omi1{1,i} = single_trial_rew.("day"+string(i)).omi1_lick{:,"omi_lick_rate"};
            sts_a.omi2{1,i} = single_trial_rew.("day"+string(i)).omi2_lick{:,"omi_lick_rate"};
        else
            sts_a.omi1{1,i} = [];
            sts_a.omi2{1,i} = [];
        end
        if ~isempty(single_trial_rew.("day"+string(i)).unpred_lick)
            sts_a.unpred{1,i} = single_trial_rew.("day"+string(i)).unpred_lick{:,"rew_lick_rate"};
        else
            sts_a.unpred{1,i} = [];
        end
    end
    rewomiunpre_ntrial = [max(cellfun(@(x) size(x,1),sts_a.rew1),[],'all'),max(cellfun(@(x) size(x,1),sts_a.rew2),[],'all'),...
        max(cellfun(@(x) size(x,1),sts_a.omi1),[],'all'),max(cellfun(@(x) size(x,1),sts_a.omi2),[],'all'),...
        max(cellfun(@(x) size(x,1),sts_a.unpred),[],'all')];

    sts_a.rew1 = cellfun(@(x) addnan(x,rewomiunpre_ntrial(1)), sts_a.rew1,UniformOutput=0);
    sts_a.rew2 = cellfun(@(x) addnan(x,rewomiunpre_ntrial(2)), sts_a.rew2,UniformOutput=0);
    sts_a.omi1 = cellfun(@(x) addnan(x,rewomiunpre_ntrial(3)), sts_a.omi1,UniformOutput=0);
    sts_a.omi2 = cellfun(@(x) addnan(x,rewomiunpre_ntrial(4)), sts_a.omi2,UniformOutput=0);
    sts_a.unpred = cellfun(@(x) addnan(x,rewomiunpre_ntrial(5)), sts_a.unpred,UniformOutput=0);
    
    sts_across.rew1 = cell2mat(sts_a.rew1(1,:));
    sts_across.rew2 = cell2mat(sts_a.rew2(1,:));
    sts_across.omi1 = cell2mat(sts_a.omi1(1,:));
    sts_across.omi2 = cell2mat(sts_a.omi2(1,:));
    sts_across.unpred = cell2mat(sts_a.unpred(1,:));

    single_trial_rew.across = sts_across;
    clear sts_across;

    output = struct;
    output.single_trial_struct = single_trial_struct;
    output.single_trial_taking_all_ITI_frame = single_trial_taking_all_ITI_frame;
    output.single_trial_rew = single_trial_rew;
    
    function arr_out = specific_merge_across_day_1(arr_cell_in)
        max_n_trial = max(cellfun(@(x) size(x,1),arr_cell_in),[],"omitmissing");
        arr_out = cellfun(@(x) addnan(x,max_n_trial),arr_cell_in,UniformOutput=false);
        arr_out = cell2mat(arr_out);
    end

    function out = addnan(arr,sz)
        out = cat(1,arr,nan(sz-size(arr,1),1));
    end
    
    function lick_index_fhandle = get_lick_index_fhandle(lick_index_type)
        % return lick index algorithm based on lick index type
        % 1:= (pre-iti)/(pre+iti); the standard index 
        % 2:= log(pre) * (pre-iti)/(pre+iti) % should probably do log(pre+1) but anyway...
        % 3:= (pre-"iti")/(pre+"iti"); the standard index where the "iti" is the mean iti across each single days
        % 4:= log(pre) * (pre-"iti")/(pre+"iti"); similar to 3 but the "iti" is the mean iti across each single days
        fhandles = {@fhandle1,@fhandle2,@fhandle3,@fhandle4};
        lick_index_fhandle = fhandles{lick_index_type};

        function index = fhandle1(pre,iti)
            index = (pre-iti)./(pre+iti);
        end
        function index = fhandle2(pre,iti)
            index = log(pre).*(pre-iti)./(pre+iti);
            index(isinf(log(pre))) = 0;
        end
        function index = fhandle3(pre,iti)
            mean_iti = mean(iti,"omitmissing");
            index = (pre-mean_iti)./(pre+mean_iti);
        end
        function index = fhandle4(pre,iti)
            mean_iti = mean(iti,"omitmissing");
            index = log(pre).*(pre-mean_iti)./(pre+mean_iti);
            index(isinf(log(pre))) = 0;
        end
    end

    function out = find_previous_zero(bits,frame)
        out = find(bits(1:frame)==0,1,"last")+1;
        if isempty(out)
            out = 1;
        end
    end
    
    function out = lick_count_rate(lick_count,leading,ending,sr)
        out = nan(1,2);
        out(1) = sum(lick_count(leading:ending));
        out(2) = mean(lick_count(leading:ending))*sr;
    end

end

function data_out = build_plotdata_by_session(data_in,session_info,null_window_s,varargin)
    % combine data into sessions based on session_info
    flag_rebaseline = true;

    ip = inputParser;
    ip.addParameter("flag_rebaseline",true)
    ip.parse(varargin{:})
    for j = fields(ip.Results)'
        eval([j{1}, '= ip.Results.', j{1}, ';']);
    end

    mouse_names = string(fields(data_in)');
    
    mouse_most_fr = dictionary;
    mouse_most_fr(mouse_names) = 18;
    mouse_most_fr(["G12","G15"]) = [30,30];

    % manual
    session_names = string(fields(session_info)');

    for m_name = mouse_names
        if contains(m_name,"G20")
            continue
        end
        mouse_sr = mouse_most_fr(m_name);
        if isempty(session_info)
            data_out.(m_name) = cell2table(data_in.(m_name).data(:,[1,3]),VariableNames=["sample_rate","data_1"]);
            data_out.(m_name).Properties.RowNames = "session_"+string(1:size(data_out.(m_name),1));
            continue
        end
        this_sessions = session_names(contains(session_names,m_name));
        n_session = length(this_sessions);
        table_cell = cell(n_session,3);
        for i = 1:n_session
            session_ids = session_info.(this_sessions(i));
            to_merge = data_in.(m_name).data(session_ids,[1,3]);
            to_interp = [to_merge{:,1}]~=mouse_sr;
            if any(to_interp)
                to_merge(to_interp,2) = cellfun(@(x) interp1(1:size(x,1),x,linspace(1,size(x,1),2*mouse_sr)),... % this 2*mouse_sr is sloppy, cuz I assume the signal window is always 2 second
                    to_merge(to_interp,2),UniformOutput=0);
            end
            merged = [];
            for j = 1:size(to_merge,1)
                merged = cat(3,merged,to_merge{j,2});
            end
            % get null
            null_window = 1:round(null_window_s*mouse_sr);
            null_populations = merged(null_window,:,:);
            null_std = std(mean(null_populations,3,"omitmissing"),[],1,"omitmissing");
            % rebaseline if flag is on
            if flag_rebaseline
                size_ori = size(merged);
                null_mu = mean(null_populations(:,4:end,:),[1,3],"omitmissing");
                merged(:,4:end,:) = merged(:,4:end,:) - repmat(null_mu,[size_ori(1),1,size_ori(3)]);
            end
            table_cell{i,1} = mouse_sr;
            table_cell{i,2} = merged;
            table_cell{i,3} = {null_std};
        end
        data_out.(m_name) = cell2table(table_cell,VariableNames=["sample_rate","data_1","data_1_null"],RowNames=this_sessions);
    end
end

function cwa_struct = cwa_table_to_struct(cwa_table)
    % convert cwa table to classical cwa struct
    cwa_struct = struct;

    mnames = string(fields(cwa_table)');
    for mname = mnames
        snames = arrayfun(@(x) string(x{1}(4:end)),cwa_table.(mname).Properties.RowNames)';
        for i = 1:length(snames)
            cwa_struct.(snames(i)).(mname).sr = cwa_table.(mname){i,1};
            cwa_struct.(snames(i)).(mname).activity = cwa_table.(mname){i,2}{1};
            cwa_struct.(snames(i)).(mname).mu = mean(cwa_table.(mname){i,2}{1},3,"omitmissing");
            cwa_struct.(snames(i)).(mname).sem = std(cwa_table.(mname){i,2}{1},[],3,"omitmissing")/sqrt(size(cwa_table.(mname){i,2}{1},3));
            cwa_struct.(snames(i)).(mname).null = cwa_table.(mname){i,3}{1};
        end
    end
end

function tas_out = specialized_add_to_tri_avg_place_holder()
    m_list1 = ["G12","G15","G17","G19","G21","G22","G23","G24","G25","G26","G27"];
    % m_list2 = ["G12","G15","G17","G19","G21","G22","G23","G24"];

    tas_out = struct;    
    for m = m_list1
        tas_out.late_early_diff.(m) = "placeholder";
        tas_out.LED_omission_late_diff.(m) = "placeholder";
        if ~contains(m,["G25","G26","G27"])
            tas_out.Tone_omission_late_diff.(m) = "placeholder";
        end
    end
end

function [pk_pks,pk_locs,pk_pks_single,pk_significance] = get_peak(roi_traces,roi_traces_mu,pk_window,min_amp_threshold,varargin)
    exclude_peaks = [];

    ip = inputParser;
    ip.addParameter('exclude_peaks',[]); % if this is not empty, it means to find rebounds
    ip.parse(varargin{:})
    for j = fields(ip.Results)'
        eval([j{1},'=ip.Results.',j{1},';']);
    end

    n_trial = size(roi_traces,3);
    n_rois = size(roi_traces,2)-3;
    
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

function [output,data_path] = ROI_loader_DA(path,day)
    % load data by mouse_name and day_id
    dir_info = dir(path);
    file_names = string({dir_info.name});
    record_dates = arrayfun(@get_date ,file_names,UniformOutput=false);
    tmp = ~cellfun(@isempty,record_dates);
    file_names = file_names(tmp);
    record_dates = string(record_dates(tmp)); clear("tmp");
    % assuming files are ordered in ascending alphabetical order (i.e. NTFS style)
    % if not, sort based on record_dates (to be implemented if needed)
    data_path = [char(path),'\',char(file_names(day))];
    output = load(data_path);
    function out = get_date(x)
        out = char(x);
        tmp = regexp(out,'\d{6}','once');
        out = out(tmp:tmp+5);
    end
end

function roidata = frequency_pass(roidata,Fourier_frequency,sample_rate)
    if any(isnan(roidata.Fc),"all") % pad nan is nan exist
        warning("ROI.Fc contains NaN, padded with mean.")
        n_rois  = size(roidata.Fc,2);
        nan_mask = isnan(roidata.Fc);
        for i=1:n_rois
            roidata.Fc(nan_mask(:,i),i) = mean(roidata.Fc(:,i),"omitmissing");
        end
    end
    Fourier_flag = any(Fourier_frequency);
    if Fourier_flag
        padsize7 = [round(size(roidata.Fc,1)/2),0]; % pad data before any filtering with half_length symmetric
        Fc_padded7 = padarray(roidata.Fc,padsize7,"symmetric","both");
        if ~isnan(Fourier_frequency(1)) && ~isnan(Fourier_frequency(2)) % bandpass
            Fc_filtered7 = bandpass(Fc_padded7,Fourier_frequency,sample_rate);
        elseif ~isnan(Fourier_frequency(1)) && isnan(Fourier_frequency(2)) % highpass
            Fc_filtered7 = highpass(Fc_padded7,Fourier_frequency(1),sample_rate);
        elseif isnan(Fourier_frequency(1)) && ~isnan(Fourier_frequency(2)) % lowpass
            Fc_filtered7 = lowpass(Fc_padded7,Fourier_frequency(2),sample_rate);
        else
            Fc_filtered7 = Fc_padded7;
        end
        roidata.Fc = Fc_filtered7(padsize7(1)+1:end-padsize7(1),:);
    end
    if any(isnan(roidata.Fc),"all") % reassign nan back
        roidata.Fc(nan_mask) = nan;
    end
end

function output = extract_task_info(task_info)
    % get event time from task info
    for ci = 1:2
        this_cue_bit = logical(task_info.cue==ci);
        output.("cue"+ci) = task_info.cueStart(this_cue_bit);
        output.("rew"+ci) = task_info.rewOn(this_cue_bit);
        output.("rewconsump"+ci) = task_info.rewLick(this_cue_bit);
    end
    output.unpred = task_info.unpredRewOn;
    output.unpredconsump = task_info.unpredRewLick;
end

function new = merge_pav2cue_data(old,new,main_sr)
    % modified from function with same name in task_TA_Glu.m
    % merge old pav2cue data struct into new ones
    for c_i = 1:2
        c = "cue"+c_i;

        % get cue
        try % see if activity exist
            size(new.cueOn.(c).activity);
        catch
            new.cueOn.(c).activity = [];
            new.cueOn.(c).session_num_trials = [];
        end
        n_trial = size(old.(c).activity,3);
        old_act = old.(c).activity;
        old.sample_rate = round(length(old.cue1.windowIdx)/4); % sloppy
        if old.sample_rate ~= main_sr
            to_cat = random_functions.interp_ta(old_act,old.sample_rate,main_sr);
        else
            to_cat = old_act;
        end
        new.cueOn.(c).activity = cat(3,new.cueOn.(c).activity,to_cat);
        new.cueOn.(c).session_num_trials = cat(2,new.cueOn.(c).session_num_trials,n_trial);
        
        % % get rew (ok reward is removed from the paper, ignore this)
        % try % see if activity exist
        %     size(new.rewOn.(c).rew.activity);
        % catch
        %     new.rewOn.(c).rew.activity = [];
        %     new.rewOn.(c).rew.session_num_trials = [];
        % end
        % n_trial = size(old.cueOn.(c).rew.activity,3);
        % new.rewOn.(c).rew.activity = cat(3,new.rewOn.(c).rew.activity,old.rewOn.(c).rew.activity);
        % new.rewOn.(c).rew.session_num_trials = cat(2,new.rewOn.(c).rew.session_num_trials,n_trial);
    end
end

