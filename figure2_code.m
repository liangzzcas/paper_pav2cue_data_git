%% figure 2B
% top
close all;clear;clc;
cf = [pwd,'\'];
mouse_name = "G23";
this_ta = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']).cwa_raw.(mouse_name);
session_info = common_functions.get_training_info();
this_pavday = session_info{[session_info{:,1}]==mouse_name,2};

imagesc_ta = nan(72,this_pavday);
for di = 1:this_pavday
    imagesc_ta(:,di) = this_ta.("file"+di).cueOn.cue1.mu(:,1);
end

fig = figure(Position=[100,100,800,600]);
ax=axes(fig);
imagesc(ax,(1:72)/18-1,1:size(imagesc_ta,2),imagesc_ta'/max(imagesc_ta,[],"all","omitmissing"));
colorbar(ax);clc

xline(ax,0,LineWidth=1.5,Color=[1,1,1]);
xlabel(ax,"Time (s)"); ylabel(ax, "Days")
saveas(fig,[cf,'fig2B_top.png'],'png')

% bottom
raw_path = [pwd,'\raw_data\'];
mouse_name = ["G23"];

mouse_path = [raw_path,char(mouse_name),'\'];
tmp = dir(mouse_path);
tmp = string({tmp.name});
dates = tmp(~contains(tmp,[".",".."])); clear("tmp");
load_ids = 1:this_pavday;
% the following steps take a long time since it involve a complete processing of raw data
output = get_triavg(mouse_path, mouse_name, dates, load_ids, load_ft_fr=[0.3,nan], load_ft_fr_1=[nan,nan],get_405=0,cue2Tone_shift_forward_flag=0);
imagesc_ta = nan(72,this_pavday);
for di = 1:this_pavday
    imagesc_ta(:,di) = output.("file"+di).cueOn.cue2.mu(:,1);
end

fig = figure(Position=[100,100,800,600]);
ax=axes(fig);
imagesc(ax,(1:72)/18-1,1:size(imagesc_ta,2),imagesc_ta'/max(imagesc_ta,[],"all","omitmissing"));
colorbar(ax);clc

xline(ax,0,LineWidth=1.5,Color=[1,1,1]);
xlabel(ax,"Time (s)"); ylabel(ax, "Days")
saveas(fig,[cf,'fig2B_bot.png'],'png')
delete(fig)

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 2C
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G22","G21","G23","G24"];
lick_data = load([cf,'processed_and_organized_data\across_mice_lick_index_data_whole_ITI.mat']);
training_info = common_functions.get_training_info();

fig = figure(position=[100,100,800,600]);
ax = axes(fig);
hold(ax,"on")
for mouse_name = mouse_names
    last_pav_day = training_info{[training_info{:,1}]==mouse_name,2};
    plot_data = lick_data.(mouse_name).single_trial_struct.across_1s.cue1_index_across(:,1:last_pav_day);
    datamu = mean(plot_data,1,"omitmissing");
    datasem = std(plot_data,[],1,"omitmissing")/sqrt(size(plot_data,1));
    errorbar(ax,1:last_pav_day,datamu,datasem);
end
hold(ax,"on")
saveas(fig,[cf,'fig2C.png'],'png')
delete(fig)

fig = figure(position=[100,100,800,600]);
ax = axes(fig);
hold(ax,"on")
for mouse_name = mouse_names
    last_pav_day = training_info{[training_info{:,1}]==mouse_name,2};
    plot_data = lick_data.(mouse_name).single_trial_struct.across_1s.cue2_index_across(:,1:last_pav_day);
    datamu = mean(plot_data,1,"omitmissing");
    datasem = std(plot_data,[],1,"omitmissing")/sqrt(size(plot_data,1));
    errorbar(ax,1:last_pav_day,datamu,datasem);
end
hold(ax,"on")
saveas(fig,[cf,'fig2C.png'],'png')
delete(fig)

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 2FH
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);

plot_info_1 = {"G12",39,30;"G21",27,18;"G19",32,18;};
plot_imagesc_clim = {0.05;0.06;0.03;};
plot_info_2 = ["cue1early","cue1late"];
plot_colors = lines(2);

fig = figure(Position=[100,100,1200,1000]);
tiled = tiledlayout(fig,3,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,3*size(plot_info_1,1));
for i=1:3*size(plot_info_1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axes([pi,pi+size(plot_info_1,1),pi+2*size(plot_info_1,1)]);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    
    plot_x = (1:sr*4)/sr-1;
    tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
    imagesc(this_axes(2),plot_x,1:nt1,permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
    this_axes(2).YDir = "reverse";
    colormap(this_axes(2),common_functions.redblue());
    colorbar(this_axes(2));
    tmp = cwa.(plot_info_2(2)).(this_info{1}).activity(:,this_info{2}+3,:); nt2 = size(tmp,3);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));
    imagesc(this_axes(3),plot_x,1:nt2,permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
    this_axes(3).YDir = "reverse";
    colormap(this_axes(3),common_functions.redblue());
    colorbar(this_axes(3));
    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-0.5,2]);
    for ax=this_axes([2,3])
        xline(ax,0,'--');
        xlim(ax,plot_x([1,end]));
    end
    ylim(this_axes(2),[1,nt1]);
    ylim(this_axes(3),[1,nt2]);
end
hold(axes,"off");
saveas(fig,'fig2F.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 2GIJK
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_normal",{["cue1late","cue1early"]},mouse_names,'fig2JK',tabular_sheet_name="fig2GI");


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 2L
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.gradient_test(cf,"tas_normal", {{"cue1late_cue1early",1,1;"cue1late_cue1early",2,2;},{"cue2late_cue2early",1,1;"cue2late_cue2early",2,2;},},mouse_names,...
    {{"cue1_early_pk";"cue2_early_pk";},{"cue1_dip";"cue2_dip";},},"cue_learning",'fig2L');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 2M
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.gradient_test(cf,"tas_normal", {{"cue1late_cue1early",1,1;"cue2late_cue2early",1,1;},{"cue1late_cue1early",2,2;"cue2late_cue2early",2,2;},{"cue1late_cue1early",3,1;"cue2late_cue2early",3,1;},},mouse_names,...
   {{"cue1_pk";"cue1_dip";"cue1_re";},{"cue2_pk";"cue2_dip";"cue2_re"},},"cue_learning_1vs2",'fig2M');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

% ███████╗██╗░░░██╗███╗░░██╗░█████╗░████████╗██╗░█████╗░███╗░░██╗░██████╗
% ██╔════╝██║░░░██║████╗░██║██╔══██╗╚══██╔══╝██║██╔══██╗████╗░██║██╔════╝
% █████╗░░██║░░░██║██╔██╗██║██║░░╚═╝░░░██║░░░██║██║░░██║██╔██╗██║╚█████╗░
% ██╔══╝░░██║░░░██║██║╚████║██║░░██╗░░░██║░░░██║██║░░██║██║╚████║░╚═══██╗
% ██║░░░░░╚██████╔╝██║░╚███║╚█████╔╝░░░██║░░░██║╚█████╔╝██║░╚███║██████╔╝
% ╚═╝░░░░░░╚═════╝░╚═╝░░╚══╝░╚════╝░░░░╚═╝░░░╚═╝░╚════╝░╚═╝░░╚══╝╚═════╝░

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
            sr = behav_functions.get_sr(behav7i);
            flag_session_with_405=0;
            loco = motion_analysis.get_overall_vel_acc(behav7i,...
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
                % if string(mouse)=="G17" && r == 41
                %     keyboard
                % end
                this_model = GLM.mdls{r};
                if flag_session_with_405 && any(contains(string(this_model.CoefficientNames),"fc405"))
                    this_glm_input_value = [ROI405i.Fc(:,r),loco.vel,loco.acc];
                    this_glm_input_name = ["fc405","vel_zero","acc_zero"];
                else
                    this_glm_input_value = [loco.vel,loco.acc];
                    this_glm_input_name = ["vel_zero","acc_zero"];
                end
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
                % coef_values_bit = ~isnan(this_model.Coefficients{:,"pValue"});
                % coef_values = this_model.Coefficients{coef_values_bit,"Estimate"};
                % coef_names = string(this_model.Coefficients.Properties.RowNames(coef_values_bit))';
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

        ta.(['file' num2str(f)]).sample_rate = behav_functions.get_sr(behav);
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

end % get_triavg end

function act_out = permute_reshape_activity(act_in)
    % reshape activity dimension from f*r*t into (f*t)*r
    n_rois = size(act_in,2);
    act_out = permute(act_in,[1,3,2]);
    act_out = reshape(act_out,[],n_rois);
end
