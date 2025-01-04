%% figure 4A-E
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
            saveas(fig_1,[cf,'sup_fig4A_E_upper.png'])
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
                saveas(fig,[cf,'sup_fig4A_E_lower.png'])
                delete(fig)
            end
        end
    end
end



%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 4F-G
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
    saveas(fig,[cf,'sup_fig4F_G_',char(phase_name),'.png'])
    delete(fig)
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 4H
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G25","G26","G27"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue1late","cue1early"]},mouse_names,'sup_fig4H',cmap_limit={[-1,1]*0.018,[-1,1]*0.017,[-1,1]*0.014});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% helper function
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
