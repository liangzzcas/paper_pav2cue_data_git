%% figure 4A-E
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
CT = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
main_sr = common_functions.get_main_samplerate();



phase_names = ["cue1early","cue1late","cue2early","cue2late"];
mouse_names = ["G17","G19","G21","G22","G23","G24"];
plot_rois = [];
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
        if ri == 1
            ylabel(axs_1(1),pattern_phase(1));
            ylabel(axs_1(2),pattern_phase(2));
        end
        k = k+1;
        if k == nr * nc + 1 || ri == length(r_to_plot)
            k = 0;
            saveas(fig_1,[cf,'decel_scatter_',char(mouse_name),'.fig'])
            saveas(fig_1,[cf,'decel_scatter_',char(mouse_name),'.png'])
            delete([fig_1])
        end
    end
end

% --------trace plot-------- %
plot_rois = [];

output_bkup = struct;
seperated_struct = struct;

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

    if contains(phase_name,"cue")
        vel_window = [0,0.6];
    elseif contains(phase_name,["rew","unpred"])
        vel_window = [-0.2,0.6];
    end
    high_low_decel_thres = vel_thres_quantile.(phase_name).across;

    % seperate high low decel trials
    for mouse_name = mouse_names
        seperated_struct.(mouse_name).(phase_name) = this_cwa.(mouse_name).activity;
    end
end

    % plot seperated trials
    plot_colors = parula(length(phase_names));
    for mouse_name = mouse_names
        this_sr = main_sr(mouse_name);
        plot_x = [1:size(seperated_struct.(mouse_name).(phase_names(1)),1)]/this_sr-1;

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
            lgs = gobjects(1,length(phase_names));
            for i = 1:length(phase_names)
                phase_name = phase_names(i);
                if r~=3
                    tmp_act_hi = seperated_struct.(mouse_name).(phase_name)(:,r,:);
                else
                    tmp_act_hi = squeeze(seperated_struct.(mouse_name).(phase_name)(:,2,:));
                    tmp_act_hi = diff(tmp_act_hi,1,1)*main_sr(mouse_name);
                    tmp_act_hi = cat(1,tmp_act_hi(1,:),tmp_act_hi);
                    tmp_act_hi = smoothdata(tmp_act_hi,1,vel_smooth_info{1},round(vel_smooth_info{2}*this_sr));
                end
                if ~isempty(tmp_act_hi)
                    lgs(i) = common_functions.plot_data_single(ax,plot_x,tmp_act_hi,plot_color = plot_colors(i,:)*0.8,sem_color = plot_colors(i,:));
                end
            end
            hold(ax,"off")
            xline(ax,0)
            yline(ax,0)
            legend(lgs,phase_names,AutoUpdate="off")
            title(ax,title_text(ri))
            k = k+1;
            if k == nr * nc + 1 || ri == length(r_to_plot)
                k = 0;
                saveas(fig,[cf,'decel_trace_',char(mouse_name),'.fig'])
                saveas(fig,[cf,'decel_trace_',char(mouse_name),'.png'])
                delete(fig)
            end
        end
    end