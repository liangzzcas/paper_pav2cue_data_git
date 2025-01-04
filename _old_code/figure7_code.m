%% figure 7B-C
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\DA_components_window_activity_filtered.mat']); % this is not rebased
plot_info_1 = {"DL21",-2,30;"DL21",25,30;};
plot_imagesc_clim = {[-1,1]*0.05;[-1,1]*0.05;};

plot_info_2 = ["post","LED_omi","Tone_omi"];
plot_info_2_text = ["post","LED omi","Tone omi"];
plot_colors = lines(length(plot_info_2));

fig = figure(Position=[100,100,600,900]);
nr=length(plot_info_2)+1;
nc = size(plot_info_1,1);
tiled = tiledlayout(fig,nr,size(plot_info_1,1),TileSpacing="tight");
axs=gobjects(1,nr*size(plot_info_1,1));
for i=1:nr*size(plot_info_1)
    axs(i) = nexttile(tiled,i);
end
hold(axs,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axs((0:nr-1)*nc+pi);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    plot_x = (1:sr*4)/sr-1;
    nts = nan(1,nr-1);
    for nri=1:nr-1
        tmp = cwa.(this_info{1}).(plot_info_2(nri)).cueOn.cue1.activity(:,this_info{2}+3,:);
        tmp = rebase_act(tmp,sr);
        if this_info{2}+3 == 1 % if it's licking, normalize it
            tmp = tmp/max(tmp,[],"all","omitmissing");
            redblue_flag = false;
            imageclim = [0,1];
        else
            redblue_flag = true;
            imageclim = plot_imagesc_clim{pi};
        end
        common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(nri,:)*0.8,sem_color=plot_colors(nri,:));
        if nri~=1
            tmp = cwa.(this_info{1}).(plot_info_2(nri)+"_complete").cueOn.cue1.activity(:,this_info{2}+3,:);
            tmp = rebase_act(tmp,sr);
        end
        nts(nri) = size(tmp,3);
        imagesc(this_axes(nri+1),plot_x,1:nts(nri),permute(tmp,[3,1,2]),imageclim);
        this_axes(nri+1).YDir = "reverse";
        if redblue_flag
            colormap(this_axes(nri+1),common_functions.redblue());
        end
        colorbar(this_axes(nri+1));
    end

    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-1,2]);
    for axi=2:nr
        ax=this_axes(axi);
        xline(ax,0,'-');
        xlim(ax,plot_x([1,end]));
        ylim(ax,[1,nts(axi-1)]);
    end
end
cla(axs(1));

% replace plot 1 with Tone (only for paper figure purpose)
lgs = gobjects(1,3);
for pi=1:1
    this_axes = axs((0:nr-1)*nc+pi);
    this_info = plot_info_1(2,:);
    sr = this_info{3};
    plot_x = (1:sr*4)/sr-1;
    for nri=1:nr-1
        tmp = cwa.(this_info{1}).(plot_info_2(nri)).cueOn.cue2.activity(:,this_info{2}+3,:);
        tmp = rebase_act(tmp,sr);
        lgs(nri)=common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(nri,:)*0.8,sem_color=plot_colors(nri,:));
    end

    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-1,2]);
end
tmp = plot_info_2;
legend(axs(1),lgs,plot_info_2_text,AutoUpdate="off");
title(axs(1),"Tone Response");
title(axs(2),"LED Response");
hold(axs,"off");
saveas(fig,'fig7B_C.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7D
close all;clear;clc;
cf = [pwd,'\'];
CT = readtable([cf,'raw_data\CT_across_DLXX_mice.xlsx']); % rebound_region_threshold: AP>0,ML<2,DV<4
cwa = load([cf,'processed_and_organized_data\DA_components_window_activity_filtered.mat']);
tas = load([cf,'processed_and_organized_data\DA_tri_avg_single_filtered_aDMs_only.mat']);

CT_AP_ML_DV = [0,2,4];
include_mask = CT{:,"fiber_bottom_AP"}>CT_AP_ML_DV(1) & CT{:,"fiber_bottom_ML"}<CT_AP_ML_DV(2) & CT{:,"fiber_bottom_DV"}<CT_AP_ML_DV(3);
CT(:,"fiber_include_bitmask") = num2cell(include_mask);
include_info = [sum(include_mask)/size(CT,1),sum(include_mask),size(CT,1)];
main_sr = 30;
n_rois = sum(include_mask);
comp_text = ["pk","dp"];
p_text = ["post","10%reward","90%reward"];

%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = cell(1,6);
%--------------------------------------------------- output file ---------------------------------------------------%

fig = figure(position=[100,100,1600,900]);
tiled = tiledlayout(fig,2,2,TileSpacing="tight");
axs = gobjects(2,2);
for i=1:4
    axs(i) = nexttile(tiled,i);
end

for ci = 1:2
    c = "cue"+ci;
    if ci == 1
        p_names = ["post","LED_omi","Tone_omi"];
    else
        p_names = ["post","Tone_omi","LED_omi"];
    end
    
    test_values_pk = [];
    test_values_dp = [];
    for pi = 1:length(p_names)
        p_name = p_names(pi);
        mu1 = tas.(p_name).(c).pk_loc(1,:,1);
        sig1 = tas.(p_name).(c).pk_loc(1,:,3);
        % mu1(~sig1) = nan;
        test_values_pk = cat(2,test_values_pk,mu1');
        mu2 = tas.(p_name).(c).pk_loc(2,:,1);
        sig2 = tas.(p_name).(c).pk_loc(2,:,3);
        % mu2(~sig2) = nan;
        test_values_dp = cat(2,test_values_dp,mu2');
        %--------------------------------------------------- output file ---------------------------------------------------%
        output_tbl(:,pi+3*(ci-1)) = {mu2'};
        %--------------------------------------------------- output file ---------------------------------------------------%
    end
    plot_x = repmat([1,2,3],[n_rois,1]);
    % hypothesis test
    test_fhandle = @kruskalwallis;
    [p_pk,tbl_pk,stat_pk] = test_fhandle(test_values_pk(:),plot_x(:),'off');
    [p_dp,tbl_dp,stat_dp] = test_fhandle(test_values_dp(:),plot_x(:),'off');
    tmp = multcompare(stat_pk,Display="off"); multcomp_pk = tmp(:,6);
    tmp = multcompare(stat_dp,Display="off"); multcomp_dp = tmp(:,6);

    ax1 = axs(2*(ci-1)+1); ax2 = axs(2*(ci-1)+2);
    title(ax1,"cue"+ci+" pk");title(ax2,"cue"+ci+" dp");
    hold([ax1,ax2],"on")
    b1 = boxchart(ax1,test_values_pk,markerstyle="none");
    b2 = boxchart(ax2,test_values_dp,markerstyle="none");
    s1 = scatter(ax1,plot_x,test_values_pk,20,'k','filled','o');
    s2 = scatter(ax2,plot_x,test_values_dp,20,'k','filled','o');
    yline(ax1,0);yline(ax2,0);
    ax1.XTickLabel = p_text;ax2.XTickLabel = p_text;
    x_pos = {[1.1,1.9],[2.1,2.9],[1.1,2.9]};
    y_pos = [1,1,1.05];
    common_functions.add_significance_to_ax(ax1,x_pos,y_pos,common_functions.p_to_asterisk(multcomp_pk([1,3,2])));
    common_functions.add_significance_to_ax(ax2,x_pos,y_pos,common_functions.p_to_asterisk(multcomp_dp([1,3,2])));
    hold([ax1,ax2],"off")
    delete(ax1)
end
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = table(output_tbl{:},VariableNames=["LED_"+["post","extinct_10","extinct_90"],"Tone_"+["post","extinct_10","extinct_90"]]);
writetable(output_tbl,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot_2\fig7_stats.xlsx'],Sheet="7d_dip_amp",WriteRowNames=true,WriteMode="overwritesheet");
%--------------------------------------------------- output file ---------------------------------------------------%
saveas(fig,[cf,'fig7D.png'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7E
close all;clear;clc;
cf = [pwd,'\'];
Ach_tas_ori = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
Ach_cwa_ori = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
CT_ach = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

DA_cwa_ori = load([cf,'processed_and_organized_data\DA_components_window_activity_filtered.mat']);
DA_tas_ori = load([cf,'processed_and_organized_data\DA_tri_avg_single_filtered_aDMs_only.mat']);
CT_da = readtable([cf,'raw_data\CT_across_DLXX_mice.xlsx']);

ach_mouse_names = ["G15","G12","G17","G19","G22","G21","G23","G24"];
ach_omi_texts = ["LEDomi","Toneomi"];
da_mouse_names = ["DL18","DL20","DL21","DL23"];
da_omi_texts = ["LED_omi","Tone_omi"];

Ach_data = struct;
DA_data = struct;
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = cell(1,2);
%--------------------------------------------------- output file ---------------------------------------------------%
for ci=1:2
    c="cue"+ci;
    % ACh
    Ach_data.(c).tas = [];
    Ach_data.(c).cwa = [];
    for mouse_name = ach_mouse_names
        this_CT = CT_ach(CT_ach{:,"mouse_name"}==mouse_name,:);
        include_mask = common_functions.get_aDMS_fiber(this_CT);
        this_tas = Ach_tas_ori.(c+ach_omi_texts(ci)).(mouse_name);
        tmp = cat(3,this_tas.mu_location_value,this_tas.mu_location_value_significance);
        Ach_data.(c).tas = cat(2,Ach_data.(c).tas,tmp(:,include_mask,:));
        this_mu = Ach_cwa_ori.(c+ach_omi_texts(ci)).(mouse_name).mu(:,logical([0;0;0;include_mask]));
        if size(this_mu,1) ~= 18*4
            x1 = 1:size(this_mu,1);
            x2 = linspace(1,size(this_mu,1),18*4);
            this_mu = interp1(x1,this_mu,x2);
        end
        Ach_data.(c).cwa = cat(2,Ach_data.(c).cwa,this_mu);
    end

    % DA
    DA_data.(c).tas = [];
    DA_data.(c).cwa = [];
    DA_data.(c).tas = DA_tas_ori.(da_omi_texts(ci)).(c).pk_loc;
    DA_mouse_names = [];
    for mouse_name = da_mouse_names
        this_CT = CT_da(CT_da{:,"mouse_name"}==mouse_name,:);
        include_mask = common_functions.get_aDMS_fiber(this_CT);
        tmp = DA_cwa_ori.(mouse_name).(da_omi_texts(ci)).cueOn.(c).activity(:,logical([0;0;0;include_mask]),:);
        tmp = mean(tmp,3,"omitmissing");
        DA_data.(c).cwa = cat(2,DA_data.(c).cwa,tmp);
        DA_mouse_names = [DA_mouse_names,repmat(mouse_name,[1,size(tmp,2)])];
    end
end

% plot data
title_texts = ["cue1","cue1","cue2","cue2"];
fig = figure(position=[100,100,1600,900]);
tiled = tiledlayout(fig,2,2,TileSpacing="compact");
sgtitle(tiled,"task")
axs=gobjects(1,4);
for i=1:4
    axs(i)=nexttile(tiled,i);
    title(axs(i),title_texts(i));
end
hold(axs,"on")

% find largest DA larges dip rois
this_da_amp = DA_data.cue1.tas(2,:,1);
[x,y] = sort(this_da_amp,"ascend");
largest_n = 5;
largest_id = y(1:largest_n);

for ci=1:2
    c="cue"+ci;
    ax1=axs(2*(ci-1)+1);ax2=axs(2*(ci-1)+2);

    this_ach_tas = Ach_data.(c).tas(3,:,2)-1;
    this_ach_cwa = Ach_data.(c).cwa;

    this_da_tas = DA_data.(c).tas(2,:,2)/30-1;
    this_da_cwa = DA_data.(c).cwa;
    % % filter sig
    % sig_bit = DA_data.(c).tas(2,:,3);
    % this_da_tas(~sig_bit) = nan;
    % this_da_cwa(:,~sig_bit) = nan;
    % % filter amplitude
    this_da_tas = this_da_tas(largest_id);
    this_da_cwa = this_da_cwa(:,largest_id);
    %--------------------------------------------------- output file ---------------------------------------------------%
    if ci==1
    output_tbl{1} = this_da_tas';
    output_tbl{2} = this_ach_tas';
    end
    %--------------------------------------------------- output file ---------------------------------------------------%
    this_da_mouse_names = DA_mouse_names(largest_id);

    DA_Ach_lag = mean(this_ach_tas,"all","omitmissing") - mean(this_da_tas,"all","omitmissing");
    ax1.Title.String = ax1.Title.String+" "+sprintf("%0.2fs",DA_Ach_lag);

    test_data = [this_da_tas,this_ach_tas];
    test_data_tag = [repmat("DA",[1,length(this_da_tas)]),...
        repmat("ACh",[1,length(this_ach_tas)])];
    
    test_fhandle = @kruskalwallis;
    [p_pk,tbl_pk,stat_pk] = test_fhandle(test_data,test_data_tag,'off');
    tmp = multcompare(stat_pk,Display="off"); multcomp_pk = tmp(:,6);
    boxchart_order = ["DA","ACh"];
    x_pos = {[1.1,1.9],[2.1,2.9],[1.1,2.9]};
    y_pos = [1,1,1.05];
    

    test_data_tag = test_data_tag(1:end-2);
    boxchart_order = boxchart_order(1:2);
    test_data = test_data(1:end-2);
    x_pos = {[1.1,1.9]};
    y_pos = 1.05;
    multcomp_pk = multcomp_pk(1);
    b1 = boxchart(ax1,categorical(test_data_tag,boxchart_order),test_data,markerstyle="none");
    s1 = scatter(ax1,1,this_da_tas,20,'k','filled','o');
    s2 = scatter(ax1,2,this_ach_tas,20,'k','filled','o');
    common_functions.add_significance_to_ax(ax1,x_pos,y_pos,common_functions.p_to_asterisk(multcomp_pk));
    
    sensor_max = nan(1,2);
    sensor_max(1) = max(this_ach_cwa(1:18*3,:),[],"all");
    sensor_max(2) = max(this_da_cwa(1:30*3,:),[],"all");
    sensor_scale =  1./(sensor_max/sensor_max(2));
    
    tmp = this_ach_cwa(1:18*3,:); tmp = tmp-mean(tmp(1:18,:),"all","omitmissing"); tmp = tmp*sensor_scale(1);
    yyaxis(ax2,"right");
    p1=common_functions.plot_data_single(ax2,(1:18*3)/18-1,tmp,plot_color=[.8,0,0],sem_color=[1,0,0]);
    yline(ax2,0,'r--'); ylim(ax2,[-0.02,0.02]);
    yyaxis(ax2,"left");
    tmp = this_da_cwa(1:30*3,:); tmp = tmp-mean(tmp(1:30,:),"all","omitmissing"); tmp = tmp*sensor_scale(2);
    p2=common_functions.plot_data_single(ax2,(1:30*3)/30-1,tmp,plot_color=[0,.8,0],sem_color=[0,1,0]);
    yline(ax2,0,'g--'); ylim(ax2,[-0.03,0.045]);
    legend(ax2,[p1,p2],"normalized "+["Ach","DA"],AutoUpdate="off");
    ylabel(ax2,"random unit")
    xline(ax2,0);
end
hold(axs,"off")
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = random_functions.pad_cell_for_table(output_tbl);
output_tbl = table(output_tbl{:},VariableNames=["DA","ACh"]);
writetable(output_tbl,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot_2\fig7_stats.xlsx'],Sheet="7e_DAACh_timing",WriteRowNames=true,WriteMode="overwritesheet");
%--------------------------------------------------- output file ---------------------------------------------------%
delete(axs([3,4]))
saveas(fig,[cf,'fig7E.png'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7F
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\DA_components_window_activity_filtered.mat']);
tas = load([cf,'processed_and_organized_data\DA_tri_avg_single_filtered_aDMs_only.mat']);
CT = readtable([cf,'raw_data\CT_across_DLXX_mice.xlsx']); % rebound_region_threshold: AP>0,ML<2,DV<4
lick_info = load([cf,'DA_across_mice_lick_index_data_whole_ITI.mat']);
mouse_names_ori = ["DL18","DL20","DL21","DL23"];

include_mask = CT{:,"fiber_bottom_AP"}>0 & CT{:,"fiber_bottom_ML"}<2 & CT{:,"fiber_bottom_DV"}<4;
CT(:,"fiber_include_bitmask") = num2cell(include_mask);
include_info = [sum(include_mask)/size(CT,1),sum(include_mask),size(CT,1)];

pnames_to_plot = ["pre","post","LED_omi","Tone_omi"];
main_sr = 30;

n_rois = sum(include_mask);
pk_dp_window_s = common_functions.DA_get_comp_timewindow();
pk_dp_window_frame = (pk_dp_window_s+1)*main_sr;
dp_window_frame = pk_dp_window_frame(2,1):pk_dp_window_frame(2,2);

% put all ROIs togather
mouse_names = string(fields(cwa))';
for ci=1:2
    c = "cue"+ci;
    if ci==1
        pname="LED_omi";
    else
        pname="Tone_omi";
    end

    across_mouse_frame_diff = struct;
    for mouse_name = mouse_names
        this_CT = CT(CT{:,"mouse_name"}==mouse_name,:);
        this_CT_bitmask = this_CT{:,"fiber_include_bitmask"};
        this_act = cwa.(mouse_name).(pname+"_complete").cueOn.(c).activity;
        this_act = this_act(:,logical([0;0;0;this_CT_bitmask]),:);
        this_n_rois = size(this_act,2); this_n_trial = size(this_act,3);
        mouse_data = nan(this_n_rois,this_n_trial);
        for r=1:this_n_rois
            for t=1:this_n_trial
                this_single = this_act(:,r,t);
                [dip,loc] = findpeaks(-this_single(dp_window_frame),dp_window_frame,NPeaks=1,SortStr='descend');
                if isempty(loc)
                    loc = nan;
                end
                mouse_data(r,t) = loc;
            end
        end

        trial_batch_size=1;
        batch_lick_info = [];
        this_m = mouse_name;
        mouse_day_info = common_functions.get_include_days_DA(this_m);
        phase = pname;
        this_day = mouse_day_info{phase+"_complete"};

        batch_lick_across_day = [];
        batch_lick_across_day_n = nan(1,length(this_day));
        for d_i = 1:length(this_day)
            d = this_day(d_i);
            this_lick = lick_info.(this_m).single_trial_taking_all_ITI_frame.("day"+d).("cue"+ci+"_lick");
            n_batch = ceil(size(this_lick,1)/trial_batch_size);
            batch_trial_id = [0:n_batch]*trial_batch_size; 
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
        batch_lick_info = cat(1,batch_lick_info,batch_lick_across_day);
        
        this_rois = 1:this_n_rois;

        plot_color = lines(4);
        lick_index_type = 3;
        smooth_detail = {30,'lowess',0};
        cusum_detail = {5,5};
        diff_chang_point_across = {};
        
        re_single = mouse_data';
        lick_single = batch_lick_info(:,3);
        % smooth data
        lick_smoothed = smooth(lick_single,smooth_detail{:});
        re_smoothed = nan(size(re_single));
        for ri = 1:size(re_single,2)
            re_smoothed(:,ri) = smooth(re_single(:,ri),smooth_detail{:});
        end
        
        % change detection
        real_control = 5;
        cusum_lick_mean = mean(lick_smoothed(1:real_control),"omitmissing");
        cusum_lick_std = std(lick_smoothed(1:real_control),[],"omitmissing");
        cusum_re_mean = mean(re_smoothed(1:real_control,:),1,"omitmissing");
        cusum_re_std = std(re_smoothed(1:real_control,:),1,"omitmissing");
        % if no previous sessions are used, no need for rigid threshold
        real_control = 5;
        control_n_trial = 0;
        
        [~,~,lick_uppersum,lick_lowersum] = cusum(lick_smoothed,cusum_detail{:},cusum_lick_mean,cusum_lick_std);
        cusum_lick_control = min([lick_lowersum(1:real_control);-cusum_detail{1}*cusum_lick_std],[],"all");
        lick_smooth_ilower = find(lick_lowersum<cusum_lick_control,1,"first");

        change_point = cell(1,length(this_rois)+1); % {lick_change,roi_changes}
        change_point{1} = lick_smooth_ilower;
        for ri = 1:length(this_rois)
            r = ri;
            [~,~,re_uppersum,re_lowersum] = cusum(re_smoothed(:,ri),3,1,cusum_re_mean(ri),cusum_re_std(ri));
            re_cusum_control = min([re_lowersum(1:real_control);-cusum_detail{1}*cusum_re_std(ri)],[],"all");
            re_smooth_iupper = find(re_lowersum<re_cusum_control,1,"first");
            change_point{ri+1} = re_smooth_iupper;
        end
        diff_chang_point_across = cat(2,diff_chang_point_across,cellfun(@(x) x-change_point{1},change_point(2:end),UniformOutput=false));
        across_mouse_frame_diff.(mouse_name) = diff_chang_point_across;
    end

    tmps = string(fields(across_mouse_frame_diff))';
    hist_data = [];
    for tmp=tmps
        tmp1 = across_mouse_frame_diff.(tmp);
        tmp1 = tmp1(~cellfun(@isempty,tmp1));
        hist_data = cat(2,hist_data,cell2mat(tmp1));
    end

    % histogram plot
    fig = figure(Position=[100,100,800,600]);
    sgtitle(fig,c+" dip vs lick");
    ax = axes(fig);
    histogram(ax,hist_data,Normalization="count",BinWidth=10);
    xline(ax,median(hist_data,"omitmissing"),'k-',LineWidth=1.5)
    thisxlim = xlim(ax);
    xlim(ax,[-1,1]*max(abs(thisxlim)));
    saveas(fig,[cf,'fig7E_',char(c),'.png'])
    delete(fig)
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7G
close all;clear;clc;
cf = [pwd,'\'];
DA_task = load([cf,'processed_and_organized_data\DA_components_window_activity_filtered.mat']);

mouse_names = ["DL18","DL20","DL21","DL23"];
sensor_names = ["DA","DA","DA","DA"];

phase_names = ["post","LED_omi","Tone_omi"];
phase_texts = ["post","LED ext","Tone ext"];

% parameters
vel_window = [0,1];
mouse_struct = struct;

for m_i = 1:length(mouse_names)
    mouse_name = mouse_names(m_i);
    sensor_name = sensor_names(m_i);

    lick_data = DA_task;
    main_sr = 30;

    vel_window_frame = (vel_window+1)*main_sr;

    mouse_struct.(mouse_name) = nan(2,3,2); % 2 cues * 3 phase * 2 types(acc, diff_v)
    for ci=1:2
        for pi=1:3
            phase_name = phase_names(pi);
            if ~isfield(lick_data.(mouse_name),phase_name)
                continue
            end
            vel_acc = mean(lick_data.(mouse_name).(phase_name).cueOn.("cue"+ci).activity(:,2:3,:),3,"omitmissing");
            vel_window_frames = vel_window_frame(1):vel_window_frame(2);
            [M,I] = max(abs(vel_acc(vel_window_frames,2)),[],"omitmissing");
            I = I+vel_window_frame(1)-1;
            mouse_struct.(mouse_name)(ci,pi,1) = M;
            mouse_struct.(mouse_name)(ci,pi,2) = max(vel_acc(vel_window_frames,1))-min(vel_acc(vel_window_frames,1));
        end
    end
end

%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = cell(1,12);
output_tbl_varnames = reshape(["peak_decel_light_","peak_decel_tone_","vel_change_light_","vel_change_tone_",]+["post";"ex_10";"ex_20"],1,[]);
%--------------------------------------------------- output file ---------------------------------------------------%

test_name = "kruskalwallis";
test_fh = @kruskalwallis;
movement_comp_names = ["acc","diff vel"];
fig = figure(position=[100,100,1600,900]);
tiled = tiledlayout(fig,2,2,TileSpacing="tight");
sgtitle(tiled,"testing acceleration velocity ("+test_name+")")
axs = gobjects(1,4);
for i=1:4
    axs(i)=nexttile(tiled,i);
end
hold(axs,"on")
for ci=1:2
    for axi=1:2
        ax = axs(2*(ci-1)+axi);
        this_data = cell2mat(structfun(@(x) {x(ci,:,axi)},mouse_struct));
        %--------------------------------------------------- output file ---------------------------------------------------%
        axitexts = ["peak_decel_","vel_change_"];
        citexts = ["light_","tone_"];
        this_col_ids = contains(output_tbl_varnames,axitexts(axi)+citexts(ci));
        output_tbl(this_col_ids) = num2cell(this_data,1);
        %--------------------------------------------------- output file ---------------------------------------------------%
        [p,tbl,stats] = test_fh(this_data,[],"off");
        c = multcompare(stats,display="off");
        multi_p = c(:,6)';
        asters = common_functions.p_to_asterisk(multi_p);

        this_test_data_tag = categorical([repmat(phase_texts(1),length(mouse_names),1);repmat(phase_texts(2),length(mouse_names),1);repmat(phase_texts(3),length(mouse_names),1)],phase_texts);
        boxchart(ax,this_test_data_tag,reshape(this_data,[],1),MarkerStyle='o');
        plot(ax,[1,2,3],this_data','ko-');
        common_functions.add_significance_to_ax(ax,{[1.1,1.9],[2.1,2.9],[1.1,2.9]},[1,1,1.05],asters([1,3,2]));
        title(ax,"cue "+ci+" "+movement_comp_names(axi)+" p="+p);
    end    
end
hold(axs,"off")
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = table(output_tbl{:},VariableNames=output_tbl_varnames);
writetable(output_tbl,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot_2\fig7_stats.xlsx'],Sheet="7g_vel_decel",WriteRowNames=true,WriteMode="overwritesheet");
%--------------------------------------------------- output file ---------------------------------------------------%
saveas(fig,[cf,'fig7G.png'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7H
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\DA_ITI_lick_ta.mat']).filtered;

main_sr = 30;
m_name = "DL21";
plot_x = (1:2*main_sr)/main_sr-1;
el_tag = ["pre","post"];
    
% trace plot
plot_color = lines(2);
lgs = gobjects(1,2);
fig = figure(Position=[100,100,1600,500]);
tiled = tiledlayout(1,3,Parent=fig,TileSpacing="tight");
sgtitle(tiled,m_name+" ITI licking average")

ax = nexttile(tiled,1);
hold(ax,"on")
for el = 1:2
    plot_act = cwa.(m_name).(el_tag(el))(:,2+3,:);
    lgs(el)=common_functions.plot_data_single(ax,plot_x,plot_act,plot_color=plot_color(el,:)*0.8,sem_color=plot_color(el,:));
end
hold(ax,"off")
yline(ax,0);xline(ax,0,'k--');
legend(ax,lgs,["pre-learning","post-learning"]);
ylabel(ax,"Mean \DeltaF/F");xlabel(ax,"Time(s) from lick onset");

% imagesc plot
ax3 = nexttile(tiled,2); ax4 = nexttile(tiled,3);
hold([ax3,ax4],"off")
act1 = permute(cwa.(m_name).pre(:,2+3,:),[1,3,2]); act2 = permute(cwa.(m_name).post(:,2+3,:),[1,3,2]);
n_trial = [size(act1,2),size(act2,2)];

tmp_clim = max(abs([prctile(act1,[2,98],"all"),prctile(act2,[2,98],"all")]),[],"all");
imagesc(ax3,plot_x,1:n_trial(1),act1',[-1,1]*tmp_clim);
imagesc(ax4,plot_x,1:n_trial(2),act2',[-1,1]*tmp_clim);
colormap(ax3,"redblue");colormap(ax4,"redblue");
xline(ax3,0);xline(ax4,0);
hold([ax3,ax4],"off")
colorbar(ax3);colorbar(ax4);
title(ax3,"pre-learning"); title(ax4,"post-learning");
saveas(fig,[cf,'fig7H.png'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7I
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\DA_ITI_lick_ta.mat']).filtered;
CT = readtable([cf,'raw_data\CT_across_DLXX_mice.xlsx']); % rebound_region_threshold: AP>0,ML<2,DV<4

include_mask = CT{:,"fiber_bottom_AP"}>0 & CT{:,"fiber_bottom_ML"}<2 & CT{:,"fiber_bottom_DV"}<4;
CT(:,"fiber_include_bitmask") = num2cell(include_mask);
include_info = [sum(include_mask)/size(CT,1),sum(include_mask),size(CT,1)];
n_rois = sum(include_mask);
el_tag = ["pre","post"];

mouse_names_ori = string(fields(cwa)');
main_sr = 30;
plot_x = (1:2*main_sr)/main_sr-1;
plot_color = lines(2);

across_mouse_data = struct;
for el = 1:2
    across_mouse_data.(el_tag(el)) = {};
end
for mi = 1:length(mouse_names_ori)
    m_name = mouse_names_ori(mi);
    this_CT = CT(CT{:,"mouse_name"}==m_name,:);
    this_bitmask = this_CT{:,"fiber_include_bitmask"};
    this_bitmask_3_roi_only = cat(1,[false;false;false],this_bitmask);
    
    n_rois = size(cwa.(m_name).pre,2)-3;
    for el = 1:2
        tmp = cwa.(m_name).(el_tag(el));
        across_mouse_data.(el_tag(el)) = cat(2,across_mouse_data.(el_tag(el)),num2cell(tmp(:,this_bitmask_3_roi_only,:),3));
    end
end

% plot tri_ave across all ROIs
across_all = cell(2,1);
for el=1:2
    tmp = across_mouse_data.(el_tag(el));
    tmp = cellfun(@(x) mean(x,"all","omitmissing"),tmp);
    across_all(el) = {tmp};
end

fig = figure(position=[100,100,800,600]);
sgtitle(fig,"average across all aDS ROIs")
ax = axes(fig);
hold(ax,"on")
p1=common_functions.plot_data_single(ax,plot_x,across_all{1},plot_color=plot_color(1,:)*0.8,sem_color=plot_color(1,:));
p2=common_functions.plot_data_single(ax,plot_x,across_all{2},plot_color=plot_color(2,:)*0.8,sem_color=plot_color(2,:));
hold(ax,"off")
xline(ax,0,'k--');yline(ax,0);
legend(ax,[p1,p2],el_tag+"-learning",AutoUpdate="off")
ylabel(ax,"Mean \DeltaF/F");xlabel(ax,"Time(s) from lick onset")
saveas(fig,[cf,'fig7I.png'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7J
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\DA_ITI_lick_tas.mat']);
CT = readtable([cf,'raw_data\CT_across_DLXX_mice.xlsx']); % rebound_region_threshold: AP>0,ML<2,DV<4
include_mask = CT{:,"fiber_bottom_AP"}>0 & CT{:,"fiber_bottom_ML"}<2 & CT{:,"fiber_bottom_DV"}<4;
CT(:,"fiber_include_bitmask") = num2cell(include_mask);
include_info = [sum(include_mask)/size(CT,1),sum(include_mask),size(CT,1)];
n_rois = sum(include_mask);

mouse_names = ["DL18","DL20","DL21","DL23"];
el_text = ["pre","post"];

tas_across = struct;
for el = el_text
    tas_across.(el).mu = [];
    tas_across.(el).sig = [];
end
for mi = 1:length(mouse_names)
    mouse_name = mouse_names(mi);
    this_CT = CT(CT{:,"mouse_name"}==mouse_name,:);
    for el = el_text
        tmp = tas.(mouse_name).(el).mu_location_value(:,this_CT{:,"fiber_include_bitmask"},:);
        tas_across.(el).mu = cat(2,tas_across.(el).mu,tmp);
        tmp = tas.(mouse_name).(el).significance(:,this_CT{:,"fiber_include_bitmask"},:);
        tas_across.(el).sig = cat(2,tas_across.(el).sig,tmp);
    end
end

%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = cell(1,2);
output_tbl_varnames = ["pre","post"];
%--------------------------------------------------- output file ---------------------------------------------------%

title_text = ["pk","dip"];
fig = figure(position=[100,100,1200,600]);
tiled = tiledlayout(fig,1,2,TileSpacing="tight");
axs = gobjects(1,2);
for i=1:2
    axs(i)=nexttile(tiled,i);
    title(axs(i),title_text(i));
end
plot_x = repmat([1,2],[n_rois,1]);
mu_pk_e = tas_across.pre.mu(1,:,1)';
mu_pk_l = tas_across.post.mu(1,:,1)';
mu_pk = [mu_pk_e,mu_pk_l];
mu_dp_e = tas_across.pre.mu(2,:,1)';
mu_dp_l = tas_across.post.mu(2,:,1)';
mu_dp = [mu_dp_e,mu_dp_l];
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = num2cell(mu_dp,1);
%--------------------------------------------------- output file ---------------------------------------------------%
test_fhandle = @kruskalwallis;
[p_pk,tbl_pk,stat_pk] = test_fhandle(mu_pk(:),plot_x(:),'off');
[p_dp,tbl_dp,stat_dp] = test_fhandle(mu_dp(:),plot_x(:),'off');
tmp = multcompare(stat_pk,Display="off"); multcomp_pk = tmp(:,6);
tmp = multcompare(stat_dp,Display="off"); multcomp_dp = tmp(:,6);

ax1 = axs(1); ax2 = axs(2);
hold([ax1,ax2],"on")
b1 = boxchart(ax1,mu_pk,markerstyle="none");
b2 = boxchart(ax2,mu_dp,markerstyle="none");
s1 = scatter(ax1,plot_x,mu_pk,20,'k','filled','o');
s2 = scatter(ax2,plot_x,mu_dp,20,'k','filled','o');
ax1.XTickLabel = el_text+"-learning"; ax2.XTickLabel = el_text+"-learning";
x_pos = {[1.1,1.9]};
y_pos = [1.1];
common_functions.add_significance_to_ax(ax1,x_pos,y_pos,random_functions.p_to_asterisk(multcomp_pk));
common_functions.add_significance_to_ax(ax2,x_pos,y_pos,random_functions.p_to_asterisk(multcomp_dp));
yline(ax1,0);yline(ax2,0);
hold([ax1,ax2],"off")
ylabel(ax1,"Mean \DeltaF/F dip"); ylabel(ax2,"Mean \DeltaF/F dip");
%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = table(output_tbl{:},VariableNames=output_tbl_varnames);
writetable(output_tbl,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot_2\fig7_stats.xlsx'],Sheet="7j_dip_amp_prepost",WriteRowNames=true,WriteMode="overwritesheet");
%--------------------------------------------------- output file ---------------------------------------------------%
delete(ax1);
saveas(fig,[cf,'fig7J.png'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 7k
close all;clear;clc;
cf = [pwd,'\'];
Ach_tas_ori = load([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat']);
Ach_cwa_ori = load([cf,'processed_and_organized_data\ITI_licking_components_window_activity_filtered.mat']);
CT_ach = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);


% averaging aDMs ROIs for all DA mice and get timing of components
cwa = load([cf,'processed_and_organized_data\DA_ITI_lick_ta.mat']);
tas = load([cf,'processed_and_organized_data\DA_ITI_lick_tas.mat']);
CT_DA = readtable([cf,'raw_data\CT_across_DLXX_mice.xlsx']); % rebound_region_threshold: AP>0,ML<2,DV<4
include_mask = CT_DA{:,"fiber_bottom_AP"}>0 & CT_DA{:,"fiber_bottom_ML"}<2 & CT_DA{:,"fiber_bottom_DV"}<4;
CT_DA(:,"fiber_include_bitmask") = num2cell(include_mask);
include_info = [sum(include_mask)/size(CT_DA,1),sum(include_mask),size(CT_DA,1)];
n_rois = sum(include_mask);
el_tag = ["pre","post"];

mouse_names_ori = string(fields(cwa.filtered)');
main_sr = 30;

across_mouse_data = struct;
for el = 1:2
    across_mouse_data.(el_tag(el)) = [];
    across_mouse_data.(el_tag(el)+"timing") = [];
end
for mi = 1:length(mouse_names_ori)
    m_name = mouse_names_ori(mi);
    this_CT = CT_DA(CT_DA{:,"mouse_name"}==m_name,:);
    this_bitmask = this_CT{:,"fiber_include_bitmask"};
    this_bitmask_3 = cat(1,[true;true;true],this_bitmask);
    this_bitmask_3_roi_only = cat(1,[false;false;false],this_bitmask);
    plot_x = (1:2*main_sr)/main_sr-1;

    % ROI plot
    n_rois = size(cwa.filtered.(m_name).pre,2)-3;
    n_trial = nan(1,2);
    for el = 1:2
        tmp = cwa.filtered.(m_name).(el_tag(el));
        n_trial(el) = size(tmp,3);
        across_mouse_data.(el_tag(el)) = cat(2,across_mouse_data.(el_tag(el)),...
            mean(tmp(:,this_bitmask_3_roi_only,:),3,"omitmissing"));
        
        across_mouse_data.(el_tag(el)+"timing") = cat(2,across_mouse_data.(el_tag(el)+"timing"),...
            cat(3,tas.(m_name).(el_tag(el)).mu_location_value(:,this_bitmask,:),tas.(m_name).(el_tag(el)).significance(:,this_bitmask)));
    end
end
DA_tas_ori = across_mouse_data;


% averaging aDMs ROIs for all ACh mice and get timing of components
Ach_data = Ach_tas_ori.early;
Ach_cwa = Ach_cwa_ori.early;
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
Ach_timing = [];
Ach_trace = [];
for mouse_name = mouse_names
    this_CT = CT_ach(CT_ach{:,"mouse_name"}==mouse_name,:);
    include_mask = this_CT{:,"fiber_bottom_AP"}>0 & this_CT{:,"fiber_bottom_ML"}<2 & this_CT{:,"fiber_bottom_DV"}<4;
    tmp = cat(3,Ach_data.(mouse_name).mu_location_value,Ach_data.(mouse_name).mu_location_value_significance);
    Ach_timing = cat(2,Ach_timing,tmp(:,include_mask,:));
    
    tmp = Ach_cwa.(mouse_name).mu(:,logical([0;0;0;include_mask]));
    if Ach_cwa.(mouse_name).sr ~= 18
        x1 = 1:size(tmp,1);
        x2 = linspace(1,size(tmp,1),18*2);
        tmp = interp1(x1,tmp,x2);
    end
    Ach_trace = cat(2,Ach_trace,tmp);
end
DA_timing = DA_tas_ori.pretiming;
DA_trace = DA_tas_ori.pre;

%
plot_colors = lines(2);

DA_dip_timing = DA_timing(2,:,2);
DA_dip_sig = DA_timing(2,:,3);
Ach_dip_timing = Ach_timing(3,:,2);
Ach_dip_sig = Ach_timing(3,:,3);
% Ach_dip_timing(~Ach_dip_sig) = nan;
% DA_dip_timing(~DA_dip_sig) = nan;
boxchart_order = ["DA dip","ACh late peak"];
test_data = [DA_dip_timing,Ach_dip_timing];
test_data_tag = [repmat(boxchart_order(1),[1,length(DA_dip_timing)]),repmat(boxchart_order(2),[1,length(Ach_dip_timing)])];

test_fhandle = @kruskalwallis;
[p_pk,tbl_pk,stat_pk] = test_fhandle(test_data,test_data_tag,'off');
tmp = multcompare(stat_pk,Display="off"); multcomp_pk = tmp(:,6);

%--------------------------------------------------- output file ---------------------------------------------------%
output_tbl = cell(1,2);
output_tbl_varnames = ["DA dip","ACh re"];
output_tbl{1} = DA_dip_timing';
output_tbl{2} = Ach_dip_timing';
output_tbl = random_functions.pad_cell_for_table(output_tbl);
output_tbl = table(output_tbl{:},VariableNames=output_tbl_varnames);
writetable(output_tbl,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot_2\fig7_stats.xlsx'],Sheet="7k_DAACh_timing",WriteRowNames=true,WriteMode="overwritesheet");
%--------------------------------------------------- output file ---------------------------------------------------%

fig = figure(position=[100,100,1200,600]);
tiled = tiledlayout(fig,1,2,TileSpacing="compact");
sgtitle(tiled,"ITI licking dip timing")
ax1=nexttile(tiled,1); ax2=nexttile(tiled,2);
hold([ax1,ax2],"on")
b1 = boxchart(ax1,categorical(test_data_tag,boxchart_order),test_data,markerstyle="none");
s1 = scatter(ax1,1,DA_dip_timing,20,'k','filled','o');
s2 = scatter(ax1,2,Ach_dip_timing,20,'k','filled','o');
x_pos = {[1.1,1.9]};
y_pos = [1.05];
common_functions.add_significance_to_ax(ax1,x_pos,y_pos,random_functions.p_to_asterisk(multcomp_pk));
ylabel(ax1,"Latency (s)")

yyaxis(ax2,"left");
p1=common_functions.plot_data_single(ax2,(1:36)/18-1,Ach_trace,plot_color=[.8,0,0],sem_color=[1,0,0]);
yline(ax2,0,Color=[.8,0,0],LineStyle='--');
yyaxis(ax2,"right");
p2=common_functions.plot_data_single(ax2,(1:60)/30-1,DA_trace,plot_color=[0,.8,0],sem_color=[0,1,0]);
yline(ax2,0,Color=[0,.8,0],LineStyle='--');
hold([ax1,ax2],"off")
xline(ax2,0,'k--');
xlabel(ax2,"Time (s) from lick onset")
legend(ax2,[p1,p2],["Ach","DA"],AutoUpdate="off");
saveas(fig,[cf,'fig7K.png'])
delete(fig)









%% functions
function out = rebase_act(act_in,main_sr)
    tmp = mean(act_in(1:main_sr,:,:),[1,3],"omitmissing");
    tmp = repmat(tmp,[size(act_in,1),1,size(act_in,3)]);
    out = act_in-tmp;
end
        