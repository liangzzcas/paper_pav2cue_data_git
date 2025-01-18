%% sup figure 6A
close all;clear;clc;
cf = [pwd,'\'];

ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
main_sr = common_functions.get_main_samplerate();

comp_time_window = common_functions.get_ITI_licking_comp_timewindow();
lick_data = load([cf,'processed_and_organized_data\ITI_licking_aligned_highpass03.mat']);

null_window_s = 0.7;
include_days = common_functions.get_include_days();
plot_data = build_plotdata_by_session(lick_data,include_days,null_window_s);

axs_title_texts = ["early","late","LEDomi","Toneomi"];
component_texts = ["peak","dip"];

component_colors = lines(5);
mouse_names = string(fields(plot_data)');
mouse_names = mouse_names(~contains(mouse_names,["G25","G26","G27"]));
time_maxmin = cell(length(mouse_names),3*4);
plot_colors = lines(8);
for mi = 1:length(mouse_names)
    mouse_name = mouse_names(mi);
    this_table = ct_table(ct_table{:,"mouse_name"}==mouse_name,:);
    sig_rois = find(logical(this_table{:,"significance"})');
    this_tbl = plot_data.(mouse_name);
    session_names = string(this_tbl.Properties.RowNames');
    for i = 1:4
        s_i = find(contains(session_names,axs_title_texts(i)));
        if isempty(s_i)
            continue
        end
        roi_sr = this_tbl{s_i,1};
        roi_traces = mean(cell2mat(this_tbl{s_i,2}),3,"omitmissing");
        roi_traces = roi_traces(:,[1,2,3,sig_rois+3]);

        dip_window = round(([0,0.6]+1)*roi_sr); dip_window = dip_window(1):dip_window(2);
        peak_window = round(([-0.27,0.8]+1)*roi_sr); peak_window = peak_window(1):peak_window(2);

        [~,dipmin] = min(roi_traces(dip_window,4:end),[],1,"omitmissing"); dipmin = dipmin + dip_window(1) - 1;
        dipmin_mean = round(median(dipmin,"omitmissing"));

        [~,pkmax] = max(roi_traces(peak_window(1):dipmin_mean,4:end),[],1,"omitmissing"); pkmax = pkmax + peak_window(1) - 1;
        [~,remax] = max(roi_traces(dipmin_mean:peak_window(end),4:end),[],1,"omitmissing"); remax = remax + dipmin_mean - 1;
        time_maxmin(mi,3*i-2:3*i) = num2cell(([pkmax;dipmin;remax]')/roi_sr-1,1);
    end
end

% plot histogram of time of components for each phase
phase_xlim = [-0.5,1];
binwidth = 0.05;
cum_x = phase_xlim(1):binwidth:phase_xlim(2);
if cum_x(end) ~= phase_xlim(2)
    cum_x(end+1) = phase_xlim(2);
end

% histogram of component time ad hoc (early_late in one fig, LEDTONEomi in another fig)
fig_3s = gobjects(1,2);
tiled_3s = gobjects(1,2);
for i = 1:1
    fig_3 = figure(Position=[100,100,1600,900]);
    tiled_3 = tiledlayout(fig_3,2,2,TileSpacing="tight");
    sgtitle(tiled_3,axs_title_texts(2*i-1)+" "+axs_title_texts(2*i));
    fig_3s(i) = fig_3;
    tiled_3s(i) = tiled_3;
end

for i = 1:2
    tiled_3 = tiled_3s(ceil(i/2));
    axs = gobjects(1,2);
    for ax_i = 1:2
        axs(ax_i) = nexttile(tiled_3,2-mod(i,2)+2*ax_i-2);
        title(axs(ax_i),axs_title_texts(i)+" "+component_texts(ax_i))
        xlabel(axs(ax_i),"Time since lick onset (second)")
        ylabel(axs(ax_i),"component count")
    end
    hold(axs,"on")

    comm_frames = cell2mat(time_maxmin(:,3*i-2:3*i));
    pk_re_frames = [comm_frames(:,1);comm_frames(:,3)];
    dp_frames = comm_frames(:,2);
    tmp1 = quantile(comm_frames,0.05,1);
    tmp2 = quantile(comm_frames,0.95,1);
    
    histogram(axs(1),pk_re_frames,BinWidth=binwidth,Normalization="count",FaceColor=component_colors(1,:));
    histogram(axs(2),dp_frames,BinWidth=binwidth,Normalization="count",FaceColor=component_colors(2,:));

    if ~isempty(comp_time_window)
        %%%%%%
        % xl1_1 = xregion(axs(1),comp_time_window(1,1),comp_time_window(1,2),FaceColor=component_colors(5,:));
        % xl3_1 = xregion(axs(1),comp_time_window(3,1),comp_time_window(3,2),FaceColor=component_colors(3,:));
        % xl2_1 = xregion(axs(2),comp_time_window(2,1),comp_time_window(2,2),FaceColor='k');
        %%%%%%
        tmp = xline(axs(1),[comp_time_window(1,1),comp_time_window(1,2)],Color=plot_colors(5,:),LineStyle='-',Linewidth=1.8);
        xl1_1 = tmp(1);
        tmp = xline(axs(1),[comp_time_window(3,1),comp_time_window(3,2)],Color=plot_colors(3,:),LineStyle='-',Linewidth=1.8);
        xl3_1 = tmp(1);
        tmp = xline(axs(2),[comp_time_window(2,1),comp_time_window(2,2)],'k-',Linewidth=1.8);
        xl2_1 = tmp(1);
        %%%%%%
    end
    legend(axs(1),[xl1_1(1),xl3_1(1)],["pk1 window","pk2 window"],AutoUpdate="off");
    legend(axs(2),[xl2_1(1)],["dp window"],AutoUpdate="off");
    xlim(axs,phase_xlim);
    hold(axs,"off")
end
for i = 1:1
    fig_3 = fig_3s(i);
    tiled_3 = tiled_3s(i);
    axs = gobjects(1,4);
    for ax_i = 1:4
        axs(ax_i) = nexttile(tiled_3,ax_i);
    end
    linkaxes(axs,"xy")
    saveas(fig_3,[cf,'sup_fig6A.png'])
    delete(fig_3)
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 6B
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat']);
CT_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
phases = ["early"];
mouse_name_cell = {["G15","G12","G17","G19","G22","G21","G23","G24"]};

comp_text = ["pk","dp","re"];
for m_i = 1:size(mouse_name_cell,2)
    mouse_names = mouse_name_cell{m_i};
    save_suffix = char(strjoin(mouse_names,"_"));
    data_table_across_phase = struct;
    for phase = phases
        data_merged = [];
        table_merged = [];
        for mouse_name = mouse_names
            tmp_data = tas.(phase).(mouse_name).mu_location_value;
            tmp_sig = logical(tas.(phase).(mouse_name).mu_location_value_significance);
            tmp_data(~repmat(tmp_sig,[1,1,2])) = nan;
            tmp_table = CT_table(string(CT_table{:,"mouse_name"})==mouse_name,:);
            tmp_sig = tmp_table{:,"significance"}==1;
            tmp_data = tmp_data(:,tmp_sig,:);
            tmp_table = tmp_table(tmp_sig,:);
            data_merged = cat(2,data_merged,tmp_data);
            table_merged = cat(1,table_merged,tmp_table);
        end
        data_table_across_phase.(phase).data_merged = data_merged;
        data_table_across_phase.(phase).table_merged = table_merged;
    end

    % plot histogram of component ammplitude
    pd_alpha = 0.4;
    amp_thres = dictionary;
    for phase = phases
        data_merged_amp = data_table_across_phase.(phase).data_merged(:,:,1);
        for comp_i = 1:3
            if comp_i ~= 2
                tmp_y = data_merged_amp(comp_i,:);
            else
                tmp_y = -data_merged_amp(comp_i,:);
            end
            pd = fitdist(tmp_y',"InverseGaussian");
            amp_thres(phase+comp_i) = icdf(pd,pd_alpha);
        end
    end
    
    for phase = phases
        data_merged = data_table_across_phase.(phase).data_merged;
        table_merged = data_table_across_phase.(phase).table_merged;
        data_merged_amp = data_merged(:,:,1);
        data_merged_time = data_merged(:,:,2);

        no_signal_ROI = all(isnan(data_merged_amp),1)';
        
        data_merged_binary = zeros(2,size(data_merged,2));
        data_merged_binary(1,~isnan(data_merged_time(1,:))) = 1;
        data_merged_binary(2,~isnan(data_merged_time(3,:))) = 1;

        nan_style = {'MarkerFaceColor','none','MarkerFaceAlpha','0','MarkerEdgeColor',[0.65,0.65,0.65],'MarkerEdgeAlpha','1',...
            'Marker','o','LineWidth',1.5};

        fig1 = figure(Position=[100,100,1600,900]);
        tiled1 = tiledlayout(fig1,1,2);
        sgtitle(tiled1,phase+" "+strjoin(mouse_names," "))
        axs1 = gobjects(1,2);
        for ax_i = 1:2
            axs1(ax_i) = nexttile(tiled1,ax_i);
        end

        tmp = lines(3);
        plot_colors1 = [0.65,0.65,0.65;tmp(1,:);tmp(2,:);...
            0.65,0.65,0.65;tmp(3,:);tmp(2,:)];
        pks_tag = ["short_latency","long_latency"];
        pks_text = ["short latency","long latency"];
        hold(axs1,"on")
        comp_i_bijection = [1,3];
        across_phase_data_tags = struct;
        for pk_i = 1:2 % loop thru short long latency peaks to get fibers with both large components
            data_merged_tag = data_merged_binary(pk_i,:);
            this_amp_thres = amp_thres(phase+comp_i_bijection(pk_i));
            data_merged_tag(data_merged_amp(comp_i_bijection(pk_i),:)<=this_amp_thres) = 0.5;
            across_phase_data_tags.(pks_tag(pk_i)) = data_merged_tag;
        end
        both_large_roi_bit = across_phase_data_tags.(pks_tag(1)) == 1 & across_phase_data_tags.(pks_tag(2)) == 1;

        for pk_i = 1:2 % loop thru short long latency peaks
            ax1 = axs1(pk_i);
            title(ax1,pks_text(pk_i));
            data_merged_tag = data_merged_binary(pk_i,:);
            this_amp_thres = amp_thres(phase+comp_i_bijection(pk_i));
            data_merged_tag(data_merged_amp(comp_i_bijection(pk_i),:)<=this_amp_thres) = 0.5;
            data_merged_tag_1 = data_merged_tag;
            data_merged_tag_1(both_large_roi_bit) = 2;

            tag_count1 = nan(1,3);
            lgs1 = gobjects(1,3);
            lgs_text1 = ["small or no ","large ","both "] + pks_text(pk_i) + " peak";

            this_bit = data_merged_tag==0 | data_merged_tag == 0.5; % cross is the same for both fig
            nan_count = sum(this_bit);

            if ~all(this_bit==0)
                p_nan_1 = circle_UI_function.scatter_3d(ax1,nan(sum(this_bit),1),table_merged(this_bit,:),nan_style=nan_style,skipBubble=no_signal_ROI(this_bit),setuserdata=0);
            end
            
            % plot fig1
            %%%%%%
            % unique_tag = [0.5,1,2];
            %%%%%%
            unique_tag = [0.6,1,2]; % 0.6 doesn't exist, just use as a placeholder for formatting, small peaks (0.5) are merged into no peaks (nan) above
            %%%%%%
            lgs_mask1 = true(1,length(unique_tag));
            for ti = 1:length(unique_tag)
                t = unique_tag(ti);
                this_bit = data_merged_tag_1==t;
                tag_count1(ti) = sum(this_bit);
                if all(this_bit==0)
                    lgs_mask1(ti) = 0;
                    continue
                else
                    p = circle_UI_function.scatter_3d(ax1,ones(sum(this_bit),1),table_merged(this_bit,:),colormapOption=plot_colors1(length(unique_tag)*(pk_i-1)+ti,:),View=[0,90],skipBubble=no_signal_ROI(this_bit),setuserdata=0);
                    lgs1(ti) = p(3);
                end
            end
            tag_count1(1) = nan_count;
            
            % account for the change that grey are merged into nan
            lgs_mask1(1) = 1;
            lgs1(1) = p_nan_1(1);

            tag_frac1 = tag_count1/sum(tag_count1);
            for i = 1:3
                lgs_text1(i) = lgs_text1(i)+"  ("+tag_count1(i)+", "+sprintf("%0.2f%%",tag_frac1(i)*100)+")";
            end

            legend(lgs1(lgs_mask1),lgs_text1(lgs_mask1));
            colorbar(ax1,"off")
            axis(ax1,"vis3d")
        end
        hold(axs1,"off")
        saveas(fig1,[cf,'sup_fig6B_hor.png'])
        for ax = axs1
            ax.View = [-90,0];
            axis(ax,"vis3d")
        end
        saveas(fig1,[cf,'sup_fig6B_sag.png'])
        delete(fig1)
    end
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 6C
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_itilick",{["late","early"]},mouse_names,'sup_fig6C',cmap_limit={[-1,1]*0.018,[-1,1]*0.007,[-1,1]*0.007});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 6DE
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G25","G26","G27"];
common_plot_functions.smoothed_spatial_map(cf,"tas_itilick",{["ITIlate","ITIearly"]},mouse_names,'sup_fig6D_smoothed',tabular_sheet_name="sup_fig6D",excel_only=true);
common_plot_functions.phase_diff_circlemap(cf,"tas_itilick",{["late","early"]},mouse_names,'sup_fig6E',cmap_limit={[-1,1]*0.018,[-1,1]*0.007,[-1,1]*0.007});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 6FG
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.modality_between_context(cf,"extinction_vs_lick",mouse_names,'sup_fig6FG')


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% functions
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
