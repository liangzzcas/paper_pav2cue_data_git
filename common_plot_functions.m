 classdef common_plot_functions
    methods(Static)

        function component_time_histogram(cf,plot_phase,save_name)
            % plot_phases = ...
            %     ["unpred","cue1late","cue1early","cue1LEDomi","cue1Toneomi",...
            %     "cue2late","cue2early","cue2LEDomi","cue2Toneomi"];
            plot_phases = plot_phase;
            if all(contains(plot_phase,["cue","unpred"]))
                cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
                tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
            elseif all(contains(plot_phase,"rew"))
                cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump_rew_merged.mat']);
                tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat']);
            else
                error("Check plot_phase, it should contains cue/unpred or rew.")
            end

            CT_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
            comp_window_cousump = common_functions.get_comp_timewindow();
            
            null_factor = 3;
            
            mouse_name_all = ["G15", "G12", "G17", "G19", "G22" ,"G21", "G23", "G24"];
            
            output = struct;
            for phase_name = plot_phase
                this_data = cwa.(phase_name);
                mouse_names = intersect(string(fields(this_data)'),mouse_name_all);
                for m_i = 1:length(mouse_names)
                    mouse_name = mouse_names(m_i);
                    n_rois = size(this_data.(mouse_name).mu,2)-3;
                    fr = size(this_data.(mouse_name).mu,1)/4;
            
                    this_mu = this_data.(mouse_name).mu(:,4:end);
                    if contains(phase_name,["unpred","rew"])
                        null_window = [1:round(0.6*fr)];
                    else
                        null_window = [1:round(1*fr)];
                    end
                    null = std(this_mu(null_window,:),[],1);
                    
                    if ~contains(phase_name,["unpred","rew"])
                        dip_window = [round(1.1*fr):round(2.5*fr)];
                        peak_window = [fr:round(2.5*fr)];
                        re_window = [fr:round(3*fr)];
                    else
                        dip_window = [round(0.8*fr):round(2.5*fr)];
                        peak_window = [round(0.8*fr):round(2.5*fr)];
                        re_window = [round(0.8*fr):round(3*fr)];
                    end
                    output_tmp = nan(3,n_rois,3);
                    for r = 1:n_rois
                        % get dip
                        [dp_pks,dp_locs] = findpeaks(-this_mu(dip_window,r),NPeaks=1,SortStr="descend",MinPeakHeight=null_factor*null(r));
                        dp_locs = dp_locs+dip_window(1)-1;
                        if ~isempty(dp_pks)
                            output_tmp(2,r,1) = -dp_pks;
                            output_tmp(2,r,2) = dp_locs/fr;
                            output_tmp(2,r,3) = dp_locs;
                        end
                        % get peak 1 & peak 2
                        if isempty(dp_pks)
                            [pk_pks,pk_locs] = findpeaks(this_mu(peak_window,r),NPeaks=1,SortStr="descend",MinPeakHeight=null_factor*null(r));
                            re_pks = []; re_locs = [];
                        else
                            if dp_locs-peak_window(1)<=1
                                pk_pks = []; pk_locs=[];
                            else
                                [pk_pks,pk_locs] = findpeaks(this_mu(peak_window(1):dp_locs,r),NPeaks=1,SortStr="descend",MinPeakHeight=null_factor*null(r));
                            end
            
                            if re_window(end)-dp_locs<=1
                                re_pks = []; re_locs=[];
                            else
                                [re_pks,re_locs] = findpeaks(this_mu(dp_locs:re_window(end),r),NPeaks=1,SortStr="descend",MinPeakHeight=null_factor*null(r));
                            end
                        end
                        pk_locs = pk_locs+peak_window(1)-1;
                        re_locs = re_locs+dp_locs-1;
                        if ~isempty(pk_pks)
                            output_tmp(1,r,1) = pk_pks;
                            output_tmp(1,r,2) = pk_locs/fr;
                            output_tmp(1,r,3) = pk_locs;
                        end
                        if ~isempty(re_pks)
                            output_tmp(3,r,1) = re_pks;
                            output_tmp(3,r,2) = re_locs/fr;
                            output_tmp(3,r,3) = re_locs;
                        end
                    end
                    output.(phase_name).(mouse_name) = output_tmp;
                end
            end
            
            % get consistent component window for early & late
            if all(contains(plot_phase,"cue")) && ~any(contains(plot_phase,"omi"))
                total_window = struct;
                pps = ["cue1","cue2"];
                for pi = 1:length(pps)
                    ps = pps(pi);
                    if ~isfield(output,ps+"early")
                        continue
                    end
                    [~,tmp_early] = merge_across_mice(output,tas,CT_table,ps+"early");
                    [~,tmp_late] = merge_across_mice(output,tas,CT_table,ps+"late");
                    tmp = cat(2,tmp_early,tmp_late);
                    total_window.(ps+"el") = [min(tmp(1,:,2),[],"all","omitmissing"),max(tmp(1,:,2),[],"all","omitmissing")]-1;
                end
            elseif all(contains(plot_phase,"cue")) && any(contains(plot_phase,"omi"))
                total_window = struct;
                pps = ["cue1","cue2"];
                for pi = 1:length(pps)
                    ps = pps(pi);
                    if ~isfield(output,ps+"LEDomi")
                        continue
                    end
                    [~,tmp_early] = merge_across_mice(output,tas,CT_table,ps+"LEDomi");
                    [~,tmp_late] = merge_across_mice(output,tas,CT_table,ps+"Toneomi");
                    tmp = cat(2,tmp_early,tmp_late);
                    total_window.(ps+"omi") = [min(tmp(1,:,2),[],"all","omitmissing"),max(tmp(1,:,2),[],"all","omitmissing")]-1;
                end
            end
            
            
            for p_i = 1:length(plot_phases)
                phase_name = plot_phases(p_i);
                mouse_names = string(fields(output.(phase_name))');
                [data_merged,~] = merge_across_mice(output,tas,CT_table,phase_name);
                
                if contains(phase_name,"unpred")
                    comp_window_phase_name = "unpred";
                else
                    comp_window_phase_name = phase_name;
                end
                input_window = comp_window_cousump.(comp_window_phase_name)-1;
                if all(contains(plot_phase,"cue"))  && ~any(contains(plot_phase,"omi"))
                    if contains(phase_name,"cue1") && ~contains(phase_name,"omi")
                        input_window(1,2) = total_window.cue1el(2);
                    elseif contains(phase_name,"cue2") && ~contains(phase_name,"omi")
                        input_window(1,2) = total_window.cue2el(2);
                    elseif contains(phase_name,"cue1") && contains(phase_name,"omi")
                        input_window(1,2) = total_window.cue1omi(2);
                    elseif contains(phase_name,"cue2") && contains(phase_name,"omi")
                        input_window(1,2) = total_window.cue2omi(2);
                    else
                        % skip
                    end
                end
            
                % set consistent ylim for early and late
                phase_ylim = [];
                if contains(phase_name,["cue1late","cue1early"])
                    phase_ylim = [0,90];
                elseif contains(phase_name,["cue1LEDomi","cue1Toneomi"])
                    phase_ylim = [0,120];
                elseif contains(phase_name,["cue2late","cue2early"])
                    phase_ylim = [0,90];
                elseif contains(phase_name,["cue2LEDomi","cue2Toneomi"])
                    phase_ylim = [0,80];
                end
            
                if contains(phase_name,["unpred","rew"])
                    phase_switch_tag = "rew";
                else
                    phase_switch_tag = "cue";
                end
            
                switch phase_switch_tag
                    case "cue"
                        phase_xlim = [0,1.5];
                        binwidth = 0.05;
                    case "rew"
                        phase_xlim = [-0.5,1.5];
                        binwidth = 0.05;
                    otherwise
                end
            
                % make shade cover a whole bin
                input_window = input_window/binwidth;
                input_window = round(input_window) .* binwidth;
                
                cum_x = phase_xlim(1):binwidth:phase_xlim(2);
                if cum_x(end) ~= phase_xlim(2)
                    cum_x(end+1) = phase_xlim(2);
                end
                
                plot_colors = lines(5);
                title_texts = ["peak","dip"];
                
                fig = figure(Position=[100,100,1600,900]);
                tiled = tiledlayout(2,1,Parent=fig,TileSpacing="tight");
                sgtitle(tiled,"components time across ROIs "+phase_name)
                axs = gobjects(1,2);
                for i = 1:2
                    axs(i) = nexttile(tiled,i);
                end
                lgs = gobjects(1,6);
                hold(axs,"on")
                data_to_plot = data_merged;
                data_to_plot(:,:,2) = data_to_plot(:,:,2)-1;
                tmp = cat(1,data_to_plot(1,:,2),data_to_plot(3,:,2));
                tmp_valid = (tmp>input_window(1,1) & tmp<input_window(1,2)) | (tmp>input_window(3,1) & tmp<input_window(3,2));
                tmp(~tmp_valid) = nan;
                h1 = histogram(axs(1),tmp,BinWidth=binwidth,FaceColor=plot_colors(1,:),Normalization="count");
                tmp1 = data_to_plot(1,:,2);
                tmp1=tmp1(tmp1>input_window(1,1) & tmp1<input_window(1,2));
                mid1 = mean(tmp1,"all","omitmissing");
                tmp1 = data_to_plot(3,:,2);
                tmp1=tmp1(tmp1>input_window(3,1) & tmp1<input_window(3,2));
                mid3 = mean(tmp1,"all","omitmissing");
                
                tmp = data_to_plot(2,:,2);
                tmp_valid = tmp>input_window(2,1) & tmp<input_window(2,2);
                tmp(~tmp_valid) = nan;
                h2 = histogram(axs(2),tmp,BinWidth=binwidth,FaceColor=plot_colors(2,:),Normalization="count");
                mid2 = mean(tmp,"all","omitmissing");
            
                if ~isempty(input_window)
                    tmp = xline(axs(1),input_window(1,:),Color=plot_colors(3,:),LineStyle='-',Linewidth=1.8);
                    lgs(2) = tmp(1);
                    tmp = xline(axs(1),input_window(3,:),Color=plot_colors(5,:),LineStyle='-',Linewidth=1.8);
                    lgs(4) = tmp(1);
                    tmp = xline(axs(2),input_window(2,:),'k-',Linewidth=1.8);
                    lgs(6) = tmp(1);
                end
                xline(axs(1),[mid1,mid3],'--')
                xline(axs(2),[mid2],'--')
                ylabel(axs(1),"component count")
                ylabel(axs(2),"component count")
                legend(axs(1),lgs([2,4]),["pk1 window","pk2 window"])
                legend(axs(2),lgs(6),["dp window"])
            
                title(axs(1),title_texts(1))
                title(axs(2),title_texts(2))
                
                xlim(axs,phase_xlim)
                if ~isempty(phase_ylim)
                    ylim(axs,phase_ylim);
                end
                linkaxes(axs,"xy")
                hold(axs,"off")
                saveas(fig,[cf,char(save_name),'_',char(phase_name),'.png'])
                delete(fig)
            end

            % helper function
            function [data_merged,data_merged_tas] = merge_across_mice(output,tas,CT_table,phase_name)
                data_merged = [];
                data_merged_tas = [];
                % if ~isfield(output,phase_name)
                %     data_merged = cat(2,data_merged,output.(phase_name).(mouse_name)(:,rois_to_plot,:));
                %     data_merged_tas = cat(2,data_merged_tas,tas.(phase_name).(mouse_name).mu_location_value(:,rois_to_plot,:));
                % end

                mouse_names = string(fields(output.(phase_name))');
                for mouse_name = mouse_names
                    rois_to_plot = logical(CT_table{string(CT_table{:,"mouse_name"})==mouse_name,"significance"}');
                    data_merged = cat(2,data_merged,output.(phase_name).(mouse_name)(:,rois_to_plot,:));
                    data_merged_tas = cat(2,data_merged_tas,tas.(phase_name).(mouse_name).mu_location_value(:,rois_to_plot,:));
                end
            end

        end
        
        function component_presence_circle(cf,plot_tag,plot_phases,save_name)
            save_name = char(save_name);
            CT_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
            switch plot_tag
                case "pav2cue task"
                    if all(contains(plot_phases,["cue","unpred"]))
                        tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                    elseif all(contains(plot_phases,"rew"))
                        tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat']);
                    else
                        error("Check plot_phase, it should contains cue/unpred or rew.");
                    end
                    % plot_phases = ["unpred","cue1early","cue1late","cue2early","cue2late","cue1LEDomi","cue2Toneomi"];
                    mouse_name_cell = {["G15","G12","G17","G19","G22","G21","G23","G24"]};
                    
                case "pav2cue ITI lick"
                    tas = load([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat']);
                    % plot_phases = ["early","late"];
                    mouse_name_cell = {["G15","G12","G17","G19","G22","G21","G23","G24"]};
            end

            for m_i = 1:size(mouse_name_cell,2)
                mouse_names = mouse_name_cell{m_i};
                for phase = plot_phases
                    data_merged = [];
                    table_merged = [];
                    for mouse_name = mouse_names
                        tmp_data = tas.(phase).(mouse_name).mu_location_value(:,:,1);
                        tmp_sig = tas.(phase).(mouse_name).mu_location_value_significance;
                        tmp_data(~logical(tmp_sig)) = nan;
                        tmp_table = CT_table(string(CT_table{:,"mouse_name"})==mouse_name,:);
                        tmp_sig = tmp_table{:,"significance"}==1;
                        tmp_data = tmp_data(:,tmp_sig);
                        tmp_table = tmp_table(tmp_sig,:);
                        data_merged = cat(2,data_merged,tmp_data);
                        table_merged = cat(1,table_merged,tmp_table);
                    end
                    % find amplitude threshold
                    pd_alpha = 0.4;
                    amp_thres = nan(1,3);
                    comp_text = ["pk","dp","re"];
                    for comp_i = 1:3
                        % skip if <= 4 data points
                        if sum(~isnan(data_merged(comp_i,:)))<=4
                            amp_thres(comp_i) = 0;
                            continue
                        end
                        if comp_i ~= 2
                            tmp_y = data_merged(comp_i,:);
                        else
                            tmp_y = -data_merged(comp_i,:);
                        end
                        pd = fitdist(tmp_y',"InverseGaussian");
                        this_thres = icdf(pd,pd_alpha);
                        amp_thres(comp_i) = this_thres;
                    end
                    
                    % apply amplitude threshold
                    for comp_i = 1:3
                        data_merged(comp_i,abs(data_merged(comp_i,:))<amp_thres(comp_i)) = nan;
                    end
                    data_merged(:,table_merged{:,"fiber_bottom_DV"}<2) = nan;
            
                    % combine data
                    data_merged_binary = ~isnan(data_merged);
                    data_merged_tag = [1,3,5]*data_merged_binary;
                    unique_tag = [1,3,4,8,9,6,5];
                    tag_count = nan(size(unique_tag));
                    
                    plot_colors = lines(12);
                    plot_colors(4,:) = plot_colors(6,:);
                    plot_colors([6,7],:) = [0.8,0.8,0.8;0.8,0.8,0.8;];
                    lgs = gobjects(1,length(unique_tag));
                    fig = figure(Position=[100,100,1600,900]);
                    ax = axes(fig);
                    title(ax,phase+" "+strjoin(mouse_names," "))
                    hold(ax,"on")
                    this_bit = data_merged_tag==0;
            
                    nan_count = sum(this_bit);
                    if ~all(this_bit==0) % this is always true...
                        nan_lgs = common_functions.scatter_3d(ax,nan(sum(this_bit),1),table_merged(this_bit,:),setuserdata=0,skip_nan=0,...
                            nan_style = {'MarkerFaceColor','none','MarkerFaceAlpha','0','MarkerEdgeColor','[.8,.8,.8]','MarkerEdgeAlpha','1',...
                            'Marker','o','LineWidth',1});
                        nan_lgs = nan_lgs(1);
                    end
                    lgs_mask = true(1,length(unique_tag));
                    for ti = 1:length(unique_tag)
                        t = unique_tag(ti);
                        this_bit = data_merged_tag==t;
                        tag_count(ti) = sum(this_bit);
                        if all(this_bit==0)
                            lgs_mask(ti) = 0;
                            continue
                        else
                            p = common_functions.scatter_3d(ax,ones(sum(this_bit),1),table_merged(this_bit,:),colormapOption=plot_colors(ti,:),View=[0,90],setuserdata=0,skip_nan=1);
                            lgs(ti) = p(3);
                        end
                    end
                    hold(ax,"off")
                    lgs_text = ["pk","dp","pk+dp","dp+re","pk+dp+re","pk+re","re"];
                    lgs_mask([6]) = any(lgs_mask([6,7]));
                    if lgs_mask([7])
                        lgs([6]) = lgs([7]);
                    end
                    lgs_mask([7]) = false;
                    lgs_text(6) = "pk+re";
                    tag_count(6) = tag_count(6)+tag_count(7);
            
                    lgs_text = lgs_text(lgs_mask);
                    tag_count = tag_count(lgs_mask);
            
                    tag_count = [tag_count,nan_count]; % include empty
                    lgs_text = [lgs_text,"empty"]; % include empty
            
                    tag_frac = tag_count/sum(tag_count);
                    for i = 1:length(lgs_text)
                        lgs_text(i) = lgs_text(i)+"  ("+tag_count(i)+", "+sprintf("%0.2f%%",tag_frac(i)*100)+")";
                    end
            
                    legend([lgs(lgs_mask),nan_lgs],lgs_text);
                    colorbar(ax,"off")
                    axis(ax,"vis3d")
                    saveas(fig,[cf,save_name,'_',char(phase),'_hor.png'])
                    ax.View = [-90,0];
                    axis(ax,"vis3d")
                    saveas(fig,[cf,save_name,'_',char(phase),'_sag.png'])
                    delete(fig)
                end
            end
        end
        
        function [this_value_out,this_ranksum_p_out] = phase_diff_data(cf,plot_tag,plot_phases,mouse_names,varargin)
            % plot_tag = ["tas_normal","tas_rewmerged","tas_itilick"];
            % mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
            % phase_names = {["cue1late","cue1early"],["cue2late","cue2early"],["cue1LEDomi","cue1late"],...
            %     ["rew1late","rew1early"],["ITIlate","ITIearly"]};
            ranksum_alpha = 0.01;
            ip = inputParser;
            ip.addParameter("ranksum_alpha",0.01)
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end            
            
            ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
            switch plot_tag
                case "tas_normal"
                    this_tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                case "tas_rewmerged"
                    this_tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat']);
                case "tas_itilick"
                    this_tas = load([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat']);
                    if ~contains(string(varargin(1:2:end)),"ranksum_alpha")
                        ranksum_alpha = 0.025;
                    end
            end
            
            
            comp_text = ["pk","dp","re"];
            
            across_data = struct; % store diff and sig for both cue in this struct
            for p_i = 1:length(plot_phases)
                phase_name = plot_phases{p_i};
                phase_name_tag = strjoin(phase_name,"_");
                phase_name_text = strjoin(phase_name,"-");
                mouse_name_tag = strjoin(mouse_names,"_");
            
                this_table = [];
                mu_1 = [];
                single_1 = {};
                sig_1 = [];
                mu_2 = [];
                single_2 = {};
                sig_2 = [];
            
                for m = mouse_names
                    this_table = cat(1,this_table,ct_table(ct_table{:,"mouse_name"}==m,:));
            
                    tmp = getfield(this_tas,phase_name(1),m,"mu_location_value");
                    mu_1 = cat(2,mu_1,tmp(:,:,1));
                    single_1 = cat(2,single_1,num2cell(getfield(this_tas,phase_name(1),m,"single_values"),3));
                    sig_1 = cat(2,sig_1,getfield(this_tas,phase_name(1),m,"mu_location_value_significance"));
            
                    tmp = getfield(this_tas,phase_name(2),m,"mu_location_value");
                    mu_2 = cat(2,mu_2,tmp(:,:,1));
                    single_2 = cat(2,single_2,num2cell(getfield(this_tas,phase_name(2),m,"single_values"),3));
                    sig_2 = cat(2,sig_2,getfield(this_tas,phase_name(2),m,"mu_location_value_significance"));
                end
                in_striatum_bit = logical(this_table{:,"significance"});
                in_striatum_rois = find(in_striatum_bit);
            
            
                % get diff
                tmp = diff_phase_ranksum(mu_1,sig_1,single_1,mu_2,sig_2,single_2,in_striatum_bit);
                across_data.(mouse_name_tag).(phase_name_tag) = tmp;
            
                stats_data.value = nan([3,sum(in_striatum_bit)]);
                stats_data.type = nan([3,sum(in_striatum_bit)]);
            
                this_value_1 = across_data.(mouse_name_tag).(phase_name_tag).diff_mu';
                this_ranksum_p_1 = across_data.(mouse_name_tag).(phase_name_tag).ranksum_p';
                this_value_out = this_value_1(:,in_striatum_bit);
                this_ranksum_p_out = this_ranksum_p_1(:,in_striatum_bit);

                for comp_i = 1:3            
                    tmp = this_value_1(comp_i,in_striatum_bit);
                    stats_data.value(comp_i,:) = tmp;
                    stats_data.type(comp_i,tmp>0) = 1;
                    stats_data.type(comp_i,tmp<0) = -1;
                    tmp = this_ranksum_p_1(comp_i,in_striatum_bit)>ranksum_alpha;
                    stats_data.type(comp_i,tmp) = 0;
                end
            end
            
            % helper function
            function diff_struct = diff_phase_ranksum(mu1,sig1,single1,mu2,sig2,single2,in_striatum_bit)
                mu1(~sig1) = nan; mu2(~sig2) = nan;
                single1(~sig1) = cellfun(@(x) nan(size(x)),single1(~sig1),UniformOutput=false);
                single2(~sig2) = cellfun(@(x) nan(size(x)),single2(~sig2),UniformOutput=false);
                
                this_value = common_functions.nan_minus(mu1,mu2);
                this_ranksum_p = common_functions.nan_ranksum_cell(single1,single2);
                diff_struct.included_bit = in_striatum_bit;
                diff_struct.diff_mu = this_value';
                diff_struct.ranksum_p = this_ranksum_p';
            end

        end
        

        function phase_diff_circlemap(cf,plot_tag,plot_phases,mouse_names,save_name,varargin)
            % plot_tag = ["tas_normal","tas_rewmerged","tas_itilick"];
            % mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
            % phase_names = {["cue1late","cue1early"],["cue2late","cue2early"],["cue1LEDomi","cue1late"],...
            %     ["rew1late","rew1early"],["ITIlate","ITIearly"]};
            save_name = char(save_name);
            cmap_limit = {[],[],[]};

            ranksum_alpha = 0.01;
            ip = inputParser;
            ip.addParameter("ranksum_alpha",0.01)
            ip.addParameter("cmap_limit",{[],[],[]})
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end            
            
            ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
            switch plot_tag
                case "tas_normal"
                    this_tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                case "tas_rewmerged"
                    this_tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat']);
                case "tas_itilick"
                    this_tas = load([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat']);
                    if ~contains(string(varargin(1:2:end)),"ranksum_alpha")
                        ranksum_alpha = 0.025;
                    end
            end
            
            
            comp_text = ["pk","dp","re"];
            
            across_data = struct; % store diff and sig for both cue in this struct
            for p_i = 1:length(plot_phases)
                phase_name = plot_phases{p_i};
                phase_name_tag = strjoin(phase_name,"_");
                phase_name_text = strjoin(phase_name,"-");
                mouse_name_tag = strjoin(mouse_names,"_");
            
                this_table = [];
                mu_1 = [];
                single_1 = {};
                sig_1 = [];
                mu_2 = [];
                single_2 = {};
                sig_2 = [];
            
                for m = mouse_names
                    this_table = cat(1,this_table,ct_table(ct_table{:,"mouse_name"}==m,:));
            
                    tmp = getfield(this_tas,phase_name(1),m,"mu_location_value");
                    mu_1 = cat(2,mu_1,tmp(:,:,1));
                    single_1 = cat(2,single_1,num2cell(getfield(this_tas,phase_name(1),m,"single_values"),3));
                    sig_1 = cat(2,sig_1,getfield(this_tas,phase_name(1),m,"mu_location_value_significance"));
            
                    tmp = getfield(this_tas,phase_name(2),m,"mu_location_value");
                    mu_2 = cat(2,mu_2,tmp(:,:,1));
                    single_2 = cat(2,single_2,num2cell(getfield(this_tas,phase_name(2),m,"single_values"),3));
                    sig_2 = cat(2,sig_2,getfield(this_tas,phase_name(2),m,"mu_location_value_significance"));
                end
                in_striatum_bit = logical(this_table{:,"significance"});
                in_striatum_rois = find(in_striatum_bit);
            
            
                % get diff
                tmp = diff_phase_ranksum(mu_1,sig_1,single_1,mu_2,sig_2,single_2,in_striatum_bit);
                across_data.(mouse_name_tag).(phase_name_tag) = tmp;
            
                stats_data.value = nan([3,sum(in_striatum_bit)]);
                stats_data.type = nan([3,sum(in_striatum_bit)]);
            
                % circle plot
                fig_raw = figure(Position=[100,100,1600,900]);
                tiled_raw = tiledlayout(fig_raw,2,3,TileSpacing="tight");
                sgtitle(tiled_raw,phase_name_text+" "+strjoin(mouse_names," "));
                axs_raw = gobjects(2,3);
            
                for i = 1:6
                    axs_raw(i) = nexttile(tiled_raw,i);
                end
                hold(axs_raw,"on")
            
                this_value_1 = across_data.(mouse_name_tag).(phase_name_tag).diff_mu';
                this_ranksum_p_1 = across_data.(mouse_name_tag).(phase_name_tag).ranksum_p';
                
            
                for comp_i = 1:3
                    if ~isempty(cmap_limit{comp_i})
                        this_cmap_lim = cmap_limit{comp_i};
                    else
                        this_cmap_lim = [];
                    end
                    common_functions.scatter_3d_with_datatip(axs_raw(comp_i),this_value_1(comp_i,:),this_table,...
                        grey_bubble=this_ranksum_p_1(comp_i,:)>ranksum_alpha,dynamic_cmap_bound_flag=1,grey_bubble_pvalue=this_ranksum_p_1(comp_i,:),...
                        cmapBounds=this_cmap_lim,viewAngle=[0,90],colormapOption="redblue",skipBubble=(~in_striatum_bit)',sorting_tag='abs');
            
                    common_functions.scatter_3d_with_datatip(axs_raw(comp_i+3),this_value_1(comp_i,:),this_table,...
                        grey_bubble=this_ranksum_p_1(comp_i,:)>ranksum_alpha,dynamic_cmap_bound_flag=1,grey_bubble_pvalue=this_ranksum_p_1(comp_i,:),...
                         cmapBounds=this_cmap_lim,viewAngle=[-90,0],colormapOption="redblue",skipBubble=(~in_striatum_bit)',sorting_tag='abs');
            
                    title(axs_raw(comp_i),comp_text(comp_i))
            
                    tmp = this_value_1(comp_i,in_striatum_bit);
                    stats_data.value(comp_i,:) = tmp;
                    stats_data.type(comp_i,tmp>0) = 1;
                    stats_data.type(comp_i,tmp<0) = -1;
                    tmp = this_ranksum_p_1(comp_i,in_striatum_bit)>ranksum_alpha;
                    stats_data.type(comp_i,tmp) = 0;
                end
            
                hold(axs_raw,"off")
                saveas(fig_raw,[cf,save_name,'_',char(phase_name_tag),'.png'])
                delete(fig_raw)
                
                % % plot stats figure
                % fig_raw = figure(Position=[100,100,1600,900]);
                % tiled_raw = tiledlayout(fig_raw,2,3,TileSpacing="tight");
                % sgtitle(tiled_raw,phase_name_text+" "+strjoin(mouse_names," "));
                % axs_raw = gobjects(2,3);
                % for i = 1:6
                %     axs_raw(i) = nexttile(tiled_raw,i);
                % end
                % hold(axs_raw,"on")
                % for comp_i=1:3
                %     roi_types = [sum(stats_data.type(comp_i,:) == 1),sum(stats_data.type(comp_i,:) == -1),sum(stats_data.type(comp_i,:) == 0),sum(isnan(stats_data.type(comp_i,:)))];
                %     tmp = length(stats_data.type(comp_i,:));
                %     roi_type_text = sprintf("red/blue/grey/x = %d/%d/%d/%d out of %d, %0.2f/%0.2f/%0.2f/%0.2f%%",roi_types(1),roi_types(2),roi_types(3),roi_types(4),tmp,...
                %         roi_types(1)/tmp*100,roi_types(2)/tmp*100,roi_types(3)/tmp*100,roi_types(4)/tmp*100);
                %     axs = axs_raw([comp_i,comp_i+3]);
                % 
                %     histogram(axs(1),stats_data.value(comp_i,:),BinWidth=0.0005,Normalization="probability")
                %     tmp = stats_data.value(comp_i,:);
                %     tmp(stats_data.type(comp_i,:)==0) = nan;
                %     histogram(axs(2),tmp,BinWidth=0.0005,Normalization="probability")
                %     title(axs(1),["including grey",roi_type_text]);title(axs(2),["excluding grey",roi_type_text]);
                %     xlabel(axs(1),"\DeltaF/F");xlabel(axs(2),"\DeltaF/F");
                % end
                % hold(axs_raw,"off")
                % saveas(fig_raw,[save_dir,'stats_',char(phase_name_tag),'_',char(mouse_name_tag),'.fig'])
                % saveas(fig_raw,[save_dir,'stats_',char(phase_name_tag),'_',char(mouse_name_tag),'.png'])
                % close(fig_raw)
            end
            
            % helper function
            function diff_struct = diff_phase_ranksum(mu1,sig1,single1,mu2,sig2,single2,in_striatum_bit)
                mu1(~sig1) = nan; mu2(~sig2) = nan;
                single1(~sig1) = cellfun(@(x) nan(size(x)),single1(~sig1),UniformOutput=false);
                single2(~sig2) = cellfun(@(x) nan(size(x)),single2(~sig2),UniformOutput=false);
                
                this_value = common_functions.nan_minus(mu1,mu2);
                this_ranksum_p = common_functions.nan_ranksum_cell(single1,single2);
                diff_struct.included_bit = in_striatum_bit;
                diff_struct.diff_mu = this_value';
                diff_struct.ranksum_p = this_ranksum_p';
            end

        end
        
        function smoothed_spatial_map(cf,plot_tag,plot_phases,mouse_names,save_name,varargin)
            % plot_tag = ["tas_normal","tas_rewmerged","tas_itilick"];
            % mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
            % phase_names = {["cue1late","cue1early"],["cue2late","cue2early"],["cue1LEDomi","cue1late"],...
            %     ["rew1late","rew1early"],["ITIlate","ITIearly"]};
            save_name = char(save_name);
            cmap_limit = {[],[],[]};
            tabular_sheet_name = '';
            excel_only = false;
            
            ranksum_alpha = 0.01;
            ip = inputParser;
            ip.addParameter("ranksum_alpha",0.01)
            ip.addParameter("cmap_limit",{[],[],[]})
            ip.addParameter("tabular_sheet_name",'')
            ip.addParameter("excel_only",false)
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            if isempty(tabular_sheet_name)
                tabular_sheet_name = save_name;
            end

            hotzone_prctile_thres=90;
            scaleUnit = 0.05;
            hotzone_color_posneg = {'y','c'};

            ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
            switch plot_tag
                case "tas_normal"
                    this_tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                case "tas_rewmerged"
                    this_tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat']);
                case "tas_itilick"
                    this_tas = load([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat']);
                    if ~contains(string(varargin(1:2:end)),"ranksum_alpha")
                        ranksum_alpha = 0.025;
                    end
            end
            comp_text = ["pk","dp","re"];

            %%
            across_data = struct; % store diff and sig for both cue in this struct
            for p_i = 1:length(plot_phases)
                phase_name = plot_phases{p_i};
                phase_name_tag = strjoin(phase_name,"_");
                phase_name_text = strjoin(phase_name,"-");
                mouse_name_tag = strjoin(mouse_names,"_");
                if contains(phase_name(1),"ITI")
                    tmp = char(phase_name(1)); phase_name(1) = string(tmp(4:end));
                    tmp = char(phase_name(2)); phase_name(2) = string(tmp(4:end));
                end

                this_table = [];
                mu_1 = [];
                single_1 = {};
                sig_1 = [];
                mu_2 = [];
                single_2 = {};
                sig_2 = [];
            
                for m = mouse_names
                    this_table = cat(1,this_table,ct_table(ct_table{:,"mouse_name"}==m,:));
            
                    tmp = getfield(this_tas,phase_name(1),m,"mu_location_value");
                    mu_1 = cat(2,mu_1,tmp(:,:,1));
                    single_1 = cat(2,single_1,num2cell(getfield(this_tas,phase_name(1),m,"single_values"),3));
                    sig_1 = cat(2,sig_1,getfield(this_tas,phase_name(1),m,"mu_location_value_significance"));
            
                    tmp = getfield(this_tas,phase_name(2),m,"mu_location_value");
                    mu_2 = cat(2,mu_2,tmp(:,:,1));
                    single_2 = cat(2,single_2,num2cell(getfield(this_tas,phase_name(2),m,"single_values"),3));
                    sig_2 = cat(2,sig_2,getfield(this_tas,phase_name(2),m,"mu_location_value_significance"));
                end
                in_striatum_bit = logical(this_table{:,"significance"});
            
            
                % get diff
                tmp = diff_phase_ranksum(mu_1,sig_1,single_1,mu_2,sig_2,single_2,in_striatum_bit);
                across_data.(mouse_name_tag).(phase_name_tag) = tmp;
            
                stats_data.value = nan([3,sum(in_striatum_bit)]);
                stats_data.type = nan([3,sum(in_striatum_bit)]);

                % get fiber stats
                this_value_1 = across_data.(mouse_name_tag).(phase_name_tag).diff_mu';
                this_ranksum_p_1 = across_data.(mouse_name_tag).(phase_name_tag).ranksum_p';
                for comp_i = 1:3
                    tmp = this_value_1(comp_i,in_striatum_bit);
                    stats_data.value(comp_i,:) = tmp;
                    stats_data.type(comp_i,tmp>0) = 1;
                    stats_data.type(comp_i,tmp<0) = -1;
                    tmp = this_ranksum_p_1(comp_i,in_striatum_bit)>ranksum_alpha;
                    stats_data.type(comp_i,tmp) = 0;
                end
                in_striatum_ct = this_table(in_striatum_bit,:);
                frac_data = [];
                frac_data_mouses = table;
                title_texts = [];
                for comp_i=1:3
                    roi_types = [sum(stats_data.type(comp_i,:) == 1),sum(stats_data.type(comp_i,:) == -1),sum(stats_data.type(comp_i,:) == 0),sum(isnan(stats_data.type(comp_i,:)))];
                    frac_data_mouse = table(Size=[length(mouse_names),4],VariableType=["double","double","double","double"],RowNames=mouse_names,VariableNames=comp_text(comp_i)+" "+["increase/red","decrease/blue","nochange/grey","x"]);
                    for mi = 1:length(mouse_names)
                        tmp = stats_data.type(comp_i,in_striatum_ct.mouse_name==mouse_names(mi));
                        frac_data_mouse(mi,:) = {sum(tmp==1),sum(tmp==-1),sum(tmp==0),sum(isnan(tmp))};
                    end
                    tmp = length(stats_data.type(comp_i,:));
                    roi_type_text = sprintf("red/blue/grey/x = %d/%d/%d/%d out of %d, %0.2f/%0.2f/%0.2f/%0.2f%%",roi_types(1),roi_types(2),roi_types(3),roi_types(4),tmp,...
                        roi_types(1)/tmp*100,roi_types(2)/tmp*100,roi_types(3)/tmp*100,roi_types(4)/tmp*100);
                    frac_data = cat(1,frac_data,roi_types(1:4)/tmp);
                    title_texts = cat(1,title_texts,roi_type_text);
                    frac_data_mouses = cat(2,frac_data_mouses,frac_data_mouse);
                end
                writetable(frac_data_mouses,[cf,'tabular_data_of_figures_autogenerated.xlsx'],Sheet=tabular_sheet_name+"_"+phase_name_tag,WriteMode="overwritesheet");
                
                if excel_only
                    continue
                end

                % circle plot
                fig_grid = figure(Position=[100,100,1600,900]);
                tiled_grid = tiledlayout(fig_grid,2,3,TileSpacing="tight");
                sgtitle(tiled_grid,phase_name_text+" "+strjoin(mouse_names," "));
                axs_grid = gobjects(2,3);
            
                for i = 1:6
                    axs_grid(i) = nexttile(tiled_grid,i);
                end
                hold(axs_grid,"on");
                
                for comp_i = 1:3
                    if ~isempty(cmap_limit{comp_i})
                        this_cmap_lim = cmap_limit{comp_i};
                    else
                        this_cmap_lim = [];
                    end
                    tmp = this_value_1(comp_i,:);
                    tmp(isnan(tmp))=0;

                    [~,~,~,~] = common_plot_functions.interp_scatter_3d(axs_grid(comp_i),tmp,this_table,scaleUnit=scaleUnit,plot_value_nan=0,cmapBounds=this_cmap_lim,hotzone_prctile_thres=hotzone_prctile_thres,...
                        insig_bit=this_ranksum_p_1(comp_i,:)>ranksum_alpha,viewplane="hor",colormapOption="redblue",skipBubble=(~in_striatum_bit)');
                    [~,~,~,~] = common_plot_functions.interp_scatter_3d(axs_grid(comp_i+3),tmp,this_table,scaleUnit=scaleUnit,plot_value_nan=0,cmapBounds=this_cmap_lim,hotzone_prctile_thres=hotzone_prctile_thres,...
                        insig_bit=this_ranksum_p_1(comp_i,:)>ranksum_alpha,viewplane="sag",colormapOption="redblue",skipBubble=(~in_striatum_bit)');
                    
                    % get gradient
                    if comp_i ~= 2
                        tmp_mdl = -this_value_1(comp_i,:);
                    else
                        tmp_mdl = this_value_1(comp_i,:);
                    end
                    
                    lme_mdl = common_plot_functions.lme_3Dgradient(tmp_mdl(in_striatum_bit),this_table(in_striatum_bit,:));
                    APMLDV_coef = double(lme_mdl.Coefficients(2:4,2));
                    common_plot_functions.plot_3Dgradient(axs_grid(comp_i),APMLDV_coef,"hor");
                    common_plot_functions.plot_3Dgradient(axs_grid(comp_i+3),APMLDV_coef,"sag");
                end
            
                hold(axs_grid,"off");
                saveas(fig_grid,[cf,save_name,'_',char(phase_name_tag),'.png'])
                delete(fig_grid)
            end
            
            % helper function
            function diff_struct = diff_phase_ranksum(mu1,sig1,single1,mu2,sig2,single2,in_striatum_bit)
                mu1(~sig1) = nan; mu2(~sig2) = nan;
                single1(~sig1) = cellfun(@(x) nan(size(x)),single1(~sig1),UniformOutput=false);
                single2(~sig2) = cellfun(@(x) nan(size(x)),single2(~sig2),UniformOutput=false);
                
                this_value = common_functions.nan_minus(mu1,mu2);
                this_ranksum_p = common_functions.nan_ranksum_cell(single1,single2);
                diff_struct.included_bit = in_striatum_bit;
                diff_struct.diff_mu = this_value';
                diff_struct.ranksum_p = this_ranksum_p';
            end

        end
        
        function gradient_test(cf,plot_tag,plot_phase_comps,mouse_names,phase_tags,phase_names_tags,save_name,varargin)
            % plot_tag = ["tas_normal","tas_rewmerged","tas_itilick"];
            % mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
            % plot_phase_comps = {...
            %     {{"cue1late_cue1early",1,1;"cue1late_cue1early",2,2;},{"cue2late_cue2early",1,1;"cue2late_cue2early",2,2;},};...
            %     {{"cue1LEDomi_cue1late",1,2;"cue1LEDomi_cue1late",2,1;},{"cue1LEDomi_cue1late",1,2;"cue1LEDomi_cue1late",3,1;},{"cue1LEDomi_cue1late",2,1;"cue1LEDomi_cue1late",3,1;},};...
            %     {{"cue2Toneomi_cue2late",1,2;"cue2Toneomi_cue2late",2,1;},{"cue2Toneomi_cue2late",1,2;"cue2Toneomi_cue2late",3,1;},{"cue2Toneomi_cue2late",2,1;"cue2Toneomi_cue2late",3,1;},};...
            %     {{"cue1late_cue1early",1,1;"cue2late_cue2early",1,1;},{"cue1late_cue1early",2,2;"cue2late_cue2early",2,2;},{"cue1late_cue1early",3,1;"cue2late_cue2early",3,1;}};...
            %     {{"cue1LEDomi_cue1late",1,2;"cue2Toneomi_cue2late",1,2;},{"cue1LEDomi_cue1late",2,1;"cue2Toneomi_cue2late",2,1;},{"cue1LEDomi_cue1late",3,1;"cue2Toneomi_cue2late",3,1;},};...
            %     };
            % phase_tags = {...
            %     {{"cue1_early_pk";"cue2_early_pk";},{"cue1_dip";"cue2_dip";},};...
            %     {{"cue1_pk1";"cue1_pk2";"cue1_dp3";},{"cue1_dp1";"cue1_re2";"cue1_re3";},};...
            %     {{"cue2_pk1";"cue2_pk2";"cue2_dp3";},{"cue2_dp1";"cue2_re2";"cue2_re3";},};...
            %     {{"cue1_pk";"cue1_dip";"cue1_re";},{"cue2_pk";"cue2_dip";"cue2_re"},};...
            %     {{"cue1_pk";"cue1_dip";"cue1_re";},{"cue2_pk";"cue2_dip";"cue2_re"},};...
            %     };
            % phase_names_tags = [...
            %     "cue_learning",...
            %     "cue1_extinction",...
            %     "cue2_extinction",...
            %     "cue_learning_1vs2",...
            %     "cue_extinction_1vs2",...
            %     ];

            save_name = char(save_name);
            tabular_sheet_name = '';

            ranksum_alpha = 0.01;
            ip = inputParser;
            ip.addParameter("ranksum_alpha",0.01)
            ip.addParameter("tabular_sheet_name",'')
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            if isempty(tabular_sheet_name)
                tabular_sheet_name = save_name;
            end

            ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
            switch plot_tag
                case "tas_normal"
                    this_tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                case "tas_rewmerged"
                    this_tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat']);
                case "tas_itilick"
                    this_tas = load([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat']);
                    if ~contains(string(varargin(1:2:end)),"ranksum_alpha")
                        ranksum_alpha = 0.025;
                    end
            end
            comp_tag = ["pk","dp","re"]; posneg_tag = ["pos","neg"];

            % build data first
            across_data = struct; % store diff and sig for both cue in this struct
            tmp = cat(1,plot_phase_comps{:});
            plot_phases_combined = unique([tmp{:,1}]);
            for p_i = 1:length(plot_phases_combined)
                phase_name = strsplit(plot_phases_combined(p_i),"_");
                phase_name_tag = strjoin(phase_name,"_");
                phase_name_text = strjoin(phase_name,"-");
                mouse_name_tag = strjoin(mouse_names,"_");
                if contains(phase_name(1),"ITI")
                    tmp = char(phase_name(1)); phase_name(1) = string(tmp(4:end));
                    tmp = char(phase_name(2)); phase_name(2) = string(tmp(4:end));
                end

                this_table = [];
                mu_1 = [];
                single_1 = {};
                sig_1 = [];
                mu_2 = [];
                single_2 = {};
                sig_2 = [];
            
                for m = mouse_names
                    this_table = cat(1,this_table,ct_table(ct_table{:,"mouse_name"}==m,:));
            
                    tmp = getfield(this_tas,phase_name(1),m,"mu_location_value");
                    mu_1 = cat(2,mu_1,tmp(:,:,1));
                    single_1 = cat(2,single_1,num2cell(getfield(this_tas,phase_name(1),m,"single_values"),3));
                    sig_1 = cat(2,sig_1,getfield(this_tas,phase_name(1),m,"mu_location_value_significance"));
            
                    tmp = getfield(this_tas,phase_name(2),m,"mu_location_value");
                    mu_2 = cat(2,mu_2,tmp(:,:,1));
                    single_2 = cat(2,single_2,num2cell(getfield(this_tas,phase_name(2),m,"single_values"),3));
                    sig_2 = cat(2,sig_2,getfield(this_tas,phase_name(2),m,"mu_location_value_significance"));
                end
                in_striatum_bit = logical(this_table{:,"significance"});

                % get diff
                tmp = diff_phase_ranksum(mu_1,sig_1,single_1,mu_2,sig_2,single_2,in_striatum_bit);
                across_data.(mouse_name_tag).(phase_name_tag) = tmp;
            end

            % test gradient
            this_ct = ct_table(ct_table.significance==1 & contains(ct_table.mouse_name,mouse_names),:);
            nc = size(plot_phase_comps,2);
            for compi=1:nc
                gradient_table = table(Size=[7,3],VariableType=repmat("double",[1,3]),...
                    RowNames=[phase_tags{1}{compi}+[" mean"," sem"," p"],phase_tags{2}{compi}+[" mean"," sem"," p"],phase_tags{1}{compi}+" vs "+phase_tags{2}{compi}+" p"],...
                    VariableNames=["AP","ML","DV",]);
                this_phase = plot_phase_comps{compi};
                APMLDV_coef = [];
                APMLDV_se = [];
                APMLDV_sig = [];
                test_value = [];
                for pi = 1:size(this_phase,1)
                    phase_name_tag = this_phase{pi,1};
                    phase_name_comp = this_phase{pi,2};
                    tmp = across_data.(mouse_name_tag).(phase_name_tag).included_bit;
                    this_value_1 = across_data.(mouse_name_tag).(phase_name_tag).diff_mu(tmp,phase_name_comp)';
                    % get gradient
                    %   set specific gradient direction
                    switch phase_names_tags
                        case "cue1_extinction"
                            if phase_name_comp == 1
                                tmp_mdl = -this_value_1;
                            else
                                tmp_mdl = this_value_1;
                            end
                        case "cue2_extinction"
                            if phase_name_comp == 1
                                tmp_mdl = -this_value_1;
                            else
                                tmp_mdl = this_value_1;
                            end
                        case "cue_extinction"
                            if phase_name_comp == 1
                                tmp_mdl = -this_value_1;
                            else
                                tmp_mdl = this_value_1;
                            end
                        otherwise
                            if phase_name_comp ~= 2
                                tmp_mdl = -this_value_1;
                            else
                                tmp_mdl = this_value_1;
                            end
                    end
                    lme_mdl = common_plot_functions.lme_3Dgradient(tmp_mdl,this_ct);
                    APMLDV_coef = cat(1,APMLDV_coef,double(lme_mdl.Coefficients(2:4,2))');
                    APMLDV_se = cat(1,APMLDV_se,double(lme_mdl.Coefficients(2:4,3))');
                    APMLDV_sig = cat(1,APMLDV_sig,double(lme_mdl.Coefficients(2:4,6))');
                    test_value = cat(2,test_value,tmp_mdl');
                end
                lme_mdl = common_plot_functions.lme_3Dgradient_interact(test_value,this_ct);
                gradient_table{7,:} = double(lme_mdl.Coefficients(6:8,6))';
                gradient_table{1:6,:} = [APMLDV_coef(1,:);APMLDV_se(1,:);APMLDV_sig(1,:);APMLDV_coef(2,:);APMLDV_se(2,:);APMLDV_sig(2,:)];
                writetable(gradient_table,[cf,'tabular_data_of_figures_autogenerated.xlsx'],Sheet=strjoin([tabular_sheet_name,phase_names_tags,string(compi)],"_"),WriteMode="overwritesheet",WriteRowNames=true);
            end
            
            % helper function
            function diff_struct = diff_phase_ranksum(mu1,sig1,single1,mu2,sig2,single2,in_striatum_bit)
                mu1(~sig1) = nan; mu2(~sig2) = nan;
                single1(~sig1) = cellfun(@(x) nan(size(x)),single1(~sig1),UniformOutput=false);
                single2(~sig2) = cellfun(@(x) nan(size(x)),single2(~sig2),UniformOutput=false);
                
                this_value = common_functions.nan_minus(mu1,mu2);
                this_ranksum_p = common_functions.nan_ranksum_cell(single1,single2);
                diff_struct.included_bit = in_striatum_bit;
                diff_struct.diff_mu = this_value';
                diff_struct.ranksum_p = this_ranksum_p';
            end
        end
        
        function [this_plot,imagesc_data,interp_f,hotzone_hull,hotzone_com] = interp_scatter_3d(plot_ax,values,fib_loc_ori,varargin)
            interp_method = ["interp","3dconv"];

            skipBubble = [];
            insig_bit = [];
            flag_insig_zero = true;
            flag_include_x_as_zero = false;
            cmapBounds = [];
            colormapOption = '';
            plot_value_nan = nan;
            color_nan = [];
            
            hotzone_prctile_thres = nan;

            conv_radius = 1.2;

            viewplane = "sag";
            APlim = [];
            MLlim = [];
            DVlim = [];            
            scaleUnit = 0.1;
            
            plot_image = true;
            plot_hotzone = true;
            output_contour_path = [];

            ip = inputParser;
            ip.addParameter('interp_method',"3dconv"); % 3dconv or interp

            ip.addParameter('skipBubble',false(size(values)));
            ip.addParameter('insig_bit',false(size(values)));
            ip.addParameter('flag_insig_zero',true);
            ip.addParameter('flag_include_x_as_zero',true);

            ip.addParameter('cmapBounds',[]); % manually set the color bounds (note RedBlue will automatically make it symmetric)
            ip.addParameter('colormapOption','parula'); % colormap (either the string to a map, or a matrix)
            ip.addParameter('plot_value_nan',nan);
            ip.addParameter('color_nan',[1,1,1]*0.4);
            % method 1: interp params
            ip.addParameter('hotzone_prctile_thres',70); % 75 is too small for 3dconv method, could increase
            % method 2: 3dconv params
            ip.addParameter('conv_radius',1.2);
            ip.addParameter('exp_factor',6);
            %%% axes options
            ip.addParameter('viewplane',"sag"); % default view is sagittal
            ip.addParameter('APlim',[-2.1,1.6]); % AP limits: if blank, it will be determined automatically. I like [-2 -2].
            ip.addParameter('MLlim',[0.5,4.2]); % ML limits: if blank, it will be determined automatically. I like [-.5 3.75] for R implants
            ip.addParameter('DVlim',[2,6]); % DV limits: if blank, it will be determined automatically. I like [-5.5, -2].
            
            ip.addParameter('scaleUnit',0.05);
            % other hyperparams
            ip.addParameter('plot_image',true)
            ip.addParameter('plot_hotzone',true)
            ip.addParameter('output_contour_path',[])
            %%% parser
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            this_plot=[];

            apmldv_end = [APlim;MLlim;DVlim;];
            insig_bit = logical(insig_bit);
            if istable(fib_loc_ori)
                fib_loc_ori = fib_loc_ori{:,["fiber_bottom_AP","fiber_bottom_ML","fiber_bottom_DV"]};
            end
            values = reshape(values,[],1);
            
            values = values(~skipBubble);
            fib_loc_ori = fib_loc_ori(~skipBubble,:);
            insig_bit = insig_bit(~skipBubble);
            
            hotzone_hull = cell(1,2);
            hotzone_com = cell(1,2);

            hotzone_color_posneg = {'y','c'};
            tick_APdim=apmldv_end(1,1):scaleUnit:apmldv_end(1,2); % just make start and end divisible by 0.1
            tick_MLdim=apmldv_end(2,1):scaleUnit:apmldv_end(2,2);
            tick_DVdim=apmldv_end(3,1):scaleUnit:apmldv_end(3,2);
            [mesh_ap,mesh_ml,mesh_dv] = meshgrid(tick_APdim,tick_MLdim,tick_DVdim);

            % get fiber contours
            tmp = fib_loc_ori;
            %   manually use 5.9 as max DV
            in_str_filter = tmp(:,1)>=APlim(1) & tmp(:,1)<=APlim(2) & tmp(:,2)>=MLlim(1) & tmp(:,2)<=MLlim(2) & tmp(:,3)>=DVlim(1) & tmp(:,3) <= 5.9;
            fib_loc_ori = fib_loc_ori(in_str_filter,:);
            values = values(in_str_filter);
            fib_loc_hull = fib_loc_ori;
            %   manually add dummy vertices to make contour more realistic
            fib_loc_hull = cat(1,fib_loc_hull,[-1.3,3.1,4.7],[-0.4,3,5.2],[0.5,2,5.5]);
            clear("tmp");
            hor_convhull = convhull(fib_loc_hull(:,1:2));
            sag_convhull = convhull(fib_loc_hull(:,[3,1]));

            % output outline of straitum if specified
            if ~isempty(output_contour_path)
                if exist(output_contour_path,"file")
                    outputfile = load(output_contour_path);
                else
                    outputfile = struct;
                end
                outputfile.straitum_contour = {[fib_loc_hull(hor_convhull,2),fib_loc_hull(hor_convhull,1)];[fib_loc_hull(sag_convhull,1),fib_loc_hull(sag_convhull,3)];};
                save(output_contour_path,'-struct',"outputfile");
            end

            % build scatterinterpolant
            %   remove nan
            if flag_include_x_as_zero
                values(isnan(values)) = 0;
                fib_loc = fib_loc_ori;
            else
                tmp = ~isnan(values);
                values = values(tmp);
                fib_loc = fib_loc_ori(tmp,:);
                insig_bit = insig_bit(tmp);
            end
            %   deal with insignificant fibers
            if flag_insig_zero
                values(insig_bit) = 0;
            end
            %   generate smoothed maps
            switch interp_method
                case "interp" % method1: interpolatoin
                    interp_f = scatteredInterpolant(fib_loc,values,'natural','none');
                    imagesc_data = interp_f(mesh_ap,mesh_ml,mesh_dv);
                case "3dconv"
                    interp_f = [];% no need

                    conv_scale_mask = common_plot_functions.get_sphere_mask(conv_radius,scaleUnit,exp_factor=exp_factor);
                    tmp_kernal = cell([length(tick_APdim),length(tick_MLdim),length(tick_DVdim)]);
                    for ri = 1:size(fib_loc,1)
                        subid = round((fib_loc(ri,:) - apmldv_end(:,1)')/scaleUnit)+1;
                        tmp_kernal{subid(1),subid(2),subid(3)} = cat(1,tmp_kernal{subid(1),subid(2),subid(3)},values(ri));
                    end
                    tmp_kernal = cellfun(@(x) mean(x,"all","omitmissing"),tmp_kernal);
                    mask = ~isnan(tmp_kernal);
                    tmp_kernal(isnan(tmp_kernal)) = 0;
                    imagesc_data = convn(tmp_kernal,conv_scale_mask,"same");
                    nan_bit = convn(mask,conv_scale_mask,"same")==0;
                    imagesc_data(nan_bit) = nan;

                    imagesc_data = permute(imagesc_data,[2,1,3]);
            end
            %---------------------debug chunk---------------------%
            % % debug plot
            % clim = max(abs(prctile(tmp,[5,95],"all")))*[-1,1];
            % bubblePlot3d(this_values,this_loc,viewAngle=[0,90],newFig=1,cmapBounds=clim,colormapOption='redblue');
            % fig=figure();
            % ax=axes(fig);
            % imagesc(ax,mean(tmp,3,"omitmissing")',clim)
            % ax.YDir="normal";
            % colormap(ax,'redblue');
            % colorbar(ax)
            % keyboard
            %---------------------debug chunk---------------------%

            % set up plotting params and plot
            if isempty(cmapBounds) && strcmp(colormapOption,"redblue")
                cmapBounds = max(abs(prctile(imagesc_data,[3,97],"all")))*[-1,1];
            elseif isempty(cmapBounds) && ~strcmp(colormapOption,"redblue")
                cmapBounds = prctile(imagesc_data,[3,97],"all")';
            end

            switch viewplane
                case "sag"
                    tmp = permute(mean(imagesc_data,1,"omitmissing"),[3,2,1]);
                    % set nan color to be grey
                    [tmp_x,tmp_y] = meshgrid(tick_APdim,tick_DVdim);
                    tmp_x = reshape(tmp_x,[],1); tmp_y = reshape(tmp_y,[],1);
                    tmp_in = inpolygon(tmp_x,tmp_y,fib_loc_hull(sag_convhull,1),fib_loc_hull(sag_convhull,3));
                    nan_mask = reshape(~tmp_in,size(tmp));
                    tmp(nan_mask) = plot_value_nan;                    
                    nan_mask = isnan(tmp);
                    % plot image
                    if plot_image
                        this_plot = imagesc(plot_ax,tick_APdim,tick_DVdim,tmp,cmapBounds);
                        image(plot_ax,tick_APdim,tick_DVdim,repmat(nan_mask,[1,1,3])*color_nan(1),AlphaData=nan_mask)
                        plot(plot_ax,fib_loc_hull(sag_convhull,1),fib_loc_hull(sag_convhull,3),'k-',LineWidth=1.5);
                        plot_ax.YDir = "reverse";plot_ax.XDir="reverse";
                        xlim(plot_ax,tick_APdim([1,end]));
                        ylim(plot_ax,tick_DVdim([1,end]));
                    end
                    % find and plot hotzone
                    posneg_bound = common_plot_functions.find_2D_hotzone(tmp,hotzone_prctile_thres);
                    [pixel_x,pixel_y] = meshgrid(tick_APdim,tick_DVdim);
                    for pni=1:2
                        bd = posneg_bound{pni};
                        hotzone_hull_tmp = cell(size(bd));
                        hotzone_com_tmp = cell(size(bd));
                        for i=1:numel(bd)
                            hotzone_hull_tmp{i} = [tick_APdim(bd{i}(:,2))',tick_DVdim(bd{i}(:,1))'];
                            % find center of mass
                            pixel_inside = inpolygon(pixel_x,pixel_y,tick_APdim(bd{i}(:,2))',tick_DVdim(bd{i}(:,1))');
                            pixel_inside_amp = tmp; pixel_inside_amp(~pixel_inside) = nan;
                            com_coord = [sum(pixel_x.*abs(pixel_inside_amp),"all","omitmissing")/sum(abs(pixel_inside_amp),"all","omitmissing"),sum(pixel_y.*abs(pixel_inside_amp),"all","omitmissing")/sum(abs(pixel_inside_amp),"all","omitmissing"),sum(abs(pixel_inside_amp),"all","omitmissing")];
                            hotzone_com_tmp{i} = com_coord;
                            if plot_hotzone
                                plot(plot_ax,tick_APdim(bd{i}(:,2)),tick_DVdim(bd{i}(:,1)),hotzone_color_posneg{pni},LineWidth=1.5);
                                scatter(plot_ax,com_coord(1),com_coord(2),100,hotzone_color_posneg{pni},'x',LineWidth=2);
                            end
                        end
                        hotzone_hull{pni} = hotzone_hull_tmp;
                        hotzone_com{pni} = hotzone_com_tmp;
                    end
                    
                case "hor"
                    tmp = mean(imagesc_data,3,"omitmissing")';
                    % set nan color to be grey
                    [tmp_x,tmp_y] = meshgrid(tick_MLdim,tick_APdim);
                    tmp_x = reshape(tmp_x,[],1); tmp_y = reshape(tmp_y,[],1);
                    tmp_in = inpolygon(tmp_x,tmp_y,fib_loc_hull(hor_convhull,2),fib_loc_hull(hor_convhull,1));
                    nan_mask = reshape(~tmp_in,size(tmp));
                    tmp(nan_mask) = plot_value_nan;
                    nan_mask = isnan(tmp);
                    % plot image
                    if plot_image
                        this_plot = imagesc(plot_ax,tick_MLdim,tick_APdim,tmp,cmapBounds);
                        image(plot_ax,tick_MLdim,tick_APdim,repmat(nan_mask,[1,1,3])*color_nan(1),AlphaData=nan_mask)
                        plot(plot_ax,fib_loc_hull(hor_convhull,2),fib_loc_hull(hor_convhull,1),'k-',LineWidth=1.5);
                        xlim(plot_ax,tick_MLdim([1,end]));
                        ylim(plot_ax,tick_APdim([1,end]));
                    end
                    % find and plot hotzone
                    posneg_bound = common_plot_functions.find_2D_hotzone(tmp,hotzone_prctile_thres);
                    [pixel_x,pixel_y] = meshgrid(tick_MLdim,tick_APdim);
                    for pni=1:2
                        bd = posneg_bound{pni};
                        hotzone_hull_tmp = cell(size(bd));
                        hotzone_com_tmp = cell(size(bd));
                        for i=1:numel(bd)
                            hotzone_hull_tmp{i} = [tick_MLdim(bd{i}(:,2))',tick_APdim(bd{i}(:,1))'];
                            % find center of mass
                            pixel_inside = inpolygon(pixel_x,pixel_y,tick_MLdim(bd{i}(:,2))',tick_APdim(bd{i}(:,1))');
                            pixel_inside_amp = tmp; pixel_inside_amp(~pixel_inside) = nan;
                            com_coord = [sum(pixel_x.*abs(pixel_inside_amp),"all","omitmissing")/sum(abs(pixel_inside_amp),"all","omitmissing"),sum(pixel_y.*abs(pixel_inside_amp),"all","omitmissing")/sum(abs(pixel_inside_amp),"all","omitmissing"),sum(abs(pixel_inside_amp),"all","omitmissing")];
                            hotzone_com_tmp{i} = com_coord;
                            if plot_hotzone
                                plot(plot_ax,tick_MLdim(bd{i}(:,2)),tick_APdim(bd{i}(:,1)),hotzone_color_posneg{pni},LineWidth=1.5);
                                scatter(plot_ax,com_coord(1),com_coord(2),100,hotzone_color_posneg{pni},'x',LineWidth=2);
                            end
                        end
                        hotzone_hull{pni} = hotzone_hull_tmp;
                        hotzone_com{pni} = hotzone_com_tmp;
                    end
                    
            end

            if plot_image
                cmap = colormap(plot_ax,colormapOption);
                colorbar(plot_ax,Ticks=linspace(cmapBounds(1),cmapBounds(2),7),TickLabels=arrayfun(@(x) sprintf("%0.4f",x),linspace(cmapBounds(1),cmapBounds(2),7)));
            end
        end
        
        function posneg = find_2D_hotzone(input,hotzone_prctile_thres,varargin)
            min_area = 4;

            ip = inputParser;
            ip.addParameter('min_area',16);
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            % find largest positive and negative hotzones
            posneg = cell(1,2);
            for ri=1:2
                if ri==1 % positive hotzone
                    this_input = input;
                else % negative hotzone
                    this_input = -input;
                end
                pos_thres = prctile(this_input(this_input>0),hotzone_prctile_thres,"all");
                high_region = this_input>pos_thres;
                [b_pos,l_pos,n_pos] = bwboundaries(high_region,8,"noholes");
                area = nan(1,n_pos);
                for i=1:n_pos
                    area(i) = sum(l_pos==i,"all");
                end
                [area_sorted,tmp] = sort(area,"descend");
                tmp = tmp(area_sorted>=min_area);
                if ~isempty(tmp)
                    posneg{ri} = b_pos(tmp);
                else
                    posneg{ri} = {double.empty(0,2)};
                end
            end
        end

        function output = get_sphere_mask(radius,unit,varargin)
            exp_factor = 3.5;
            ip = inputParser;
            ip.addParameter('exp_factor',3.5)
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            d1 = -radius:unit:radius;
            [m1,m2,m3] = meshgrid(d1,d1,d1);
            tmp = arrayfun(@(x,y,z) sqrt(x^2+y^2+z^2),m1,m2,m3);
            output = arrayfun(@(x) exp(-exp_factor*x), tmp);
        end
            
        function mdl = lme_3Dgradient(value,ct_table)
            value = normalize(value);
            model_tbl = table(Size=[length(value),5],VariableType=["double","double","double","double","string"],VariableNames=["value","AP","ML","DV","mouse"]);
            model_tbl.value = reshape(value,[],1);
            model_tbl(:,["AP","ML","DV","mouse"]) = ct_table(:,["fiber_bottom_AP","fiber_bottom_ML","fiber_bottom_DV","mouse_name"]);
            % model_tbl{:,["AP","ML","DV"]} = normalize(model_tbl{:,["AP","ML","DV"]},1);
            mdl = fitlme(model_tbl,'value~AP+ML+DV+(AP|mouse)+(ML|mouse)+(DV|mouse)+(1|mouse)');
        end

        function mdl = lme_3Dgradient_interact(value,ct_table)
            assert(size(value,2)==2,"For now only compare two groups.")
            model_tbls = {};
            value = normalize(value,1);
            for i=1:size(value,2)    
                model_tbl = table(Size=[size(value,1),6],VariableType=["double","double","double","double","string","string"],VariableNames=["value","AP","ML","DV","mouse","trialType"]);
                model_tbl.value = value(:,i);
                model_tbl(:,["AP","ML","DV","mouse"]) = ct_table(:,["fiber_bottom_AP","fiber_bottom_ML","fiber_bottom_DV","mouse_name"]);
                % model_tbl{:,["AP","ML","DV"]} = normalize(model_tbl{:,["AP","ML","DV"]},1);
                model_tbl.trialType = repmat(string(i),size(value,1),1);
                model_tbls =cat(1,model_tbls,model_tbl);
            end
            mdl = fitlme(model_tbls,'value~trialType*(AP+ML+DV)+(trialType+AP+ML+DV|mouse)+(1|mouse)');
        end
        
        function plot_3Dgradient(ax,APMLDV_coef,view,varargin)
            lcolor = 'black';

            ip = inputParser;
            ip.addParameter('lcolor','black') % change this to change plot density
            
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            APMLDV_coef = reshape(APMLDV_coef,1,[]);
            % fakez=1;fakex=0;
            this_center = [0,2,4];
            switch view
                case "hor"
                    length = 1;
                    this_coef = APMLDV_coef([2,1]);
                    this_coef = this_coef/sqrt(sum(this_coef.^2));
                    endpoint = [2,0] + [-this_coef;this_coef]*length;
                    plot(ax,endpoint(:,1),endpoint(:,2),Color=lcolor,LineWidth=2);
                    scatter(ax,endpoint(2,1),endpoint(2,2),100,lcolor,'filled','o');
                case "sag"
                    length = 1;
                    this_coef = APMLDV_coef([1,3]);
                    this_coef = this_coef/sqrt(sum(this_coef.^2));
                    endpoint = [0,4] + [-this_coef;this_coef]*length;
                    plot(ax,endpoint(:,1),endpoint(:,2),Color=lcolor,LineWidth=2);
                    scatter(ax,endpoint(2,1),endpoint(2,2),100,lcolor,'filled','o');
            end
        end
        
        function modality_between_context(cf,plot_tag,mouse_names,save_name,varargin)
            % plot_tag = ["cue_vs_cue_learning","cue_vs_cue_extinction","cue_vs_rew","extinction_vs_lick"];
            % mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
            save_name = char(save_name);
            tabular_sheet_name = '';

            ranksum_alpha = 0.01;
            ip = inputParser;
            ip.addParameter("ranksum_alpha",0.01)
            ip.addParameter("tabular_sheet_name",'')
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            if isempty(tabular_sheet_name)
                tabular_sheet_name = save_name;
            end
            
            ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
            c_names = ["LED","Tone"];
            comp_names = ["pk","dp","re"];
            omi_fname = ["LEDomi","Toneomi"];
        
            nr=2;nc=3;
            fig = figure(Position = [1,41,2560,1323]);
            tiled = tiledlayout(fig,nr,nc);
            sgtitle(tiled,plot_tag);
            axs = gobjects(1,nr*nc); for i=1:nr*nc;axs(i)=nexttile(tiled,i);end
            hold(axs,"on")
            for c_i = 1:2
                switch plot_tag
                    case "cue_vs_cue_learning"
                        tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                        cue_rews = ["cue1","cue2"];                    
                        tas12 = {tas,tas};
                        test_fieldname = ["cue"+c_i+"early","cue"+c_i+"late";...
                            "cue2early","cue2late";];
                        cue_rew_fix1 = ["cue1","cue2"];
                        cue_rew_fix2 = "cue2";
                        test_ranksum_alphas = [0.01,0.01];
                        comp_sig_direction = {@gt,@lt,@gt;@gt,@lt,@gt;};
                        fisher_data_comp = [1,2];
                    case "cue_vs_cue_extinction"
                        tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                        cue_rews = ["cue1","cue2"];
                        tas12 = {tas,tas};
                        test_fieldname = ["cue"+c_i+"late","cue"+c_i+omi_fname(c_i);...
                            "cue2late","cue2Toneomi";];
                        cue_rew_fix1 = ["cue1","cue2"];
                        cue_rew_fix2 = "cue2";
                        test_ranksum_alphas = [0.01,0.01];
                        comp_sig_direction = {@lt,@gt,@gt;@lt,@gt,@gt;};
                        fisher_data_comp = [1,2];
                    case "cue_vs_rew"
                        % tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                        tas = load("F:\Safa_Processed\components_GLM_highpass03_wGLM_v2\across\_components_window\tri_avg_single_filtered_rebase_rew_consump.mat");
                        tas_rew = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump_rew_merged.mat']);
                        cue_rews = ["cue","rew"];
                        tas12 = {tas,tas_rew};
                        test_fieldname = ["cue"+c_i+"early","cue"+c_i+"late";...
                            "rew1early","rew1late";];
                        cue_rew_fix1 = ["cue1","cue2"];
                        cue_rew_fix2 = "rew";
                        test_ranksum_alphas = [0.01,0.01];
                        comp_sig_direction = {@gt,@lt,@gt;@lt,@gt,@gt;};
                        fisher_data_comp = [1,2];
                    case "extinction_vs_lick"
                        tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
                        tas_lick = load([cf,'processed_and_organized_data\ITI_licking_tri_avg_single_filtered.mat']);
                        cue_rews = ["extinct","lick"];
                        tas12 = {tas,tas_lick};
                        test_fieldname = ["cue"+c_i+"late","cue"+c_i+omi_fname(c_i);...
                            "early","late";];
                        cue_rew_fix1 = ["cue1","cue2"];
                        cue_rew_fix2 = "ITI_lick";
                        test_ranksum_alphas = [0.01,0.025];
                        % comp_sig_direction = {@lt,@gt,@gt;@lt,@gt,@lt;};
                        comp_sig_direction = {@gt,@gt,@gt;@lt,@gt,@lt;};
                        fisher_data_comp = [1,3];
                end
        
                fig_cc = figure(Position = [1,41,2560,1323]);
                tiled_cc = tiledlayout(fig_cc,nr,nc,TileSpacing="compact");
                sgtitle(tiled_cc,plot_tag+" "+c_names(c_i));
                axs_cc = gobjects(1,nr*nc); for i=1:nr*nc;axs_cc(i)=nexttile(tiled_cc,i);end
                hold(axs_cc,"on")   
        
                late_early_data = struct;
                late_early_data.across.(cue_rews(1)) = cell(3,1);
                late_early_data.across.(cue_rews(2)) = cell(3,1);
                late_early_data.across.(cue_rews(1)+"_sig") = cell(3,1);
                late_early_data.across.(cue_rews(2)+"_sig") = cell(3,1);
                for mouse_name = mouse_names
                    n_rois = size(tas.cue1early.(mouse_name).mu_location_value,2);
                    for cue_rew_i = 1:2
                        cue_rew = cue_rews(cue_rew_i);
                        this_tas = tas12{cue_rew_i};
                        test_ranksum_alpha = test_ranksum_alphas(cue_rew_i);
                        for compi=1:3
                            value_early = this_tas.(test_fieldname(cue_rew_i,1)).(mouse_name).mu_location_value(compi,:,1);
                            sig_early = this_tas.(test_fieldname(cue_rew_i,1)).(mouse_name).mu_location_value_significance(compi,:);
                            single_early = num2cell(this_tas.(test_fieldname(cue_rew_i,1)).(mouse_name).single_values(compi,:,:),3);
        
                            value_late = this_tas.(test_fieldname(cue_rew_i,2)).(mouse_name).mu_location_value(compi,:,1);
                            sig_late = this_tas.(test_fieldname(cue_rew_i,2)).(mouse_name).mu_location_value_significance(compi,:);
                            single_late =  num2cell(this_tas.(test_fieldname(cue_rew_i,2)).(mouse_name).single_values(compi,:,:),3);
                            
                            diff_struct = diff_phase_ranksum(value_late,sig_late,single_late,value_early,sig_early,single_early,[]);
        
                            % get significant of differences
                            ranksum_p = diff_struct.ranksum_p';
                            value_late_early = diff_struct.diff_mu';
                            
                            late_early = value_late_early;
                            late_early(~logical(sig_early)&~logical(sig_late)) = nan; % they are already all nan for this version of tas
                            
                            late_early_data.(mouse_name).(cue_rew){compi} = late_early;
                            late_early_data.(mouse_name).(cue_rew+"_sig"){compi} = ranksum_p <= test_ranksum_alpha;
                
                            late_early_data.across.(cue_rew){compi} = cat(2,late_early_data.across.(cue_rew){compi},late_early);
                            late_early_data.across.(cue_rew+"_sig"){compi} = cat(2,late_early_data.across.(cue_rew+"_sig"){compi},ranksum_p <= test_ranksum_alpha);
                        end
                    end
                end
                position_tag = "";
                position_filter_bit = struct;
                position_filter_bit.across = [];
                position_filter_bit.ct_across = {};
                for mouse_name = mouse_names
                    tmp = ct_table{ct_table{:,"mouse_name"}==mouse_name,"significance"}==1;
                    position_filter_bit.(mouse_name) = tmp';
                    position_filter_bit.across = cat(2,position_filter_bit.across,tmp');
                    position_filter_bit.ct_across = cat(1,position_filter_bit.ct_across,ct_table(ct_table.mouse_name==mouse_name&ct_table.significance==1,:));
                end
        
                % plot cue_dip_change vs rew_dip_change
                mouse_name = "across";
                fisher_data = nan(2,2);
                for compi=1:3
                    cue_total = late_early_data.(mouse_name).(cue_rews(1)){compi};
                    cue_sig_total = late_early_data.(mouse_name).(cue_rews(1)+"_sig"){compi};
                    rew_total = late_early_data.(mouse_name).(cue_rews(2)){compi};
                    rew_sig_total = late_early_data.(mouse_name).(cue_rews(2)+"_sig"){compi};
                    % apply spatial filter
                    cue_total = cue_total(logical(position_filter_bit.(mouse_name)));
                    cue_sig_total = cue_sig_total(logical(position_filter_bit.(mouse_name)));
                    rew_total = rew_total(logical(position_filter_bit.(mouse_name)));
                    rew_sig_total =  rew_sig_total(logical(position_filter_bit.(mouse_name)));
                    
                    % plot catergorized circle plot
                    tmp_fhandle_1 = comp_sig_direction{1,compi}; tmp_fhandle_2 = comp_sig_direction{2,compi};
                    circle_cmap = [0.7,0.7,0.7;0,0,1;1,0,0;0,1,0;];
                    circle_value = zeros(1,size(position_filter_bit.ct_across,1));
        
                    sig_combinations_total = [...
                        cue_sig_total & tmp_fhandle_1(cue_total,0);...
                        rew_sig_total & tmp_fhandle_2(rew_total,0);...
                        ~cue_sig_total & ~rew_sig_total;...
                        ];
                    tmp = sig_combinations_total(1,:) & sig_combinations_total(2,:);
                    sig_combinations_total_rpe = [...
                        tmp;...
                        sig_combinations_total(1,:) & ~tmp;...
                        sig_combinations_total(2,:) & ~tmp;...
                        sig_combinations_total(3,:);...
                        ];
        
                    % ------------------ debug -------------------%
                    % if c_i==1 && compi == 1
                    %     test_v11 = sig_combinations_total(1,:)& tmp_fhandle_1(cue_total,0) & tmp_fhandle_2(rew_total,0);
                    %     test_v12 = sig_combinations_total(3,:)& tmp_fhandle_2(rew_total,0);
                    %     test_v13 = sig_combinations_total;
                    %     test_v14 = rew_total;
                    %     % keyboard
                    % elseif c_i==2 && compi == 1
                    %     test_v21 = sig_combinations_total(1,:)& tmp_fhandle_1(cue_total,0) & tmp_fhandle_2(rew_total,0);
                    %     test_v22 = sig_combinations_total(3,:)& tmp_fhandle_2(rew_total,0);
                    %     test_v23 = sig_combinations_total;
                    %     test_v24 = rew_total;
                    %     keyboard
                    % end
                    % ------------------ debug -------------------%
        
                    circle_value(sig_combinations_total_rpe(1,:)) = 1;
                    circle_value(sig_combinations_total_rpe(2,:)) = 2;
                    circle_value(sig_combinations_total_rpe(3,:)) = 3;
                    circle_value(sig_combinations_total_rpe(4,:)) = nan;
                    for ci=0:3
                        tmp = circle_value==ci;
                        if sum(tmp)==0
                            continue;
                        end
                        plot_functions.scatter_3d_with_datatip(axs_cc(compi),ones(1,sum(tmp)),position_filter_bit.ct_across(tmp,:),...
                            viewAngle=[0,90],colormapOption=circle_cmap(ci+1,:),sorting_tag="location");
                        plot_functions.scatter_3d_with_datatip(axs_cc(compi+nc),ones(1,sum(tmp)),position_filter_bit.ct_across(tmp,:),...
                            viewAngle=[-90,0],colormapOption=circle_cmap(ci+1,:),sorting_tag="location");
                    end
                    tmp = isnan(circle_value);
                    plot_functions.scatter_3d_with_datatip(axs_cc(compi),nan(1,sum(tmp)),position_filter_bit.ct_across(tmp,:),...
                        viewAngle=[0,90],colormapOption=circle_cmap(ci+1,:),sorting_tag="location");
                    plot_functions.scatter_3d_with_datatip(axs_cc(compi+nc),nan(1,sum(tmp)),position_filter_bit.ct_across(tmp,:),...
                        viewAngle=[-90,0],colormapOption=circle_cmap(ci+1,:),sorting_tag="location");
        
                    title(axs_cc(compi+nc),["red:"+cue_rew_fix1(c_i)+" green:"+cue_rew_fix2+" blue:both"])
                    
                    ax = axs(nc*(c_i-1)+compi);
                    tmp = sum(sig_combinations_total_rpe,2)';tmp=[tmp([2,1,3]),size(sig_combinations_total_rpe,2)-sum(tmp(1:3))];
                    if compi==fisher_data_comp(1)
                        fisher_data(1,:) = tmp(1:2);
                    elseif compi==fisher_data_comp(2)
                        fisher_data(2,:) = tmp(1:2);
                    end
                    % histogram(ax,Categories={'cue only','cue & rew','rew only','opposite or nan'},BinCounts=tmp,Normalization="pdf")
                    pie_prc = tmp/sum(tmp)*100;
                    formattedLabels = strcat({[char(cue_rew_fix1(c_i)),' only'],[char(cue_rew_fix1(c_i)),' and ',char(cue_rew_fix2)],[char(cue_rew_fix2),' only'],'opposite or nan'}, " (", arrayfun(@(x,y) sprintf("%d->%0.1f",y,x),pie_prc,tmp), "%)"); 
                    ph=bar(ax,pie_prc/100,FaceColor="flat"); ax.XTick=1:4; ax.XTickLabel=formattedLabels;
                    ph.CData = [1,0,0;0,0,1;0,1,0;0.7,0.7,0.7;];
                    % axis(ax,"vis3d");axis(ax,"off");
                    ylim(ax,[0,1]);
                    title(ax,c_names(c_i)+" "+comp_names(compi));
        
                end
                [~, p] = fishertest(fisher_data);
                title(axs(nc*(c_i-1)+2),[c_names(c_i)+" "+comp_names(2);"fishertest p = "+sprintf("%0.1e",p)+" between "+strjoin(comp_names(fisher_data_comp)," and ")]);
        
                hold(axs_cc,"off")
                saveas(fig_cc,[cf,save_name,'_',char(plot_tag),'_',char(c_names(c_i)),'.png'])
                delete(fig_cc)
            end
            hold(axs,"off")
            saveas(fig,[cf,save_name,'_',char(plot_tag),'_histogram.png'])
            delete(fig)

            % helper function
            function diff_struct = diff_phase_ranksum(mu1,sig1,single1,mu2,sig2,single2,in_striatum_bit)
                mu1(~sig1) = nan; mu2(~sig2) = nan;
                single1(~sig1) = cellfun(@(x) nan(size(x)),single1(~sig1),UniformOutput=false);
                single2(~sig2) = cellfun(@(x) nan(size(x)),single2(~sig2),UniformOutput=false);
                
                this_value = common_functions.nan_minus(mu1,mu2);
                this_ranksum_p = common_functions.nan_ranksum_cell(single1,single2);
                diff_struct.included_bit = in_striatum_bit;
                diff_struct.diff_mu = this_value';
                diff_struct.ranksum_p = this_ranksum_p';
            end
        end

    end % Static methods end
end