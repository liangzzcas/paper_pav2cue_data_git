%% figure 3A
close all;clear;clc;
cf = [pwd,'\'];
% cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump_rew_merged.mat']);
cwa = load("E:\Safa_Processed\#paper\#paper_figures\_data\highpass03_wGLM_v2\delivery\components_window_activity_filtered_rebase_rew_merged.mat");
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
%% 
% plot_info_1 = {"G19",49,18;"G24",23,18;"G15",6,30;"G15",5,30;"G22",40,18;"G22",28,18;};
% plot_info_1 = {"G24",23,18;"G22",40,18;"G22",28,18;"G22",58,18;"G22",72,18;};
% plot_info_1 = {"G24",23,18;"G22",40,18;"G22",28,18;};
% plot_trace_clim = {0.015;0.025;0.025;0.02;0.015;...
%     0.025;0.025;};
% plot_imagesc_clim = {0.02;0.04;0.04;0.04;0.025;...
%     0.04;0.04;};

% plot_info_1 = {"G12",39,30;"G12",21,30;"G21",27,18;"G23",48,18;"G22",28,18;"G24",23,18;};
% plot_trace_clim = {[-0.01,0.02];[-0.01,0.045];[-0.025,0.025];[-0.01,0.015];[-0.025,0.03];[-0.015,0.015];};
% plot_imagesc_clim = {0.03;0.06;0.04;0.02;0.04;0.02;};
% plot_info_2 = ["rew1early","rew1late","rew1LEDomi","rew1Toneomi","unpredlate"];
% plot_colors = lines(length(plot_info_2));

% plot_info_1 = {"G17",31,18;"G17",42,18;"G21",33,18;"G23",13,18;};
plot_info_1 = {"G23",36,18;"G23",54,18;"G21",27,18;"G21",42,18;"G23",55,18;};

plot_trace_clim = {[-0.01,0.02];[-0.01,0.045];[-0.025,0.025];[-0.01,0.015];[-0.025,0.03];[-0.015,0.015];[-0.015,0.015];[-0.015,0.015];[-0.015,0.015];};
plot_imagesc_clim = {0.03;0.06;0.04;0.02;0.04;0.02;0.04;0.04;0.04;};
plot_info_2 = ["cue2late","cue2Toneomi"];
plot_colors = lines(length(plot_info_2));

fig = figure(Position=[40,40,1800,1310]);
tiled = tiledlayout(fig,length(plot_info_2)+1,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,(length(plot_info_2)+1)*size(plot_info_1,1));
for i=1:(length(plot_info_2)+1)*size(plot_info_1,1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axes([pi+size(plot_info_1,1).*[0:length(plot_info_2)]]);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    
    plot_x = (1:sr*4)/sr-1;
    nts = nan(1,length(plot_info_2));
    lgs = gobjects(1,length(plot_info_2));
    for ni=1:length(plot_info_2)
        tmp = cwa.(plot_info_2(ni)).(this_info{1}).activity(:,this_info{2}+3,:); nts(ni) = size(tmp,3);
        lgs(ni)=common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(ni,:)*0.8,sem_color=plot_colors(ni,:));
        imagesc(this_axes(ni+1),plot_x,1:nts(ni),permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
        this_axes(ni+1).YDir = "reverse";
        colormap(this_axes(ni+1),common_functions.redblue());
        colorbar(this_axes(ni+1));
        coord={ct_table{ct_table.mouse_name==this_info{1} & ct_table.ROI_original==this_info{2},["fiber_bottom_AP","fiber_bottom_ML","fiber_bottom_DV"]}};
        title(this_axes(1),this_info{1}+" ROI "+this_info{2}+sprintf(" AP/ML/DV=%0.2f/%0.2f/%0.2f",coord{:}))
    end

    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-0.5,2]);
    % ylim(this_axes(1),[-1,1]*plot_trace_clim{pi});
    ylim(this_axes(1),plot_trace_clim{pi});
    for axi=1:length(plot_info_2)
        ax=this_axes(axi+1);
        xline(ax,0,'--');
        xlim(ax,plot_x([1,end]));
        ylim(this_axes(axi+1),[1,nts(axi)]);
        ylabel(this_axes(axi+1),plot_info_2(axi))
    end
end
hold(axes,"off");
legend(axes(1),lgs,plot_info_2);
% saveas(fig,[cf,'fig3A.fig']);
saveas(fig,'F:\Safa_Processed\#paper_figure\#update_review\_revision_plot_2\test.fig');
saveas(fig,'F:\Safa_Processed\#paper_figure\#update_review\_revision_plot_2\test.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 3B
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_time_histogram(cf,["rew1early","rew1late"],'fig3B');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 3C
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_presence_circle(cf,"pav2cue task",["rew1early","rew1late"],'fig3C');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 3D
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_rewmerged",{["rew1late","rew1early"]},mouse_names,'fig3D');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 3E
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
c_names = ["LED","Tone"];
ranksum_alpha = 0.05;
for c_i = 1:2
    dip_late_early_data = struct;
    dip_late_early_data.across.cue = [];
    dip_late_early_data.across.rew = [];
    dip_late_early_data.across.cue_sig = [];
    dip_late_early_data.across.rew_sig = [];
    for mouse_name = mouse_names
        n_rois = size(tas.cue1early.(mouse_name).mu_location_value,2);
        for cue_rew = ["cue","rew"]
            dip_value_early = tas.(cue_rew+c_i+"early").(mouse_name).mu_location_value(2,:,1);
            dip_sig_early = tas.(cue_rew+c_i+"early").(mouse_name).mu_location_value_significance(2,:);
            
            dip_value_late = tas.(cue_rew+c_i+"late").(mouse_name).mu_location_value(2,:,1);
            dip_sig_late = tas.(cue_rew+c_i+"late").(mouse_name).mu_location_value_significance(2,:);

            % get significant of differences
            ranksum_p = nan(1,n_rois);
            for r = 1:n_rois
                ranksum_p(r) = common_functions.nan_ranksum(...
                    squeeze(tas.(cue_rew+c_i+"early").(mouse_name).single_values(2,r,:))*dip_sig_early(r),...
                    squeeze(tas.(cue_rew+c_i+"late").(mouse_name).single_values(2,r,:))*dip_sig_late(r));
            end

            dip_value_late_early = common_functions.nan_minus(dip_value_late.*dip_sig_late,dip_value_early.*dip_sig_early);
            dip_late_early = dip_value_late_early;
            dip_late_early(~logical(dip_sig_early)&~logical(dip_sig_late)) = nan; % they are already all nan for this version of tas
            
            dip_late_early_data.(mouse_name).(cue_rew) = dip_late_early;
            dip_late_early_data.(mouse_name).(cue_rew+"_sig") = ranksum_p <= ranksum_alpha;

            dip_late_early_data.across.(cue_rew) = cat(2,dip_late_early_data.across.(cue_rew),dip_late_early);
            dip_late_early_data.across.(cue_rew+"_sig") = cat(2,dip_late_early_data.across.(cue_rew+"_sig"),ranksum_p <= ranksum_alpha);
        end
    end

    % get position filter
    anterior_only_bit = struct;
    anterior_only_bit.across = [];
    for mouse_name = mouse_names
        tmp = ct_table{ct_table{:,"mouse_name"}==mouse_name,"fiber_bottom_AP"} >= -0.5;
        anterior_only_bit.(mouse_name) = tmp';
        anterior_only_bit.across = cat(2,anterior_only_bit.across,tmp');
    end

    % plot cue_dip_change vs rew_dip_change
    % for mouse_name = [mouse_names,"across"]
    for mouse_name = ["across"]
        cue = dip_late_early_data.(mouse_name).cue;
        cue_sig = dip_late_early_data.(mouse_name).cue_sig;
        rew = dip_late_early_data.(mouse_name).rew;
        rew_sig = dip_late_early_data.(mouse_name).rew_sig;
        % apply spatial filter
        cue = cue(logical(anterior_only_bit.(mouse_name)));
        cue_sig = cue_sig(logical(anterior_only_bit.(mouse_name)));
        rew = rew(logical(anterior_only_bit.(mouse_name)));
        rew_sig =  rew_sig(logical(anterior_only_bit.(mouse_name)));
        % repalce nan by 0
        not_neither_nan = ~(isnan(cue) & isnan(rew));
        cue = cue(not_neither_nan); cue(isnan(cue)) = 0;
        cue_sig = cue_sig(not_neither_nan);
        rew = rew(not_neither_nan); rew(isnan(rew)) = 0;
        rew_sig = rew_sig(not_neither_nan);

        counts = [sum(cue<0 & rew>0),sum(cue==0 & rew>0),sum(cue>0 & rew>0),...
            sum(cue<0 & rew==0),sum(cue>0 & rew==0),...
            sum(cue<0 & rew<0),sum(cue==0 & rew<0),sum(cue>0 & rew<0)];
        fractions = counts/length(cue);
        c_f = counts+"("+arrayfun(@(x) sprintf("%0.2f",x*100),fractions)+"%)";
        
        sig_combinations = logical([cue_sig&rew_sig;cue_sig&~rew_sig;~cue_sig&rew_sig;~cue_sig&~rew_sig]);
        sig_combinations_style = {'bo','ro','go','kx'};
        
        % get linear fit of only blue ones
        lm_Y = rew(sig_combinations(1,:))';
        lm_X = cue(sig_combinations(1,:))';
        mdl = fitlm(lm_X,lm_Y,RobustOpts="on");
        mdl_R2 = mdl.Rsquared.Ordinary; mdl_p = anova(mdl,"summary"); mdl_p = mdl_p{"Model","pValue"}; mdl_slope = mdl.Coefficients{"x1","Estimate"};
        xpred = linspace(min(cue,[],"all","omitmissing"),max(cue,[],"all","omitmissing"),20)';
        [ypred,yci] = predict(mdl,xpred);

        lgs = gobjects(1,5);
        fig = figure(Position = [100,100,800,600]);
        sgtitle(fig,sprintf("sig for both cue rew linear fit: R^2 = %0.2f, p = %0.3e, slope = %0.3f",mdl_R2,mdl_p,mdl_slope))
        ax = axes(fig);
        hold(ax,"on")
        for sig_comb_i = 1:3
            lgs(sig_comb_i) = scatter(ax,cue(sig_combinations(sig_comb_i,:)),rew(sig_combinations(sig_comb_i,:)),100,lines(1),"filled",sig_combinations_style{sig_comb_i});
        end
        lgs(4) = plot(ax,xpred,ypred,'b-');
        tmp = plot(ax,xpred,yci,'b--');
        lgs(5) = tmp(1);
        % sig_comb_i= 4;
        % scatter(ax,cue(sig_combinations(sig_comb_i,:)),rew(sig_combinations(sig_comb_i,:)),100,lines(1),sig_combinations_style{sig_comb_i},LineWidth=1.8);
        hold(ax,"off")
        xline(ax,0)
        yline(ax,0)
        xlabel(sprintf("%s\n%s","diff cue dip",strjoin(c_f,"; ")))
        ylabel(["diff rew dip","AP >= -0.5"])
        legend(ax,lgs(4:5),["linear fit predict of blue dots","ci of linear fit"])
        saveas(fig,[cf,'fig3E_',char(c_names(c_i)),'.png'])
        delete(fig)
    end
end