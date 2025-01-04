%% figure 6A
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_time_histogram(cf,["cue2Toneomi","cue2LEDomi"],'sup_fig6A');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 6B
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_presence_circle(cf,"pav2cue task",["cue2Toneomi"],'sup_fig6B');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 6C
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
plot_info_1 = {"G23",36,18;"G23",48,18;"G21",41,18;"G21",9,18;"G12",5,30;};
plot_info_2 = ["cue2late","cue2LEDomi","cue2Toneomi"];
plot_colors = lines(3);

fig = figure(Position=[100,100,1600,300]);
nr=1;
tiled = tiledlayout(fig,nr,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,nr*size(plot_info_1,1));
for i=1:nr*size(plot_info_1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axes(pi);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    
    plot_x = (1:sr*4)/sr-1;
    tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
    
    tmp = cwa.(plot_info_2(2)).(this_info{1}).activity(:,this_info{2}+3,:);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));

    tmp = cwa.(plot_info_2(3)).(this_info{1}).activity(:,this_info{2}+3,:);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(3,:)*0.8,sem_color=plot_colors(3,:));
   
    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-0.5,2]);
end
hold(axes,"off");
saveas(fig,'sup_fig6C.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 6D
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue2Toneomi","cue2late"]},mouse_names,'fig6D',cmap_limit={[-1,1]*0.011,[-1,1]*0.012,[-1,1]*0.011});


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 6E-F
% Incorporated into figure 4K-L
% Refer to section "figure 4K-L" in file "figure4_code.m"


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 6E-F
% Incorporated into figure 4M
% Refer to section "figure 4M" in file "figure4_code.m"


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 6H
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
cwa_data = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];

suffix_text = "Tone diff vs LED diff";
comp_names = ["pk","dp","re"];
ranksum_p_value_thres = 0.05;
cue_omission_tags = ["LEDomi","Toneomi"];

across_data = struct;
for c_i = 1:2
    across_data.("cue"+c_i).diff = [];
    across_data.("cue"+c_i).ranksum_p = [];
end

for mouse_name_i = 1:length(mouse_names)
    mouse_name = mouse_names(mouse_name_i);
    sig_rois = find(ct_table{ct_table{:,"mouse_name"}==mouse_name,"significance"})';
    cue12 = struct;
    for c_i = 1:2
        this_mu_tone_e = tas.("cue"+c_i+"late").(mouse_name).mu_location_value(:,:,1);
        this_mu_tone_single_e = tas.("cue"+c_i+"late").(mouse_name).single_values;
        this_sig_tone_e = tas.("cue"+c_i+"late").(mouse_name).mu_location_value_significance;
        this_mu_tone_e = this_mu_tone_e .* this_sig_tone_e;

        this_mu_tone_l = tas.("cue"+c_i+cue_omission_tags(c_i)).(mouse_name).mu_location_value(:,:,1);
        this_mu_tone_single_l = tas.("cue"+c_i+cue_omission_tags(c_i)).(mouse_name).single_values;
        this_sig_tone_l = tas.("cue"+c_i+cue_omission_tags(c_i)).(mouse_name).mu_location_value_significance;
        this_mu_tone_l = this_mu_tone_l .* this_sig_tone_l;

        this_diff = common_functions.nan_minus(this_mu_tone_l,this_mu_tone_e);

        diff_sig = nan(3,size(this_mu_tone_e,2));
        for r_i = 1:size(this_mu_tone_e,2)
            for comp_i = 1:3
                if xor(this_sig_tone_e(comp_i,r_i),this_sig_tone_l(comp_i,r_i))
                    diff_sig(comp_i,r_i) = -1;
                elseif this_sig_tone_e(comp_i,r_i) && this_sig_tone_l(comp_i,r_i)
                    p = ranksum(squeeze(this_mu_tone_single_e(comp_i,r_i,:)),squeeze(this_mu_tone_single_l(comp_i,r_i,:)));
                    diff_sig(comp_i,r_i) = p;
                end
            end
        end

        cue12.("cue"+c_i).diff = this_diff(:,sig_rois);
        cue12.("cue"+c_i).ranksum_p = diff_sig(:,sig_rois);

        across_data.("cue"+c_i).diff = cat(2,across_data.("cue"+c_i).diff,cue12.("cue"+c_i).diff);
        across_data.("cue"+c_i).ranksum_p = cat(2,across_data.("cue"+c_i).ranksum_p,cue12.("cue"+c_i).ranksum_p);
    end
end

% fit and plot across data

fig_fc = figure(Position=[100,100,1600,650],InvertHardcopy="off");
tiled_fc = tiledlayout(1,3,Parent=fig_fc,TileSpacing="tight");
sgtitle(tiled_fc,"across "+suffix_text)
axs = gobjects(1,3);
for i = 1:3
    axs(i) = nexttile(tiled_fc,i);
end
hold(axs,"on")
for i = 1:3
    ax = axs(i);
    xline(ax,0)
    yline(ax,0)
    title(ax,comp_names(i))
    ylabel(ax,"Tone diff")

    sig_1 = across_data.cue1.ranksum_p(i,:) <= ranksum_p_value_thres;
    sig_2 = across_data.cue2.ranksum_p(i,:) <= ranksum_p_value_thres;

    s1=scatter(ax,across_data.cue1.diff(i,sig_1 & sig_2),across_data.cue2.diff(i,sig_1 & sig_2),100,lines(1),"filled","o");
    s2=scatter(ax,across_data.cue1.diff(i,sig_1 & ~sig_2),across_data.cue2.diff(i,sig_1 & ~sig_2),100,"filled",'ro',LineWidth=2);
    s3=scatter(ax,across_data.cue1.diff(i,~sig_1 & sig_2),across_data.cue2.diff(i,~sig_1 & sig_2),100,"filled",'go',LineWidth=2);

    tmp_x = across_data.cue1.diff(i,sig_1|sig_2);
    tmp_y = across_data.cue2.diff(i,sig_1|sig_2);

    count = [sum(tmp_x<0&tmp_y>0),sum(tmp_x==0&tmp_y>0),sum(tmp_x>0&tmp_y>0),...
        sum(tmp_x<0&tmp_y==0),sum(tmp_x>0&tmp_y==0),...
        sum(tmp_x<0&tmp_y<0),sum(tmp_x==0&tmp_y<0),sum(tmp_x>0&tmp_y<0)];
    frac = count/length(tmp_x)*100;
    count_frac_1 = strjoin(arrayfun(@(x,y) sprintf("%d(%0.3f%%)",x,y),count(1:3),frac(1:3)),",");
    count_frac_2 = strjoin(arrayfun(@(x,y) sprintf("%d(%0.3f%%)",x,y),count(4:5),frac(4:5)),",");
    count_frac_3 = strjoin(arrayfun(@(x,y) sprintf("%d(%0.3f%%)",x,y),count(6:8),frac(6:8)),",");
    xlabel(ax,["LED diff",count_frac_1,count_frac_2,count_frac_3])
end
hold(axs,"off")
linkaxes(axs,"xy")
legend([s1,s2,s3],["both","LED learning","Tone learning"])
saveas(fig_fc,[cf,'sup_fig6H.png'])
delete(fig_fc)
