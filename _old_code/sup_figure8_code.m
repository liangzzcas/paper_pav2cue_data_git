%% figure 8A
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

ranksum_p_value_thres = 0.05;

sg_text = "LED diff - Tone diff";

comp_names = ["early peak","dip","late peak"];
cue_omission_tags = ["LEDomi","Toneomi"];

mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];

across_data = struct;
across_data.ct_table = {};
for c_i = 1:2
    across_data.("cue"+c_i+"_rebase_late").diff = [];
    across_data.("cue"+c_i+"_rebase_late").ranksum_p = [];
end

for mouse_name_i = 1:length(mouse_names)
    mouse_name = mouse_names(mouse_name_i);
    this_table = ct_table(ct_table{:,"mouse_name"}==mouse_name,:);
    sig_rois = find(this_table{:,"significance"})';
    across_data.ct_table = cat(1,across_data.ct_table,this_table);
    
    for c_i = 1:2
        % get post-learning
        mu_0 = tas.("cue"+c_i+"late").(mouse_name).mu_location_value(:,:,1);
        single_0 = tas.("cue"+c_i+"late").(mouse_name).single_values;
        sig_0 = tas.("cue"+c_i+"late").(mouse_name).mu_location_value_significance;
        mu_0(~sig_0) = nan;
        % get 20% reward
        mu_1 = tas.("cue"+c_i+cue_omission_tags(c_i)).(mouse_name).mu_location_value(:,:,1);
        single_1 = tas.("cue"+c_i+cue_omission_tags(c_i)).(mouse_name).single_values;
        sig_1 = tas.("cue"+c_i+cue_omission_tags(c_i)).(mouse_name).mu_location_value_significance;
        mu_1(~sig_1) = nan;
        % get 80% reward
        mu_2 = tas.("cue"+c_i+cue_omission_tags(3-c_i)).(mouse_name).mu_location_value(:,:,1);
        single_2 = tas.("cue"+c_i+cue_omission_tags(3-c_i)).(mouse_name).single_values;
        sig_2 = tas.("cue"+c_i+cue_omission_tags(3-c_i)).(mouse_name).mu_location_value_significance;
        mu_2(~sig_2) = nan;
        
        diff_struct_1 = diff_phase_ranksum(...
            mu_1,sig_1,num2cell(single_1,3),...
            mu_0,sig_0,num2cell(single_0,3),this_table{:,"significance"}');

        across_data.("cue"+c_i+"_rebase_late").diff = cat(2,across_data.("cue"+c_i+"_rebase_late").diff,...
            diff_struct_1.diff_mu(sig_rois,:)');
        across_data.("cue"+c_i+"_rebase_late").ranksum_p = cat(2,across_data.("cue"+c_i+"_rebase_late").ranksum_p,...
            diff_struct_1.ranksum_p(sig_rois,:)');
    end
end
across_data.ct_table_sig = across_data.ct_table(logical(across_data.ct_table{:,"significance"}),:);


% plot across data
this_sg_text = sg_text + " increasing only";

fig_fc_late = figure(Position=[100,100,1600,900],InvertHardcopy="off");
tiled_fc_late = tiledlayout(2,3,Parent=fig_fc_late,TileSpacing="tight");
sgtitle(tiled_fc_late,"rebase late "+this_sg_text)

axs_late = gobjects(1,6);
for i = 1:6
    axs_late(i) = nexttile(tiled_fc_late,i);
end

this_axs = axs_late;
across_data_suffix = "_rebase_late";
hold(this_axs,"on")
for i = 1:3 % loop component
    ax1 = this_axs(i); ax2 = this_axs(i+3);
    % take ROIs that are significant either for cue1 or cue2
    val1 = across_data.("cue1"+across_data_suffix).diff(i,:);
    val2 = across_data.("cue2"+across_data_suffix).diff(i,:);
    % zeroing all decreasing ROIs (i.e. late > omission)
    val1(val1<=0) = 0;
    val2(val2<=0) = 0;

    sig_1 = across_data.("cue1"+across_data_suffix).ranksum_p(i,:) <= ranksum_p_value_thres;
    sig_2 = across_data.("cue2"+across_data_suffix).ranksum_p(i,:) <= ranksum_p_value_thres;
    val1(~logical(sig_1)) = nan; val2(~logical(sig_2)) = nan;
    circle_value = common_functions.nan_minus(val1,val2);

    common_functions.scatter_3d(ax1,circle_value,across_data.ct_table_sig,colormapOption='redblue',...
        sorting_tag="abs",viewAngle=[0,90],skip_nan=1,setuserdata=1,colorful_edge_color=[0.65,0.65,0.65]);
    common_functions.scatter_3d(ax2,circle_value,across_data.ct_table_sig,colormapOption='redblue',...
        sorting_tag="abs",viewAngle=[-90,0],skip_nan=1,setuserdata=1,colorful_edge_color=[0.65,0.65,0.65]);
    
    r_b_w = [sum(circle_value>0,"all","omitmissing"),sum(circle_value<0,"all","omitmissing"),sum(circle_value==0,"all","omitmissing")]/sum(~isnan(circle_value))*100;
    title(ax2,[comp_names(i),sprintf("red/blue/white=%0.2f/%0.2f/%0.2f",r_b_w(1),r_b_w(2),r_b_w(3))]);
end
hold(this_axs,"off")
saveas(fig_fc_late,[cf,'sup_figure8A.png'])
delete(fig_fc_late)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 8B
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);

plot_info_1 = {"G23",52,18;"G17",47,18;};
plot_info_2 = ["cue1LEDomi","cue1Toneomi","cue2Toneomi","cue2LEDomi"];
plot_colors = lines(length(plot_info_2));

fig = figure(Position=[100,100,900,400]);
nr=1;
tiled = tiledlayout(fig,nr,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,nr*size(plot_info_1,1));
for i=1:nr*size(plot_info_1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
lgs=gobjects(1,4);
for pi=1:size(plot_info_1,1)
    this_axes = axes(pi);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    
    plot_x = (1:sr*4)/sr-1;
    for pi2=1:length(plot_info_2)
        tmp = cwa.(plot_info_2(pi2)).(this_info{1}).activity(:,this_info{2}+3,:);
        lgs(pi2)=common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(pi2,:)*0.8,sem_color=plot_colors(pi2,:));
    end
    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    if pi==1
        legend(this_axes(1),lgs,plot_info_2,AutoUpdate="off");
    end
    xlim(this_axes(1),[-1,2]);
end
hold(axes,"off");
saveas(fig,'sup_fig8B.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% functions 
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
