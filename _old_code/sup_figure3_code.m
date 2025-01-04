%% figure 3A
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
    plot_data = lick_data.(mouse_name).single_trial_struct.across_1s.cue2_index_across(:,1:last_pav_day);
    datamu = mean(plot_data,1,"omitmissing");
    datasem = std(plot_data,[],1,"omitmissing")/sqrt(size(plot_data,1));
    errorbar(ax,1:last_pav_day,datamu,datasem);
end
hold(ax,"on")
saveas(fig,[cf,'sup_fig3A.png'],'png')
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 3C
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_time_histogram(cf,["cue2early","cue2late"],'sup_fig3C');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 3D
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_presence_circle(cf,"pav2cue task",["cue2early","cue2late"],'sup_fig3D');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 3E
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
CT_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
plot_colors = lines(2);
main_sr = common_functions.get_main_samplerate();

plot_info = {["G12",39],["G21",27],["G21",9]};
fig = figure(Position=[100,100,1200,400]);
tiled = tiledlayout(fig,1,3,TileSpacing="tight");
axes=gobjects(1,3);
for i=1:3
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:length(plot_info)
    this_ax = axes(pi);
    this_info = plot_info{pi};
    sr = main_sr(this_info(1));
    plot_x = (1:sr*4)/sr-1;

    tmp = cwa.cue2early.(this_info(1)).activity(:,str2num(this_info(2))+3,:);
    common_functions.plot_data_single(this_ax,plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
    tmp = cwa.cue2late.(this_info(1)).activity(:,str2num(this_info(2))+3,:);
    common_functions.plot_data_single(this_ax,plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));
    xline(this_ax,0,'--');yline(this_ax,0);
    xlim(this_ax,[-0.5,2]);
end
hold(axes,"off");
saveas(fig,[cf,'sup_fig3E.fig'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 3F
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue2late","cue2early"]},mouse_names,'sup_fig3F');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 3G
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

% phase_names = ["cue1late"]; % use one column in one row to get timing
% phase_names = ["cue1early","cue1late"]; % use two column in one row to get difference value (2nd - 1st)
phase_names = {["cue2late"];};
mouse_names = {["G12","G15","G17","G19","G21","G22","G23","G24"]};

n_phase_name = size(phase_names,1);
n_mouse_name = size(mouse_names,1);

AX_TITLES = ["early peak","dip","late peak"];
histo_cmap_tags = ["cue1","cue2","rew1","rew2"];

for n_mouse_name_i = 1:n_mouse_name
    mouse_name = mouse_names{n_mouse_name_i};
    ct_table_out = ct_table(contains(string(cell2mat(ct_table{:,"mouse_name"})),mouse_name),:);
    ct_table_out{:,"ROI"} = (1:size(ct_table_out,1))';
    bad_fiber = ct_table_out{:,"significance"}~=1;
    bad_fiber = bad_fiber';
    for n_phase_name_i = 1:n_phase_name
        phase_name = phase_names{n_phase_name_i};
        data_by_phase = common_functions.get_by_mouse_phase(tas,mouse_name,phase_name);
        
        if length(phase_name) == 1
            diff_phase_flag = false;
            out = data_by_phase.(phase_name(1)).this_tas_single;
            % take invert value for dip
            out(2,:,1) = -out(2,:,1);
            out(:,:,2) = out(:,:,2) - 1; % because I took 1s before onset
            out_dim = out;
            this_percent_flag = 0; % no need to use percentage when plotting amplitude
            plot_color = lines(1);
            this_colormap = 'parula';
            cmapBounds_comps = {[0,0.5],[],[]};
        elseif length(phase_name) == 2
            diff_phase_flag = true;
            out1 = data_by_phase.(phase_name(1)).this_tas_single;
            out2 = data_by_phase.(phase_name(2)).this_tas_single;
            % take invert value for dip
            out1(2,:,1) = -out1(2,:,1);
            out2(2,:,1) = -out2(2,:,1);
            out1(:,:,2) = out1(:,:,2) - 1; % because I took 1s before onset
            out2(:,:,2) = out2(:,:,2) - 1; % because I took 1s before onset
            out = random_functions.nan_minus(out2,out1);
            out_dim = out;
            this_percent_flag = percent_flag;
            if this_percent_flag
                out_dim = sign(out_dim);
                out_dim(out_dim==-1) = 0;
            end
            plot_color = lines(2);
            plot_color = plot_color(2,:);
            this_colormap = 'redblue';
            cmapBounds_comps = {[-1,1],[-1,1],[-1,1]}; % haven't confirmed this but difference shouldn't be larger than 0.5 I guess
        end
        

        % - circle plot
        circle_fig = figure(Position = [100,100,1600,1100]);
        circle_tiled = tiledlayout(circle_fig,2,3,TileSpacing="tight");
        sgtitle(circle_tiled, "time of components " + strjoin([mouse_name,phase_name]," "))
        axs = gobjects(1,6);
        % for AX_TITLES_i = 1:3
        for AX_TITLES_i = 1:1
            axs(AX_TITLES_i) = nexttile(circle_tiled,AX_TITLES_i);
            axs(AX_TITLES_i+3) = nexttile(circle_tiled,AX_TITLES_i+3);
            ax = axs(AX_TITLES_i); ax1 = axs(AX_TITLES_i+3);
            hold([ax,ax1],"on")
            plot_data = out(AX_TITLES_i,:,2);
            common_functions.scatter_3d(ax,plot_data,ct_table_out,cmapBounds=cmapBounds_comps{AX_TITLES_i},viewAngle=[0,90],colormapOption=this_colormap, skipBubble=bad_fiber,...
                sorting_tag = "abs", setuserdata=0);
            common_functions.scatter_3d(ax1,plot_data,ct_table_out,cmapBounds=cmapBounds_comps{AX_TITLES_i},viewAngle=[-90,0],colormapOption=this_colormap, skipBubble=bad_fiber,...
                sorting_tag = "abs", setuserdata=0);
            hold([ax,ax1],"off")
            title(ax,AX_TITLES(AX_TITLES_i))
        end

        this_colormap = parula(500);
        this_colormap = cat(1,[1,1,1],this_colormap);
        for AX_TITLES_i = 2
            axs(AX_TITLES_i) = nexttile(circle_tiled,AX_TITLES_i);
            axs(AX_TITLES_i+3) = nexttile(circle_tiled,AX_TITLES_i+3);
            ax = axs(AX_TITLES_i); ax1 = axs(AX_TITLES_i+3);
            hold([ax,ax1],"on")
            plot_data = out(1,:,2);
            [~,~,~,hotzone_hor] = plot_functions.interp_scatter_3d(ax,plot_data,ct_table_out,cmapBounds=[0,0.6],viewplane="hor",colormapOption=this_colormap, skipBubble=bad_fiber,hotzone_prctile_thres=70);
            [~,~,~,hotzone_sag] = plot_functions.interp_scatter_3d(ax1,plot_data,ct_table_out,cmapBounds=[0,0.6],viewplane="sag",colormapOption=this_colormap, skipBubble=bad_fiber,hotzone_prctile_thres=70);
            hold([ax,ax1],"off")
            title(ax,AX_TITLES(1))

            %   save data to write to output file
            hor_region = hotzone_hor{1};
            sag_region = hotzone_sag{1};
            tmp = ct_table_out{~bad_fiber&~isnan(plot_data),["fiber_bottom_AP","fiber_bottom_ML","fiber_bottom_DV"]};
            output = common_functions.fiber_inpolygon(tmp,hor_region,sag_region,{},{});
            fiber_bit = output{1}&output{2};

            file_to_save = struct;
            file_to_save.description = "contour = hor/sag; fiber = only one list, hotzone_prctile_thres = 70";
            file_to_save.contour = {hotzone_hor;hotzone_sag};
            file_to_save.fiber = tmp(fiber_bit,:);
            random_functions.save_to_struct(file_to_save,"cue2late_pk_timing",target_tag="hotzone_contour");
        end
        % axis(axs,"vis3d")
        axis(axs([1,4]),"vis3d")
        % saveas(circle_fig,[cf,'sup_fig3G',char(strjoin([mouse_name,phase_name],"_")),'.png']);
        saveas(circle_fig,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot\sup_fig3G',char(strjoin([mouse_name,phase_name],"_")),'.fig']);
        saveas(circle_fig,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot\sup_fig3G',char(strjoin([mouse_name,phase_name],"_")),'.png']);
        delete(circle_fig)
    end
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 3H
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];

suffix_text = "Tone diff vs LED diff";
comp_names = ["pk","dp","re"];
ranksum_p_value_thres = 0.05;

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
        this_mu_tone_e = tas.("cue"+c_i+"early").(mouse_name).mu_location_value(:,:,1);
        this_mu_tone_single_e = tas.("cue"+c_i+"early").(mouse_name).single_values;
        this_sig_tone_e = tas.("cue"+c_i+"early").(mouse_name).mu_location_value_significance;
        this_mu_tone_e = this_mu_tone_e .* this_sig_tone_e;

        this_mu_tone_l = tas.("cue"+c_i+"late").(mouse_name).mu_location_value(:,:,1);
        this_mu_tone_single_l = tas.("cue"+c_i+"late").(mouse_name).single_values;
        this_sig_tone_l = tas.("cue"+c_i+"late").(mouse_name).mu_location_value_significance;
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
saveas(fig_fc,[cf,'sup_fig3H.png'])
delete(fig_fc)
