%% sup figure 2A
% loading % setting
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);
CT_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);
plot_info = {{["G17",35],["G17",39],["G17",32]};...
    {["G17",37],["G17",41],["G19",12],["G24",35]};...
    {["G12",51],["G12",21],["G21",26],["G17",40],["G17",48]};...
    {["G17",33],["G19",49],["G22",2],["G23",48],["G12",39]};...
    {["G21",28],["G12",19],["G22",65],["G17",28],["G22",42],["G15",5]};...
    {["G12",41],["G17",6]};...
    {["G22",27]};...
    {["G12",31],["G24",48]};...
    {["G17",23]};...
    };
main_sr = common_functions.get_main_samplerate();

% sort plot color by AP gradient
tmp = horzcat([plot_info{:}]);
tmp = cellfun(@(x) CT_table{CT_table{:,"mouse_name"}==x(1)&CT_table{:,"ROI_original"}==str2num(x(2)),"fiber_bottom_AP"},tmp);
[~,I] = sort(tmp,"descend");
plot_colors_id = nan(size(I));
for ii = 1:length(I)
    plot_colors_id(I(ii)) = ii;
end

plot_colors = winter(29); pci = 1;
fig = figure(Position=[100,100,1200,1000]);
tiled = tiledlayout(fig,3,3,TileSpacing="tight");
axes=gobjects(1,3*3);
for i=1:3*3
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:9
    this_ax = axes(pi);
    this_infos = plot_info{pi};
    for ri = 1:size(this_infos,2)
        this_info = this_infos{ri};
        tmp = cwa.cue1late.(this_info(1)).activity(:,str2num(this_info(2))+3,:);
        sr = main_sr(this_info(1));
        plot_x = (1:sr*4)/sr-1;
        common_functions.plot_data_single(this_ax,plot_x,tmp,plot_color=plot_colors(plot_colors_id(pci),:)*0.8,sem_color=plot_colors(plot_colors_id(pci),:)); pci=pci+1;
    end
    xline(this_ax,0,'--');yline(this_ax,0);
    xlim(this_ax,[-1,2]);
end
hold(axes,"off");
saveas(fig,[cf,'sup_fig2A.fig'])
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% sup figure 2B
% loading % setting
close all;clear;clc;
cf = [pwd,'\'];
tas = load([cf,'processed_and_organized_data\tri_avg_single_filtered_rebase_rew_consump.mat']);
ct_table = readtable([cf,'raw_data\CT_across_GXX_mice.xlsx']);

% phase_names = ["cue1late"]; % use one column in one row to get timing
% phase_names = ["cue1early","cue1late"]; % use two column in one row to get difference value (2nd - 1st)
phase_names = {["cue1late"];};
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
            cmapBounds_comps = {[0,0.4],[],[]};
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

        % interpolated circle plot (need modify)
        this_colormap = parula(500);
        this_colormap = cat(1,[1,1,1],this_colormap);
        for AX_TITLES_i = 2
            axs(AX_TITLES_i) = nexttile(circle_tiled,AX_TITLES_i);
            axs(AX_TITLES_i+3) = nexttile(circle_tiled,AX_TITLES_i+3);
            ax = axs(AX_TITLES_i); ax1 = axs(AX_TITLES_i+3);
            hold([ax,ax1],"on")
            plot_data = out(1,:,2);
            [~,~,~,hotzone_hor] = plot_functions.interp_scatter_3d(ax,plot_data,ct_table_out,cmapBounds=[0,0.4],viewplane="hor",colormapOption=this_colormap, skipBubble=bad_fiber,hotzone_prctile_thres=70);
            [~,~,~,hotzone_sag] = plot_functions.interp_scatter_3d(ax1,plot_data,ct_table_out,cmapBounds=[0,0.4],viewplane="sag",colormapOption=this_colormap, skipBubble=bad_fiber,hotzone_prctile_thres=70);
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
            random_functions.save_to_struct(file_to_save,"cue1late_pk_timing",target_tag="hotzone_contour");
        end

        % axis(axs,"vis3d")
        axis(axs([1,4]),"vis3d")
        axis(axs([2,5]),"vis3d")
        % saveas(circle_fig,[cf,'sup_fig2B',char(strjoin([mouse_name,phase_name],"_")),'.png']);
        saveas(circle_fig,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot\sup_fig2B',char(strjoin([mouse_name,phase_name],"_")),'.fig']);
        saveas(circle_fig,['F:\Safa_Processed\#paper_figure\#update_review\_revision_plot\sup_fig2B',char(strjoin([mouse_name,phase_name],"_")),'.png']);
        delete(circle_fig)
    end
end


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%