%% figure 2BC
close all;clear;clc;
cf = [pwd,'\'];
this_ta = load([cf,'processed_and_organized_data\event_aligned_highpass03.mat']).cwa_raw.G23;
session_info = common_functions.get_training_info();
this_pavday = session_info{[session_info{:,1}]=="G23",2};

imagesc_ta = nan(72,this_pavday);
for di = 1:this_pavday
    imagesc_ta(:,di) = this_ta.("file"+di).cueOn.cue1.mu(:,1);
end
imagesc_ta_1 = permute(this_ta.("file"+this_pavday).cueOn.cue1.activity(:,1,:),[1,3,2]);

fig = figure(Position=[100,100,800,600]);
ax=axes(fig);
imagesc(ax,(1:72)/18-1,1:size(imagesc_ta,2),imagesc_ta'/max(imagesc_ta,[],"all","omitmissing"));
colorbar(ax);
xline(ax,0,LineWidth=1.5,Color=[1,1,1]);
xlabel(ax,"Time (s)"); ylabel(ax, "Days")
% saveas(fig,[cf,'fig2B.png'],'png')
% delete(fig)

fig = figure(Position=[100,100,800,600]);
ax=axes(fig);
imagesc(ax,(1:72)/18-1,1:size(imagesc_ta_1,2),imagesc_ta_1'/max(imagesc_ta_1,[],"all","omitmissing"));
colorbar(ax);
xline(ax,0,LineWidth=1.5,Color=[1,1,1]);
xlabel(ax,"Time (s)"); ylabel(ax, "Trial #")
% saveas(fig,[cf,'fig2C.png'],'png')
% delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 2D
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
    plot_data = lick_data.(mouse_name).single_trial_struct.across_1s.cue1_index_across(:,1:last_pav_day);
    datamu = mean(plot_data,1,"omitmissing");
    datasem = std(plot_data,[],1,"omitmissing")/sqrt(size(plot_data,1));
    errorbar(ax,1:last_pav_day,datamu,datasem);
end
hold(ax,"on")
saveas(fig,[cf,'fig2D.png'],'png')
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 2F
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\components_window_activity_filtered_rebase_rew_consump.mat']);

plot_info_1 = {"G12",39,30;"G21",27,18;"G19",32,18;};
plot_imagesc_clim = {0.05;0.06;0.03;};
plot_info_2 = ["cue1early","cue1late"];
plot_colors = lines(2);

fig = figure(Position=[100,100,1200,1000]);
tiled = tiledlayout(fig,3,size(plot_info_1,1),TileSpacing="tight");
axes=gobjects(1,3*size(plot_info_1,1));
for i=1:3*size(plot_info_1)
    axes(i) = nexttile(tiled,i);
end
hold(axes,"on");
for pi=1:size(plot_info_1,1)
    this_axes = axes([pi,pi+size(plot_info_1,1),pi+2*size(plot_info_1,1)]);
    this_info = plot_info_1(pi,:);
    sr = this_info{3};
    
    plot_x = (1:sr*4)/sr-1;
    tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
    imagesc(this_axes(2),plot_x,1:nt1,permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
    this_axes(2).YDir = "reverse";
    colormap(this_axes(2),common_functions.redblue());
    colorbar(this_axes(2));
    tmp = cwa.(plot_info_2(2)).(this_info{1}).activity(:,this_info{2}+3,:); nt2 = size(tmp,3);
    common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));
    imagesc(this_axes(3),plot_x,1:nt2,permute(tmp,[3,1,2]),[-1,1]*plot_imagesc_clim{pi});
    this_axes(3).YDir = "reverse";
    colormap(this_axes(3),common_functions.redblue());
    colorbar(this_axes(3));
    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-0.5,2]);
    for ax=this_axes([2,3])
        xline(ax,0,'--');
        xlim(ax,plot_x([1,end]));
    end
    ylim(this_axes(2),[1,nt1]);
    ylim(this_axes(3),[1,nt2]);
end
hold(axes,"off");
saveas(fig,'fig2F.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 2G
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_time_histogram(cf,["cue1early","cue1late"],'fig2G');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 2H
close all;clear;clc;
cf = [pwd,'\'];
common_plot_functions.component_presence_circle(cf,"pav2cue task",["cue1early","cue1late"],'fig2H');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% fig 2I
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.phase_diff_circlemap(cf,"tas_normal",{["cue1late","cue1early"]},mouse_names,'fig2I');


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%