%% figure 6c
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["422","423"];
fnames = ["unpred"];
ftags = ["unpred"];
this_sr = 18;
across_data = load([cf,'processed_and_organized_data\TeLC_data.mat']);

plot_data = struct;
plot_data.unpred470 = [];
plot_data.unpred405 = [];
for mi=1:2
    mname = mouse_names(mi);
    mnametag = "m"+mname;
    plot_data.unpred470 = cat(2,plot_data.unpred470,mean(across_data.(mnametag).unpred470,3,"omitmissing"));
    plot_data.unpred405 = cat(2,plot_data.unpred405,mean(across_data.(mnametag).unpred405,3,"omitmissing"));
end

plot_x = (-this_sr+1:2*this_sr)/this_sr;

fig = figure(position=[100,100,800,600]);
tiled = tiledlayout(1,1,TileSpacing="tight");
axs = [nexttile(tiled,1)];
hold(axs,"on")
for ai=1:1
    tmp = plot_data.(fnames(ai)+"405");
    rebase = mean(tmp(1:this_sr,:,:),[1,3],"omitmissing");
    tmp = tmp - rebase;
    p2 = plot_functions.plot_data_single(axs(ai),plot_x,tmp,plot_color=[1,1,1]*0.8*0.4,sem_color=[1,1,1]*0.4);

    tmp = plot_data.(fnames(ai)+"470");
    rebase = mean(tmp(1:this_sr,:,:),[1,3],"omitmissing");
    tmp = tmp - rebase;
    p1 = plot_functions.plot_data_single(axs(ai),plot_x,tmp,plot_color=[1,0.7,0.5]*0.8,sem_color=[1,0.7,0.5]);
end
hold(axs,"off")
for axi=1:1
    ax = axs(axi);
    xline(ax,0);
    yline(ax,0);
    title(ax,fnames(axi));
    legend(ax,[p1,p2],["470","405"]);
    xlabel("Time since cue onset (s)")
    ylabel("\DeltaF/F")
    ylim(ax,[-1,1]*0.02)
end
saveas(fig,'fig6C.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

