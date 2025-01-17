%% figure 5AC
close all;clear;clc;
cf = [pwd,'\'];
cwa = load([cf,'processed_and_organized_data\ITI_licking_components_window_activity_filtered.mat']);

plot_info_1 = {"G23",36,18;"G23",19,18};
plot_imagesc_clim = {[-1,1]*0.02;[-1,1]*0.045;};
plot_info_2 = ["early","late"];
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
    
    plot_x = (1:sr*2)/sr-1;
    if pi~=3
        tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
        common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
        imagesc(this_axes(2),plot_x,1:nt1,permute(tmp,[3,1,2]),plot_imagesc_clim{pi});
        this_axes(2).YDir = "reverse";
        colormap(this_axes(2),common_functions.redblue());
        colorbar(this_axes(2));
        tmp = cwa.(plot_info_2(2)).(this_info{1}).activity(:,this_info{2}+3,:); nt2 = size(tmp,3);
        common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));
        imagesc(this_axes(3),plot_x,1:nt2,permute(tmp,[3,1,2]),plot_imagesc_clim{pi});
        this_axes(3).YDir = "reverse";
        colormap(this_axes(3),common_functions.redblue());
        colorbar(this_axes(3));
    else
        tmp = cwa.(plot_info_2(1)).(this_info{1}).activity(:,this_info{2}+3,:); nt1 = size(tmp,3);
        tmp = tmp/max(tmp,[],"all","omitmissing");
        common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(1,:)*0.8,sem_color=plot_colors(1,:));
        imagesc(this_axes(2),plot_x,1:nt1,permute(tmp,[3,1,2]),plot_imagesc_clim{pi});
        this_axes(2).YDir = "reverse";
        colorbar(this_axes(2));
        tmp = cwa.(plot_info_2(2)).(this_info{1}).activity(:,this_info{2}+3,:); nt2 = size(tmp,3);
        tmp = tmp/max(tmp,[],"all","omitmissing");
        common_functions.plot_data_single(this_axes(1),plot_x,tmp,plot_color=plot_colors(2,:)*0.8,sem_color=plot_colors(2,:));
        imagesc(this_axes(3),plot_x,1:nt2,permute(tmp,[3,1,2]),plot_imagesc_clim{pi});
        this_axes(3).YDir = "reverse";
        colorbar(this_axes(3));
    end

    xline(this_axes(1),0,'--');yline(this_axes(1),0);
    xlim(this_axes(1),[-1,1]);
    for ax=this_axes([2,3])
        xline(ax,0,'--');
        xlim(ax,plot_x([1,end]));
    end
    ylim(this_axes(2),[1,nt1]);
    ylim(this_axes(3),[1,nt2]);
    if pi==3
        ylim(this_axes(1),[0,1]);
    end
end
hold(axes,"off");
saveas(fig,'fig5A_C.png');
delete(fig)


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 5BDEF
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.smoothed_spatial_map(cf,"tas_itilick",{["ITIlate","ITIearly"]},mouse_names,'fig5EF',tabular_sheet_name="fig5BD");


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%--------------------------------------------------------------------------%
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%


%% figure 5G
close all;clear;clc;
cf = [pwd,'\'];
mouse_names = ["G12","G15","G17","G19","G21","G22","G23","G24"];
common_plot_functions.modality_between_context(cf,"extinction_vs_lick",mouse_names,'fig5G')


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
