classdef common_functions
    methods(Static)
        function output = ball2xy(behav,varargin)
            % 
            % By Dr. Maianh Vu @ Howelab
            %
            % output = ball2xy(behav)
            %
            % This function takes NIDAQ behavioral output, and transforms the
            % locomotion data from the ball (pitch, yaw, roll) into 2D locomotion 
            % parameters. It returns linear (XY)velocity (approximately m/s), 
            % angular velocity (approximately radians/s), as well X- and Y-
            % position coordinates (approximately in m), assuming the mouse starts 
            % at [0 0].
            %
            % It takes as input the behavioral data struct or filepath.
            %
            % It takes as optional inputs:
            %   voltage2meters:     the scaling factor from volts/s to meters/s
            %                       Basically the idea here is to calculate the 
            %                       ratio of the animal's velocity over 3.5 m/s 
            %                       (Mai-Anh found in a paper the the approximate 
            %                       maximum mouse running velocity is 3.5 m/s). Then 
            %                       multiply this ratio by 3.3V to determine how much 
            %                       voltage to send out. This means that 3.3V will be 
            %                       max velocity of 3.5 m/s, 0V will be 0 m/s
            %   ballDiameter:       the diameter of the ball in m 
            %                           (default = 8in = 0.2032m)           
            %   yawFromBoth:        whether to use both sensors to calculate yaw
            %                           (default = 0 = no)
            %   axelBall:           set this to 1 if you're using the axel ball to 
            %                       automatically ignore any yaw or roll input
            %                           (default = 0)
            %
            % It returns a struct with the following fields:
            %   roll_velocity:      roll velocity (m/s)
            %   pitch_velocity:     pitch velocity (m/s)
            %   x_velocity:         velocity along cartesian x-axis (m/s)
            %   y_velocity:         velocity along cartesian y-axis (m/s)
            %   linear_velocity:    linear velocity magnitude (m/s), sqrt(x^2+y^2)
            %   angular_velocity:   angular velocity (radians/s)
            %   angular_velocity_ms:angular velocity (m/s)
            %
            % Mai-Anh 3/25/2021
            % updated 6/28/2021 to output velocity_magnitude
            % updated 6/29/2021: velocity_magnitude is SAME as linear_velocity haha, but
            %       left both in the output for convenience
            % updated 2/8/2022: optional input for axel ball to automatically ignore yaw and roll
            % updated 6/14/2022: fix 2 calculations, oops:
            %           1) angular velocity now in radians (used to be in #rotations)
            %           2) sign calculation was wrong before and was amplifying
            %           positive magnitudes
            % updated 8/18/2022: don't need to divide sample duration because units are
            %           already m/s, and so the split/bin averaging is already m/s
            % updated 11/9/2022: also outputs angular_velocity in m/s 
            % updated 2/3/2023: re: 8/18/2022 turns out that fix was not correct. it
            %   was in fact correct to divide by sample duration. i tested this by
            %   verifying distance traveled on Brenna's linear track task, where the
            %   distance is carefully calibrated. there does need to be division by
            %   sample rate
            % updated 2/4/2023: I fixed the scaling, and verified it against Brenna's 
            %   linear track with its carefully calibrated distance, and I've also 
            %   tested this against multiple sampling rates. I wrote the code in the 
            %   raspberry pis to take the optical sensor output, which is dx and dy in 
            %   pixels. Then I convert the distance in pixels to distance in meters, 
            %   and then to a velocity by dividing by sample time. Then to determine
            %   a voltage to send out. I divide this by a theoretical max mouse 
            %   velocity of 3.5 m/s and multiply that by 3.3V. In other words, a 
            %   velocity of 3.5m/s would be 3.3V. To convert the voltage from the 
            %   sensors back to m/s, divide by 3.3V and multiply by 3.5 m/s.
        
            %%%%%%%%%%%%%%%%%%%%%%
            %%%% PARSE INPUTS %%%%
            %%%%%%%%%%%%%%%%%%%%%%
            
            ip = inputParser;
            % voltage2meters breakdown: basically the idea here is to calculate the 
            % ratio of the animal's velocity over 3.5 m/s (Mai-Anh found in a paper
            % the the approximate maximum mouse running velocity is 3.5 m/s). Then 
            % multiply this ratio by 3.3V to determine how much voltage to send out.
            % This means that 3.3V will be max velocity of 3.5 m/s, 0V will be 0 m/s
            ip.addParameter('voltage2meters',1/3.3*3.5); % 3.3V = max = 3.5 m/s
            ip.addParameter('ballDiameter',8*0.0254); % ball diameter is in meters: this is based on an 8in ball
            ip.addParameter('yawFromBoth',0); % set to 1 to use an average of sensors 1 and 2. default is 0.
            ip.addParameter('axelBall',0); % set to 1 if using axel ball, so yaw and roll go to 0. default is 0.
            
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            
            %%% load it if it isn't already a struct
            if ~isstruct(behav)
                if ~endsWith(behav,'.mat')
                    behav = [behav '.mat'];
                end
                behav = load(behav);
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% COMBINE MAGNITUDE & SIGN %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Calculating pitch, roll, and yaw: take magnitude and multiply by sign
            % Note: sign comes into nidaq as 0 or 3.33V (analog pin). So here, the 0
            % point is 1.65, which really only matters in split/bin data (previously,
            % the Boolean was >1, which is what it is in the VirMEn tasks and
            % daqSessionSave, but at the 1-2kHz sampling rate that's applicable to
            % that, the boolean is fine because the signal will be 0 or 3.33V).
            pitchRaw = behav.ballSensor1_y.*(2*double(behav.ballSensor1_ysign>1.65)-1);
            if axelBall == 0    
                rollRaw = behav.ballSensor2_y.*(2*double(behav.ballSensor2_ysign>1.65)-1);
                yawRaw = behav.ballSensor1_x.*(2*double(behav.ballSensor1_xsign>1.65)-1);
                % if we're calculating yaw from both sensors, take the mean
                if yawFromBoth == 1 
                    yawRaw2 = behav.ballSensor2_x.*(2*double(behav.ballSensor2_xsign>1.65)-1);
                    yawRaw = mean([yawRaw  yawRaw2],2);
                end
            else
                rollRaw = zeros(size(pitchRaw));
                yawRaw = zeros(size(rollRaw));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% UNIT CONVERSIONS %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % volts/s back to meters/s
            pitchMeters = pitchRaw * voltage2meters;
            rollMeters = rollRaw * voltage2meters;
            yawMeters = yawRaw * voltage2meters;
            
            % meters/s to radians/s
            yawRadians = yawMeters/(pi*ballDiameter)*2*pi; % meters to #rotations
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% VELOCITY CALCULATIONS %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % translate the XY velocity and angular velocity to position
            headingAngle = cumsum(yawRadians);
            velLinear = nan(size(pitchMeters,1),2); 
            
            % do a loop... there's probably a more efficient way, but whatever
            % distLinear = distance traveled during a frame
            for i = 1:size(pitchMeters,1)
                if i == 1
                    lastAng = 0;
                else
                    lastAng = headingAngle(i-1);
                end
                velLinear(i,:) = [cos(lastAng) -sin(lastAng); sin(lastAng) cos(lastAng)]*[rollMeters(i) pitchMeters(i)]';
            end
            
            
            %%%%%%%%%%%%%%%%
            %%%% OUTPUT %%%%
            %%%%%%%%%%%%%%%%
            
            % various velocities
            output.roll_velocity = rollMeters;
            output.pitch_velocity = pitchMeters;
            output.x_velocity = velLinear(:,1);
            output.y_velocity = velLinear(:,2);
            output.linear_velocity = (velLinear(:,1).^2 + velLinear(:,2).^2).^.5; % linear velocity
            output.angular_velocity = yawRadians; % angular velocity in radians
            output.angular_velocity_ms = yawMeters; % angular velocity in m/s calculated how we calculate pitch...
        end
        
        function [gs,cmapBounds_out] = scatter_3d(plot_ax,values,fib_loc, varargin)
            % By Liangzhu Zhang @ Howelab
            % 3d circle plot implemented in scatter3 (inputs have similar style as Maianh bubbleplot3d)
            setuserdata  = 1;
            skipBubble = zeros(size(values));
            cmapBounds = [];
            outlineID = nan(size(values));
            outlineColor = [0 0 0];
            outlineWidth = 1;
            skip_nan = 0;
            nan_style = {'MarkerFaceColor','none','MarkerFaceAlpha','0','MarkerEdgeColor',[0.8,0.8,0.8],'MarkerEdgeAlpha','1',...
                'Marker','x','LineWidth',1};
            bubbleSize = 140;
            colormapOption = 'parula';
            colormapBins = 256;
            viewAngle = [-90,0];
            DVlim = [2 6];
            APlim = [-2.1 1.6];
            MLlim = [0.5 4.2];
            dynamic_cmap_bound_flag = 0;
            sorting_tag = 'location';
            colorful_edge_color = [0.2,0.2,0.2];
            grey_bubble = zeros(size(values));
            grey_bubble_color = [0.65,0.65,0.65];
            grey_bubble_pvalue = [];
            special_color_flag = zeros(size(values));
            special_color_value = [];
            roi_label = 0 ;
            roi_label_color = [0.4667,0.6745,0.1882];

            % circle plot based on input data and CT location
            ip = inputParser;
            %%% basic settings
            ip.addParameter('setuserdata',1) % setup interactive user datatips
            ip.addParameter('skipBubble',zeros(size(values))) % do not show this bubbles (bit array)
            %%% bubble appearance
            ip.addParameter('cmapBounds',[]); % manually set the color bounds (note RedBlue will automatically make it symmetric)
            ip.addParameter('outlineID',nan(size(values))); % the indices that go with outline color or width (i.e., if you have 2 colors, you can specify whether the fiber should be color 1 or 2)
            ip.addParameter('outlineColor',[0 0 0]);% outline color: one row per color. row indices correspond with outlineID
                % if outlineColor contains nan, make it black outline only
            ip.addParameter('outlineWidth',1);% outline width: you can have multiple widths corresponding to outlineID if you want
            ip.addParameter('skip_nan',0); % if skipnan, do not plot nan
            ip.addParameter('nan_style',{'MarkerFaceColor','none','MarkerFaceAlpha','0','MarkerEdgeColor','k','MarkerEdgeAlpha','1',...
                'Marker','x','LineWidth',1.5}); % style of plotting nan
            %%% bubble size options
            ip.addParameter('bubbleSize',140); % size of bubbles
            %%% color map and colorbar options
            ip.addParameter('colormapOption','parula'); % colormap (either the string to a map, or a matrix)
            ip.addParameter('colormapBins',256); % the number of color bins for the colormap
            %%% axes options
            ip.addParameter('viewAngle',[-90 0]); % default view is sagittal
            ip.addParameter('DVlim',[]); % DV limits: if blank, it will be determined automatically. I like [-5.5, -2].
            ip.addParameter('APlim',[]); % AP limits: if blank, it will be determined automatically. I like [-2 -2].
            ip.addParameter('MLlim',[]); % ML limits: if blank, it will be determined automatically. I like [-.5 3.75] for R implants
            
            %%% Zack modulation %%%
            ip.addParameter('dynamic_cmap_bound_flag',0);

            ip.addParameter('sorting_tag','location');
            ip.addParameter('colorful_edge_color','none');
            ip.addParameter('grey_bubble',zeros(size(values))); % bubbles to grey out (bit array)
            ip.addParameter('grey_bubble_color',[0.65,0.65,0.65]);
            ip.addParameter('grey_bubble_pvalue',[]);
        
            % each flag matches with a color rgb row in special_color_value (work in the same way as outlineColor but for facecolor)
            ip.addParameter('special_color_flag',zeros(size(values))); 
            ip.addParameter('special_color_value',[]);
            
            ip.addParameter('roi_label',0);
            ip.addParameter('roi_label_color',[0.4667,0.6745,0.1882]);
            %%% parser
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end

            skipBubble = logical(skipBubble);
            grey_bubble = logical(grey_bubble);
        
            if all(viewAngle == [0,90])
                view_plane = "hor";
            elseif all(viewAngle == [-90,0])
                view_plane = "sag";
            else
                error("Other angle are to be implemented.")
            end

            no_signal_insig_alpha = 0.4;
        
            assert(max(special_color_flag)==0 || max(special_color_flag)==size(special_color_value,1),"Special color flag and value size must match.")
            if numel(outlineWidth)<size(outlineColor,1)
                outlineWidth = repmat(outlineWidth,1,size(outlineColor,1));
            end
        
            % automatically set cmapBounds
            if isempty(cmapBounds)
                value_mu = mean(values,'omitnan'); value_std = std(values,'omitnan');
                if isempty(skipBubble)
                    bound_range = true(size(values));
                else
                    bound_range = ~skipBubble;
                end
                if dynamic_cmap_bound_flag
                    lo_bound = max([min(values(bound_range)),value_mu-2.5*value_std,prctile(values,3,'all')]);
                    up_bound = min([max(values(bound_range)),value_mu+2.5*value_std,prctile(values,97,'all')]);
                else
                    lo_bound = min(values(bound_range));
                    up_bound = max(values(bound_range));
                end
                cmapBounds = [lo_bound, up_bound];
            end
            % preprocess data
            if istable(fib_loc)
                DV = fib_loc.fiber_bottom_DV;
                AP = fib_loc.fiber_bottom_AP;
                ML = fib_loc.fiber_bottom_ML;
                roiNums = 1:length(fib_loc.ROI);
            elseif isfield(fib_loc,'rc') % this is a calibration file
                DV = fib_loc.rc(:,4);
                AP = fib_loc.rc(:,2);
                ML = fib_loc.rc(:,3);
                roiNums = 1:length(fib_loc.rc(:,1));
            else % it's a fiber table
                if isfield(fib_loc,'table')
                    fib_loc = fib_loc.table;
                end
                DV = fib_loc.fiber_bottom_DV;
                AP = fib_loc.fiber_bottom_AP;
                ML = fib_loc.fiber_bottom_ML;
                roiNums = 1:length(fib_loc.ROI);
            end
        
            % colormap
            if isnumeric(colormapOption) % if a colormap matrix has been supplied
                cmapping = linspace(min(cmapBounds),max(cmapBounds),size(colormapOption,1));
                cmapcolors = colormapOption;
            else
                colormapOption = char(colormapOption);
                eval(['cmapcolors = ' colormapOption '(' num2str(colormapBins) ');']);
                if strcmp(colormapOption,'redblue') % symmetric colormap for redblue
                    cmapping = linspace(-max(abs(cmapBounds)), max(abs(cmapBounds)),colormapBins);
                else
                    cmapping = linspace(min(cmapBounds),max(cmapBounds),colormapBins);
                end
            end
            cmapBounds_out = cmapping([1,end]);
            % set up scatter3 parameters
            % x,y,z,bubbleSize,color,'filled',marker,facealpha
            roiNums = roiNums(:);
            x = ML(roiNums)';
            y = AP(roiNums)';
            depth = DV(roiNums)';
        
            not_skipped = ~skipBubble;
            
            unique_outline = unique(outlineID);
            unique_outline = unique_outline(~isnan(unique_outline));
            
            colorful_wo_outline = not_skipped & ~isnan(values)  & ~grey_bubble & isnan(outlineID);
            n_groups = 3 + length(unique_outline); % nan (no component), grey (insig), normal, outlined
            group_bit = nan(n_groups,length(roiNums));
            group_bit(1,:) = isnan(values) & not_skipped;
            group_bit(2,:) = grey_bubble & not_skipped;
            group_bit(3,:) = colorful_wo_outline;
            for i = 4:n_groups
                this_outlineID = outlineID == unique_outline(i-3);
                group_bit(i,:) = this_outlineID & not_skipped;
            end
            group_bit = logical(group_bit);

            if skip_nan
                group_bit(1,:) = false(size(group_bit(1,:)));
            end
            % get colors for colorful ones
            colorful_colors = nan(length(values),3);
            for i = 1:length(values)
                thisVal = values(i);
                if isnan(thisVal)
                    continue
                elseif thisVal<min(cmapping)
                    thisColor = cmapcolors(1,:);
                elseif thisVal>=max(cmapping)
                    thisColor = cmapcolors(end,:);
                else
                    thisColor = cmapcolors(find(cmapping<=thisVal,1,'last'),:);
                end
                colorful_colors(i,:) = thisColor;
            end
            this_bit1 = group_bit(1,:);
            this_bit2 = group_bit(2,:);
            this_bit3 = group_bit(3,:);
            % change rendering based on sorting
            amp_to_sort = values(this_bit3);
            new_dim = linspace(0,1,length(amp_to_sort));
            
            switch string(sorting_tag)
                case "location"
                    I = 1:length(amp_to_sort);
                case "ascend"
                    [~,I] = sort(amp_to_sort,"ascend");
                case "descend"
                    [~,I] = sort(amp_to_sort,"descend");
                case "abs"
                    [~,I] = sort(abs(amp_to_sort),"descend");
            end
            if ~strcmp(sorting_tag,"location")
                switch view_plane
                    case "hor"
                        % reorder_dim = "ZData";
                        depth(~this_bit3) = 6;
                        this_bit3_index = find(this_bit3);
                        this_bit3_index_sorted = this_bit3_index(I);
                        depth(this_bit3_index_sorted) = new_dim+3;
                    case "sag"
                        % reorder_dim = "XData";
                        x(~this_bit3) = 4.2;
                        this_bit3_index = find(this_bit3);
                        this_bit3_index_sorted = this_bit3_index(I);
                        x(this_bit3_index_sorted) = new_dim+1.5;
                end
            end

            % plot by group
            g_1 = scatter3(plot_ax,x(this_bit1),y(this_bit1),depth(this_bit1),bubbleSize,nan_style{:});
            if setuserdata
                set_userdata_and_datatip(g_1,this_bit1,fib_loc,values,p_value=grey_bubble_pvalue)
            end
            %%%%%%
            % g_2 = scatter3(plot_ax,x(this_bit2),y(this_bit2),depth(this_bit2),bubbleSize,grey_bubble_color,"filled",MarkerEdgeColor='k',MarkerFaceAlpha=no_signal_insig_alpha);
            %%%%%%
            g_2 = scatter3(plot_ax,x(this_bit2),y(this_bit2),depth(this_bit2),bubbleSize,grey_bubble_color,"filled",MarkerEdgeColor='none',MarkerFaceAlpha=no_signal_insig_alpha);
            %%%%%%
            if setuserdata
                set_userdata_and_datatip(g_2,this_bit2,fib_loc,values,p_value=grey_bubble_pvalue)
            end
            %%%%%%
            % g_3 = scatter3(plot_ax,x(this_bit3),y(this_bit3),depth(this_bit3),bubbleSize,colorful_colors(this_bit3,:),"filled",MarkerEdgeColor='none');
            %%%%%%
            g_3 = scatter3(plot_ax,x(this_bit3),y(this_bit3),depth(this_bit3),bubbleSize,colorful_colors(this_bit3,:),"filled",MarkerEdgeColor=colorful_edge_color);
            %%%%%%
            if setuserdata
                set_userdata_and_datatip(g_3,this_bit3,fib_loc,values,p_value=grey_bubble_pvalue)
            end
            if roi_label ~=0
                text(plot_ax,x(this_bit1|this_bit2|this_bit3),y(this_bit1|this_bit2|this_bit3),depth(this_bit1|this_bit2|this_bit3),string(find(this_bit1|this_bit2|this_bit3)),FontSize=roi_label,Color=roi_label_color);
            end
            gs = [g_1;g_2;g_3];
            for i = 4:size(group_bit,1)
                this_bit = group_bit(i,:);
                if any(isnan(outlineColor(i-3,:)))
                    g_tmp = scatter3(plot_ax,x(this_bit),y(this_bit),depth(this_bit),bubbleSize,MarkerEdgeColor = [0,0,0],MarkerFaceColor="none",LineWidth = outlineWidth(i-3));
                else
                    g_tmp = scatter3(plot_ax,x(this_bit),y(this_bit),depth(this_bit),bubbleSize,colorful_colors(this_bit,:),"filled",MarkerEdgeColor = outlineColor(i-3,:),LineWidth = outlineWidth);
                end
                if setuserdata
                    set_userdata_and_datatip(g_tmp,this_bit,fib_loc,values,p_value=grey_bubble_pvalue)
                end
                gs = cat(1,gs,g_tmp);
            end
        
            
            % change other property
            is_ct = ~any(AP>5);
            if isempty(DVlim)
                DVlim = [2,6];
            end
            if isempty(APlim)
                if is_ct
                    APlim = [-2.18,1.6];
                else
                    APlim = [1,20];
                end
            end
            if isempty(MLlim)
                if is_ct
                    MLlim = [0.5,4.2];
                else
                    MLlim = [1,14];
                end
            end
            set(plot_ax,'XLim',MLlim);
            set(plot_ax,'YLim',APlim);
            set(plot_ax,'ZLim',DVlim);

            set(plot_ax, 'Zdir', 'reverse');
            view(plot_ax,viewAngle);
        
            xlabel(plot_ax,'ML coord (mm)')
            ylabel(plot_ax,'AP coord (mm)')
            zlabel(plot_ax,'DV coord (mm)')
        
            n_ticks = 9;
            colormap(plot_ax,cmapcolors)
            cb = colorbar(plot_ax);
            cb.Ticks = linspace(0,1,n_ticks);
            cb.TickLabels = round(linspace(cmapping(1),cmapping(end),n_ticks),4);

            function set_userdata_and_datatip(this_plot,this_bit,ct,values,varargin)
                % build userdata for plots
                % recording ct table index, mouse name and roi index
                p_value = [];

                ip_1 = inputParser;
                ip_1.addParameter('p_value',0)
                ip_1.parse(varargin{:})
                for j_1 = fields(ip_1.Results)'
                    eval([j_1{1}, '= ip_1.Results.',j_1{1},';']);
                end

                userdata = ct(this_bit,["ROI","mouse_name","ROI_original"]);
                userdata.mouse_name = string(userdata.mouse_name);
                userdata = table2array(userdata);
                tmp = values(this_bit);
                tmp = reshape(tmp,[],1);
                userdata = cat(2,userdata,tmp);
                if ~isempty(p_value)
                    userdata = cat(2,userdata,p_value(this_bit)');
                end
                this_plot.UserData = userdata;
            
            %     this_plot.DataTipTemplate.DataTipRows = this_plot.DataTipTemplate.DataTipRows(1:3);
            %     this_plot.DataTipTemplate.DataTipRows(1).Label = "mouse name";
            %     this_plot.DataTipTemplate.DataTipRows(1).Value = userdata(:,2);
            %     this_plot.DataTipTemplate.DataTipRows(2).Label = "ROI"; 
            %     this_plot.DataTipTemplate.DataTipRows(2).Value = userdata(:,3);
            %     this_plot.DataTipTemplate.DataTipRows(3).Label = "plot id"; 
            %     this_plot.DataTipTemplate.DataTipRows(3).Value = userdata(:,1);
            end
        end
        
        function scatter_3d_with_datatip(plot_ax,values,fib_loc,varargin)
            % wrapper for scatter_3d from circle_UI_function, but add datatips 
            % (at the price of removing outline category just for complexity)
            [gs,cmapBounds_out] = common_functions.scatter_3d(plot_ax,values,fib_loc, varargin{:});
            for g_i = 1:length(gs)
                g = gs(g_i);
                datatip_cell = num2cell(g.UserData);
                dataTipTextRows = [dataTipTextRow("mouse",datatip_cell(:,2));dataTipTextRow("ROI",datatip_cell(:,3));dataTipTextRow("value",datatip_cell(:,4))];
                if size(g.UserData,2) == 5
                    dataTipTextRows = cat(1,dataTipTextRows,dataTipTextRow("p value",datatip_cell(:,5)));
                end
                gs(g_i).DataTipTemplate.DataTipRows = dataTipTextRows;
            end
        end

        function [ Fc, scale, center ] = FtoFc_exp(F,varargin)
            % An alternative to the FtoFc function to calculate DFF.
            % FtoFc.m calculates the baseline (B) using a sliding 8th percentile 
            % window, and then normalizes F to that B in the following way:
            % (F./B) - median(F./B)
            %
            % This function calculates the baseline B by fitting 2-term exponential
            % function to the F instead.
            %
            % Options: Instead of fitting the exponential to the F (raw fluorescence), 
            % you can insetad fit the exponential to a sliding baseline. To do that,
            % supply BOTH of the following optional inputs:
            %   'baseline_window'   - the number of frames for the sliding window
            %   'baseline_perc'     - the percentile you want 
            %                         (e.g.,8 for 8th percentile, 20 for 20th, etc) 
            % Mai-Anh Vu
            % 5/16/2023
            % modified by Zack

            %%%  parse optional inputs %%%
            ip = inputParser;
            ip.addParameter('baseline_window',[]);
            ip.addParameter('baseline_perc',[]);
            ip.addParameter('centering',1);
            ip.addParameter('fr',18);
            
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            
            % whether or not the Y variable is the F or a sliding baseline
            if ~isempty(baseline_window) && ~isempty(baseline_perc)    
                y_all = zeros(size(F));
                for i = 1:size(y_all,1)
                    y_all(i,:) = prctile(F(max(i-baseline_window,1):min(i+baseline_window,size(F,1)),:),baseline_perc);
                end
            else
                y_all = F;
            end
            
            % now fit the exponential to calculate the baseline
            scale = nan(size(y_all));
            x = 1:size(F,1);
            x = x(:);
            for i = 1:size(y_all,2)
                y = y_all(:,i);
                max_y = max(y);
                double_exp = @(intersect,a,b,c,d,x) intersect+a*exp(-b*x)+c*exp(-d*x);
                initiation = [max_y/2,max_y/4,1/(3600*fr),max_y/4,1/(360*fr)];
                upperbound = [max_y, max_y, 1/(100*fr),  max_y, 1/(100*fr)];
                lowerbound = [0,     0,     -1/(100*fr), 0,     -1/(100*fr)];
                mdl = fit(x,y,double_exp,StartPoint=initiation,Lower=lowerbound,Upper=upperbound);
                scale(:,i) = double_exp(mdl.intersect,mdl.a,mdl.b,mdl.c,mdl.d,x);
            end
            
            % DFF normalization
            Fc = F./scale;
            if centering
                center = median(Fc);
                Fc = Fc - center;
            end
        end
        
        function output = pav2cue_getTrialTimes(behav,varargin)
            % function output = pav2cue_getTrialTimes(behav)
            %
            % same output as cueRot_getTrialTimes
            % This function takes as input the path to a behavior file, or the struct
            % itself, and returns a struct, with fields:
            %   -trial          the trial number
            %   -cue            1 or 2 (see behav.experimentSetup.exp.cue_freq to know
            %                   which is which)
            %   -cue_ID         the frequency in kHz (-1 if LED)
            %   -level_analog   brightness or relative volume (actual level)
            %   -level_cat      brightness or relative volume (categorical, 
            %                   coded as 1, 2, 3, least to most salient)
            %   -cueStart       index for cue starts 
            %   -cueEnd         index for cue ends
            %   -rew            whether the trial was rewarded (1) or not (0)
            %   -rewOn          index for reward delivery
            %   -rewLick        index of first lick after reward delivery
            %   -trialTimes     trial time in s (how long from start of cue to end of cue)
            %   -unpredRewOn    onset index of unpredicted reward
            %   -unpredRewLIck  index for first lick after unpredicted reward
            %   -unpredRewSize  size of unpredicted rewards
            %   -fr             framerate (NOTE: ASSUMES INTEGER RATE)
            % 
            %
            % Mai-Anh, updated 12/10/2021
            % updated 2/1/2022 to handle the case where the imaging camera starts after
            %                      the task has already begun
            % updated 2/22/22 for more general use
            ip = inputParser;
            ip.addParameter('cue2Tone_shift_forward_flag',0)
            ip.parse(varargin{:});
            for j = fields(ip.Results)'
                eval([j{1}, '= ip.Results.', j{1}, ';']);
            end
            
            % load if it's a path
            if ~isstruct(behav)
                if ~endsWith(behav,'.mat')
                    behav = [behav '.mat'];
                end
                behav = load(behav);
            end
            trialInfo = behav.experimentSetup.exp.trials;
            rewInfo = getRewInfo(behav);
            
            % preallocate
            output.trial = trialInfo.trialNum;
            output.cue = trialInfo.cueRL;
            output.cue_ID = transpose(behav.experimentSetup.exp.cue_freq(trialInfo.cueRL));
            output.level_analog = transpose(behav.experimentSetup.exp.cue_relativeVol(trialInfo.cueRL));
            output.level_cat = nan(size(output.level_analog));
            output.cueStart = nan(size(output.trial));
            output.rew = trialInfo.rew;
            output.rewOn = nan(size(output.trial));
            output.rewLick = nan(size(output.trial));
            output.omitLick = nan(size(output.trial));
            output.cueEnd = nan(size(output.trial));
            
            % if it's a opto stim experiment, include stimulus into cue_opto_stim
            if isfield(behav,"stimDriver")
                output.cue_opto_stim = [];
                output.cue_opto_stim = behav.experimentSetup.exp.trials.stimulation;
            end
            
            cueFields = {'',''};
            for i = 1:2
                if behav.experimentSetup.exp.cue_freq(i)==-1
                    if isfield(behav.experimentSetup.exp,'stim_driver_freq') &&...
                            behav.experimentSetup.exp.stim_driver_freq >0
                        cueFields{i} = 'stimulus_led';
                    else
                        cueFields{i} = 'stimulus_led_analog';
                    end
                else
                    cueFields{i} = ['stimulus_sound' num2str(i)];
                end
            end
            
            % cue onsets and offsets & IDs so that we can match in case the first trial
            % isn't recorded
            stim_on = [];
            stim_off = [];
            stim_on_id = [];
            stim_off_id = [];
            stim_thresh = .1;
            for i = 1:2
                stim_field = behav.(cueFields{i});
                stim_on = [stim_on; find(diff(stim_field>stim_thresh)==1)+1];
                stim_on_id = [stim_on_id; ones(numel(find(diff(stim_field>stim_thresh)==1)),1)*i];
                stim_off = [stim_off; find(diff(stim_field>stim_thresh)==-1)+1];
                stim_off_id = [stim_off_id; ones(numel(find(diff(stim_field>stim_thresh)==-1)),1)*i];
            end
            % in case we have mismatch
            [stim_on,on_sort_idx]= sort(stim_on);
            stim_on_id = stim_on_id(on_sort_idx);
            [stim_off,off_sort_idx] = sort(stim_off);
            stim_off_id = stim_off_id(off_sort_idx);
            if stim_off(1)<stim_on(1)
                stim_off = stim_off(2:end);
                stim_off_id = stim_off_id(2:end);
            end
            if stim_on(end)>stim_off(end)
                stim_on = stim_on(1:end-1);
                stim_on_id = stim_on_id(1:end-1);
            end
            startTrial = strfind(trialInfo.cueRL', stim_on_id');
            output.cueStart(startTrial:(startTrial+numel(stim_on)-1),1) = stim_on;
            output.cueEnd(startTrial:(startTrial+numel(stim_off)-1),1) = stim_off;
            output.trialTimes(startTrial:(startTrial+numel(stim_off)-1),1) = behav.timestamp(stim_off)-behav.timestamp(stim_on);
            
            % let's add the reward onsets
            for i = 1:numel(output.trial)
                if sum(rewInfo.rewOn>output.cueStart(i) & rewInfo.rewOn<output.cueEnd(i)) == 1
                    output.rewOn(i) = rewInfo.rewOn((rewInfo.rewOn>output.cueStart(i) & rewInfo.rewOn<output.cueEnd(i)));
                    output.rewLick(i) = rewInfo.lickOnset((rewInfo.rewOn>output.cueStart(i) & rewInfo.rewOn<output.cueEnd(i)));
                end
            end
            
            % let's add lick time for omitted reward
            rew_deliv_delay = round(mean(output.rewOn-output.cueStart,"omitnan"));
            omitted_trial = find(isnan(output.rewOn));
            omitted_trial = omitted_trial(all(~isnan([output.cueStart(omitted_trial),output.cueEnd(omitted_trial)]),2));
            for i = omitted_trial'
                if all(~isnan([output.cueStart(i),output.cueEnd(i)]))
                    omit_lick = find(behav.lick(output.cueStart(i)+rew_deliv_delay:output.cueEnd(i)),1,"first");
                    if ~isempty(omit_lick)
                        output.omitLick(i) = omit_lick + output.cueStart(i)+rew_deliv_delay-1;
                    end
                end
            end
            
            % level category
            for i = 1:2
                thisCue = output.cue==i;
                levels = sort(unique(output.level_analog(thisCue)));
                for j = 1:numel(levels)
                    thisCueLevel = output.cue==i & output.level_analog == levels(j);
                    output.level_cat(thisCueLevel) = j;
                end
            end    
            
            % truncate uncompleted trials
            output_fields = fieldnames(output);
            keep_trials = ~isnan(output.cueEnd);
            for f = 1:numel(output_fields)
                output.(output_fields{f}) = output.(output_fields{f})(keep_trials);
            end
            
            % temporary modification
            % rewInfo.rewOn = rewInfo.rewOn(3:end);
            % rewInfo.lickOnset = rewInfo.lickOnset(3:end);
            
            % unpredicted rewards
            % include a +/-1 allowance to figure out which rewards were unpredicted
            [UR,idx] = setdiff(rewInfo.rewOn,output.rewOn);
            % I DONT KNOW WHY REWSIZE DIFFER FROM REWON, PAD IT
            if size(rewInfo.rewSize,1) < size(rewInfo.rewOn,1)
                sizemu = mean(rewInfo.rewSize);
                padsize = size(rewInfo.rewOn,1) - size(rewInfo.rewSize,1);
                pad = ones(padsize,1)*sizemu;
                rewInfo.rewSize = cat(1,rewInfo.rewSize,pad);
            end
            % try
            URsize = rewInfo.rewSize(idx);
            % catch err
            %     disp("err")
            %     rewInfo.rewSize = ones(48,1);
            %     URsize = rewInfo.rewSize(idx);
            % end
            output.unpredRewOn = UR;
            output.unpredRewLick = rewInfo.lickOnset(idx);
            output.unpredRewSize = URsize;
            
            % add frame rate (assume integer)
            output.fr = round(1/nanmean(diff(behav.timestamp)));
            
            % cue2 Tone are systemetically delayed by 365ms
            % (https://docs.google.com/presentation/d/1RkuxNwPmUD_h5ipeebS0nch6r66QI9VN5uPt_q_9ajk),
            % consider shift cueStart forward
            if cue2Tone_shift_forward_flag
                switch output.fr
                    case 18
            %             cue2shift = floor(18*0.55);
                        cue2shift = 9;
                    case 30
                        cue2shift = floor(30*0.365);
                    case 19
                        cue2shift = floor(19*0.37);
                    case 11
                        cue2shift = floor(11*0.4);
            %             cue2shift = 5; % actually I think 5 is better
                    otherwise
                        error("TBI")
                end
                output.cueStart(output.cue==2) = output.cueStart(output.cue==2) + cue2shift;
            end
        end

        function [data_by_mouse_phase,data_single_trial_by_mouse_phase] = get_by_mouse_phase(tas,mouse_list,phase_list)
            data_by_mouse_phase = struct;
            data_single_trial_by_mouse_phase = struct;
            for single_phase_name = phase_list
                this_tas_single = [];
                this_tas_trial_single = {};
                for single_name = mouse_list
                    if isfield(tas.(single_phase_name).(single_name),"mu_location_value_significance")
                        this_tas_mask_0 = tas.(single_phase_name).(single_name).mu_location_value_significance;
                        this_tas_mask_0(this_tas_mask_0==0) = nan;

                        this_tas_mask = repmat(this_tas_mask_0,[1,1,2]);
                        this_tas = tas.(single_phase_name).(single_name).mu_location_value .* this_tas_mask;
                        this_tas_single_trial = tas.(single_phase_name).(single_name).single_values;

                        this_tas_mask = repmat(this_tas_mask_0,[1,1,size(this_tas_single_trial,3)]);
                        this_tas_single_trial = this_tas_single_trial .* this_tas_mask;
                    else
                        this_tas = tas.(single_phase_name).(single_name).mu_location_value;
                        this_tas_single_trial = tas.(single_phase_name).(single_name).single_values;
                    end
                    this_tas_single = cat(2,this_tas_single,this_tas);
                    this_tas_trial_single = cat(2,this_tas_trial_single,num2cell(this_tas_single_trial,3));
                end
                data_by_mouse_phase.(single_phase_name).this_tas_single = this_tas_single;
                data_single_trial_by_mouse_phase.(single_phase_name).this_tas_single = this_tas_trial_single;
            end
        end

        function [p_mu,p_single] = plot_data_single(ax,plot_x,plot_act,varargin)
            % plot a single ROI of ta onto target ax
            ip = inputParser;
            ip.addParameter('plot_color',lines(1))
            ip.addParameter('sem_color',[])
            ip.addParameter('single_trial_color',[]) % sampling rate
            ip.addParameter('normalize_musem',false) % normalize mu sem via divide them by max([abs(mu+sem),abs(mu-sem)]); don't change single trials
            ip.addParameter('LineStyle',"-")
            ip.addParameter('Marker',"none")
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            
            p_single = []; p_mu = [];
            if ndims(plot_act) == 3
                plot_act = permute(plot_act,[1,3,2]);
            end
            n_trial = size(plot_act,2);
            sem = std(plot_act,[],2,'omitnan') / sqrt(n_trial);
            mu = mean(plot_act,2,'omitnan');

            if normalize_musem
                norm_factor = max([abs(mu+sem),abs(mu-sem)],[],"all","omitmissing");
                mu = mu/norm_factor; sem=sem/norm_factor;
            end
            
            if ~isempty(single_trial_color)
                p_single = plot(ax,plot_x,plot_act,Color=single_trial_color);
            end
            if ~isempty(sem_color)
                [sem_x,sem_y] = common_functions.build_xysem(plot_x,mu,sem);
                patch(ax,sem_x,sem_y,sem_color,FaceColor = sem_color,EdgeColor='none',Facealpha=0.7)
            end
            if ~isempty(plot_color)
                p_mu = plot(ax,plot_x,mu,Color = plot_color,LineWidth=1.5,LineStyle=LineStyle,Marker=Marker);
            end

            
        end

        function [x_patched,y_patched] = build_xysem(x,y,sem)
            % build fill shape for x,y,sem
            x = reshape(x,1,[]);
            y = reshape(y,1,[]);
            sem = reshape(sem,1,[]);
            x_patched = [x,fliplr(x)];
            y_patched = [y+sem,fliplr(y-sem)];
            x_patched = x_patched(1:length(y_patched));
        end

        function p_mu = plot_ta_single_just_mu_sem(ax,plot_x,ta,r,varargin)
            % plot a single ROI of ta onto target ax
            ip = inputParser;
            ip.addParameter('plot_color',lines(1))
            ip.addParameter('sem_color',[])
            ip.addParameter('single_trial_color',[]) % sampling rate
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end            

            if ~isfield(ta,"sem")
                sem = ta.std(:,r) / sqrt(n_trial);
            else
                sem = ta.sem(:,r);
            end
            mu = ta.mu(:,r);
            
            
            if ~isempty(single_trial_color)
                plot(ax,plot_x,act,Color=single_trial_color)
            end
            if ~isempty(sem_color)
                [sem_x,sem_y] = common_functions.build_xysem(plot_x,mu,sem);
                patch(ax,sem_x,sem_y,sem_color,FaceColor = sem_color,EdgeColor='none',Facealpha=0.7)
            end
            p_mu = plot(ax,plot_x,mu,Color = plot_color,LineWidth=1.5);
        end

        function ta_out = interp_ta(ta_in,sr1,sr2,varargin)
            % downsample or upsample given triggerd average (ta) with sr1 to data_out with sr2 using interp1
            % ta have dimensions: time * ROI * trial
            ip = inputParser;
            ip.addParameter("interp_method","spline") % spline,linear
            ip.parse(varargin{:});
            for j = fields(ip.Results)'
                eval([j{1}, '= ip.Results.', j{1}, ';']);
            end

            szs = size(ta_in);
            if length(szs) < 3
                szs = cat(2,szs,1);
            end
            x1 = 1:szs(1);
            duration_s = szs(1)/sr1;
            assert(mod(duration_s,1)==0,"Non integer duration detected.")
            x2 = linspace(1,szs(1),sr2*duration_s);
            
            ta_out = nan(length(x2),szs(2),szs(3));

            for r = 1:szs(2)
                temp1 = permute(ta_in(:,r,:),[1,3,2]);
                temp2 = interp1(x1,temp1,x2,interp_method);
                ta_out(:,r,:) = temp2;
            end

        end
        
        function sr = get_sr(behav)
            % get sampling rate from behav
            sr = round(1/mean(diff(behav.timestamp),'omitnan'));
        end

        function loco = get_overall_vel_acc(behav,varargin)
            % function loco = get_overall_vel_acc(behav,varargin)
            %
            % this function takes as input the behav file (path or struct)
            % it outputs overall velocity and overall acceleration
            %
            % note that overall velocity here is the sum of the absolute value of the
            % linear and angular velocity, where the angular velocity is
            % angular_velocity_ms (in m/s, NOT RADIANS).
            %
            % the resulting velocity is smoothed and lowpass filtered.
            % the acceleration is the diff of this velocity value, also smoothed and
            % lowpass filtered.
            %
            % the smoothing window, lowpass filter frequency, sampling freq, are all
            % optional inputs
            %
            % Mai-Anh Vu
            % 11/28/22
            % modifed by Liangzhu Zhang
            
            % parse optional inputs
            smooth_window = 0.3;
            lowpass_freq = 1.5;
            filter_flag = 1;
            vel_type = "total_vel_1"; % choose from various_total_vel, lin_vel or ang_vel
%             sgtext = '';

            ip = inputParser;
            ip.addParameter('smooth_window',0.3); % in seconds
            ip.addParameter('lowpass_freq',1.5) % lowpass frequency
            ip.addParameter('filter_flag',1) % filter algorithm
            ip.addParameter('vel_type',"total_vel_1") % filter algorithm
%             ip.addParameter('sgtext','') % sgtitle text
            ip.parse(varargin{:});
            for j=fields(ip.Results)'
                eval([j{1} '=ip.Results.' j{1} ';']);
            end
            
            
            % load behav if necessary
            if ischar(behav)
                if ~endsWith(behav,'.mat')
                    behav = [behav '.mat'];
                end
                behav = load(behav);
            end
            sampling_freq = behav_functions.get_sr(behav);
            % linear and angular velocity
            loco_tmp = ball2xy(behav);
            switch vel_type
                case "total_vel_1"
                    loco_tmp.velocity_mag_ms = loco_tmp.linear_velocity + abs(loco_tmp.angular_velocity_ms);
                case "lin_vel"
                    loco_tmp.velocity_mag_ms = loco_tmp.linear_velocity;
                case "ang_vel"
                    loco_tmp.velocity_mag_ms = loco_tmp.angular_velocity_ms;
            end
            
            % output vel and acc
            switch filter_flag
                case 1
                    vel_filtered = lowpass(loco_tmp.velocity_mag_ms,lowpass_freq,sampling_freq,Steepness=0.9);
                    vel_filtered_smoothed = smooth(vel_filtered,round(smooth_window*sampling_freq),"rlowess");
                    acc_filtered = lowpass(diff(vel_filtered_smoothed),lowpass_freq,sampling_freq,Steepness=0.9);
                    acc_filtered_smoothed = smooth(acc_filtered,round(smooth_window*sampling_freq),"rlowess");
                case 2
                    vel_filtered = cwt_filter(loco_tmp.velocity_mag_ms,lowpass_freq,sampling_freq);
                    vel_filtered_smoothed = smooth(vel_filtered,round(smooth_window*sampling_freq),"rlowess");
                    acc_filtered = cwt_filter(diff(vel_filtered_smoothed),lowpass_freq,sampling_freq);
                    acc_filtered_smoothed = smooth(acc_filtered,round(smooth_window*sampling_freq),"rlowess");
                case 3
                    vel_filtered = lowpass(loco_tmp.velocity_mag_ms,lowpass_freq,sampling_freq,Steepness=0.9);
                    vel_filtered_smoothed = smooth(vel_filtered,round(smooth_window*sampling_freq),"lowess");
                    acc_filtered = lowpass(diff(vel_filtered_smoothed),lowpass_freq,sampling_freq,Steepness=0.9);
                    acc_filtered_smoothed = smooth(acc_filtered,round(smooth_window*sampling_freq),"lowess");
            end
            loco.vel_original = loco_tmp.linear_velocity;
            loco.vel = vel_filtered_smoothed;
            loco.acc = [acc_filtered_smoothed(1); acc_filtered_smoothed] * sampling_freq;

            function xrec = cwt_filter(signal,lowpass_freq,sampling_freq)
                [wt,f] = cwt(signal,sampling_freq);
                xrec = icwt(wt,f,[min(f),lowpass_freq],SignalMean=mean(signal));
            end
        end
        
        function p = nan_ranksum_cell(x,y)
            % same as nan_ranksum, but x,y can be cell
            p = cellfun(@(xs,ys) common_functions.nan_ranksum(xs,ys),x,y);
        end

        function p = nan_ranksum(x,y)
            x = squeeze(x);
            y = squeeze(y);
            if all(isnan(x)) && all(isnan(y))
                p = nan;
            elseif all(isnan(x)) || all(isnan(y))
                p = -1;
            else
                p = ranksum(x,y);
                % [~,p] = ttest2(x,y);
            end
        end

        function out = nan_minus(in1,in2)
            assert(all(size(in1)==size(in2)),"Size of two inputs must match.")
            bothnan = isnan(in1) & isnan(in2);
            in1(logical(isnan(in1) - bothnan)) = 0;
            in2(logical(isnan(in2) - bothnan)) = 0;
            out = in1 - in2;
        end
        
        function [output,output_validate] = find_acc_start_end(vel,window,sr,varargin)
            % "time_start","time_acc_peak","time_end","acc_peak_value","vel_start","vel_peak","vel_end","vel_mean"
            delta_vel_threshold = 0.08;
            acc_smooth_info = {};
        
            ip = inputParser;
            ip.addParameter('delta_vel_threshold',0.08)
            ip.addParameter('acc_smooth_info',{})
            ip.parse(varargin{:})
            for j = fields(ip.Results)'
                eval([j{1},'=ip.Results.',j{1},';'])
            end
            
            n_trials = size(vel,2);
            output = nan(n_trials,8);
            output_validate = nan(size(vel));
            % find end first then go back and find peak
            [end_value,end_frame] = min(vel(window,:),[],1); end_frame = end_frame + window(1) - 1;
            for i = 1:n_trials
                if end_frame(i) <= window(1) || any(isnan(vel(window,i)))
                    warning("End of decel is at first frame of window, skipped.")
                    continue
                end
                this_end_frame = end_frame(i);
                this_vel = vel(:,i);
                [start_value,start_frame] = max(this_vel(window(1):this_end_frame-1)); start_frame = start_frame + window(1) - 1;
                if (start_value - end_value) < delta_vel_threshold
                    continue
                end
        
                % get acc
                acc = diff(this_vel); acc = [acc(1);acc];
                if ~isempty(acc_smooth_info)
                    acc = smoothdata(acc,1,acc_smooth_info{1},acc_smooth_info{2});
                end
                output_validate(:,i) = acc;
                target_acc = acc(start_frame:this_end_frame);
                [~, negative_parts] = separate_positive_negative(target_acc); negative_parts = negative_parts + start_frame - 1;
                
                n_decel = size(negative_parts,1);
                delta_vel = nan(1,n_decel);
                for j = 1:n_decel
                    delta_vel(j) = sum(acc(negative_parts(j,1):negative_parts(j,2)));
                end
                [delta_vel_value, delta_vel_I] = min(delta_vel);
                try
                if  isempty(delta_vel_value) || (abs(delta_vel_value) < delta_vel_threshold)
                    % warning("Probably no continouse acceleration.")
                    acc_value = nan; acc_frame = nan; vel_pk = nan;
                else
                    this_vel_id = negative_parts(delta_vel_I,:);
                    [acc_value,acc_frame] = min(acc(this_vel_id(1):this_vel_id(2))); acc_frame = acc_frame + this_vel_id(1) - 1;
                    vel_pk = this_vel(acc_frame);
                end    
                output(i,:) = [start_frame/sr,acc_frame/sr,this_end_frame/sr,acc_value,start_value,vel_pk,end_value(i),mean(vel(start_frame:this_end_frame))];
                catch err
                    keyboard
                    err.rethrow
                end
            end

            function [positive_parts, negative_parts] = separate_positive_negative(array)
                array = reshape(array,1,[]);
                % Check if the array is empty or has only one element
                if isempty(array) || numel(array) == 1
                    return;
                end
            
                % Find the indices where the sign of the elements change
                sign_changes = diff(sign(array));
            
                % Find the indices where the positive parts start and end
                positive_starts = find(sign_changes > 0) + 1;
                positive_ends = find(sign_changes < 0);
            
                % Find the indices where the negative parts start and end
                negative_starts = find(sign_changes < 0) + 1;
                negative_ends = find(sign_changes > 0);
            
                % If the array starts with a positive part, add the first positive part
                if sign(array(1)) == 1
                    positive_starts = [1, positive_starts];
                end
            
                % If the array ends with a positive part, add the last positive part
                if sign(array(end)) == 1
                    positive_ends = [positive_ends, numel(array)];
                end
            
                % If the array starts with a negative part, add the first negative part
                if sign(array(1)) == -1
                    negative_starts = [1, negative_starts];
                end
            
                % If the array ends with a negative part, add the last negative part
                if sign(array(end)) == -1
                    negative_ends = [negative_ends, numel(array)];
                end
            
                positive_parts = [positive_starts;positive_ends]';
                negative_parts = [negative_starts;negative_ends]';
            end

        end
        
        function include_mask = get_aDMS_fiber(this_CT)
            include_mask = this_CT{:,"fiber_bottom_AP"}>0 & this_CT{:,"fiber_bottom_ML"}<2 & this_CT{:,"fiber_bottom_DV"}<4;
        end

        function out = decimal_to_str(input,num_decimal)
            % convert number arr to string arr
            out = nan(size(input));
            for xi = 1:length(input)
                out(xi) = sprintf("%."+num_decimal+"f",input(xi));
            end
        end
        
        function asters = p_to_asterisk(ps)
            % convert p value to asterisk
            asters = repmat("",size(ps));
            for pi = 1:length(ps)
                p = ps(pi);

                if p > 0.05
                    ast = "nan";
                elseif p > 0.01 && p <= 0.05
                    ast = "*";
                elseif p > 0.001 && p <= 0.01
                    ast = "**";
                elseif p <= 0.001
                    ast = "***";
                else
                    ast = "nan";
                end
                asters(pi) = ast;
            end
        end
        
        function [plot_line,plot_text] = add_significance_to_ax(ax,x_pos,y_pos,asterisks)
            % add asterisks of significance to existing ax
            % y_pos is y_height in fraction
            % x_pos is a cell contains pairs of x position
            % e.g.
            % y_pos = [1,1.05,1.15,1,1.1,1];
            % x_pos = {[1.1,1.9],[1.1,2.9],[1.1,3.9],[2.1,2.9],[2.1,3.9],[3.1,3.9]};

            tmp_ys = random_functions.get_yaxis_height(ax,y_pos);
            tmp_ys1 = random_functions.get_yaxis_height(ax,y_pos+0.01);
            
            for i = 1:length(y_pos)
                plot_line = plot(ax,x_pos{i},[tmp_ys(i),tmp_ys(i)],'k-');
                plot_text = text(ax,mean(x_pos{i}),tmp_ys1(i),asterisks(i),FontSize=12);
            end

            function ys = get_yaxis_height(ax,fracs)
                % get corresponding y value of ax based on fracs
                this_ylim = ylim(ax);
                ys = this_ylim(1) + abs(diff(this_ylim))*fracs;
            end
        end
        
        function output = fiber_inpolygon(APMLDV,pos_hor,pos_sag,neg_hor,neg_sag)
            % find fibers inside the intersect of 2D polygon contours
            output = repmat({false(size(APMLDV,1),1)},1,4); % {pos_hor,pos_sag,neg_hor,neg_sag}
            input = {pos_hor,pos_sag,neg_hor,neg_sag};
        
            for pni=1:2
                for this_ploygon = input{2*pni-1}
                    for contour_i = 1:length(this_ploygon)
                        [tmp1,tmp2] = inpolygon(APMLDV(:,1),APMLDV(:,2),this_ploygon{contour_i}(:,2),this_ploygon{contour_i}(:,1));
                        output{2*pni-1} = output{2*pni-1} | tmp1 | tmp2;
                    end
                end
                for this_ploygon = input{2*pni}
                    for contour_i = 1:length(this_ploygon)
                        [tmp1,tmp2] = inpolygon(APMLDV(:,1),APMLDV(:,3),this_ploygon{contour_i}(:,1),this_ploygon{contour_i}(:,2));
                        output{2*pni} = output{2*pni} | tmp1 | tmp2;
                    end
                end
            end
        end

        % colormap function from file_exchange @ matlab
        function c = redblue(m)
            %REDBLUE    Shades of red and blue color map
            %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
            %   The colors begin with bright blue, range through shades of
            %   blue to white, and then through shades of red to bright red.
            %   REDBLUE, by itself, is the same length as the current figure's
            %   colormap. If no figure exists, MATLAB creates one.
            %
            %   For example, to reset the colormap of the current figure:
            %
            %             colormap(redblue)
            %
            %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
            %   COLORMAP, RGBPLOT.
            
            %   Adam Auton, 9th October 2009
            
            if nargin < 1, m = size(get(gcf,'colormap'),1); end
            
            if (mod(m,2) == 0)
                % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
                m1 = m*0.5;
                r = (0:m1-1)'/max(m1-1,1);
                g = r;
                r = [r; ones(m1,1)];
                g = [g; flipud(g)];
                b = flipud(r);
            else
                % From [0 0 1] to [1 1 1] to [1 0 0];
                m1 = floor(m*0.5);
                r = (0:m1-1)'/max(m1,1);
                g = r;
                r = [r; ones(m1+1,1)];
                g = [g; 1; flipud(g)];
                b = flipud(r);
            end
            
            c = [r g b]; 
        end

        % meta data
        function out = get_training_info()
            % training session infos
            % VariableNames: mouse name, last day of acquisition/devaluation 1/devaluation 2, devaluation 1 cue, devaluation 2 cue
            out = ...
                {"G12"	12	16	21	"1 Tone"	"2 LED";
                "G15"	7	12	18	"1 Tone"	"2 LED";
                "G17"	8	15	21	"1 Tone"	"2 LED";
                "G19"	9	13	18	"1 LED"	"2 Tone";
                "G22"	8	13	18	"1 LED"	"2 Tone";
                "G21"	10	15	20	"1 Tone"	"2 LED";
                "G23"	10	15	20	"1 LED"	"2 Tone";
                "G24"	10	15	20	"1 Tone"	"2 LED";
                "G25"	10	14	14	"1 LED"	"2 Tone";   % note that G25-G27 doesn't have Tone devaluation
                "G26"	15	21	21	"1 LED"	"2 Tone";   % note that G25-G27 doesn't have Tone devaluation
                "G27"	9	13	13	"1 LED"	"2 Tone";}; % note that G25-G27 doesn't have Tone devaluation
        end
        
        function output = get_training_info_DA(m_name)
            % get session info and fiber loc by mouse_name
            output = dictionary;
            switch m_name
                case "DL18"
                    output("pav") = {1:19};
                    output("LED_omi") = {25:32};
                    output("Tone_omi") = {20:24};
                case "DL20"
                    output("pav") = {1:14};
                    output("LED_omi") = {20:24};
                    output("Tone_omi") = {15:19};
                case "DL21"
                    output("pav") = {1:14};
                    output("LED_omi") = {15:19};
                    output("Tone_omi") = {20:26};
                case "DL23"
                    output("pav") = {1:14};
                    output("LED_omi") = {20:24};
                    output("Tone_omi") = {15:19};
            end
        end

        function output = get_training_info_Glu(m_name)
            % divide learning phases into pre, post and omission learnings
            output = dictionary;
            switch m_name
                case "glu924"
                    output("pre") = {1:4};
                    output("post") = {9:12};
                    output("LED_omi") = {[14:15,17]};
                    output("LED_omi_complete") = {[13:15,17]};
                case "glu926"
                    output("pre") = {1:4};
                    output("post") = {9:12};
                    % output("LED_omi") = {nan};
                    output("LED_omi") = {14:17};
                    output("LED_omi_complete") = {13:17};
            end
        end

        function include_days = get_include_days(varargin)
            % these are identified pre/post acquisition and learned devaluation days
            % based on ranksum test of licking index
            info_flag = 0;

            ip = inputParser;
            ip.addParameter("info_flag",0);
            ip.parse(varargin{:})
            for j = fields(ip.Results)'
                eval([j{1},'=ip.Results.',j{1},';'])
            end
            
            training_info = common_functions.get_training_info();

            include_days = struct;
            switch info_flag
                case 0 % original classical version
                    include_days.G12early = [1,2];include_days.G12late = [6,7,8,9,11,12];
                    include_days.G15early = [1,2];include_days.G15late = [3:7];
                    
                    include_days.G17early = [1];include_days.G17late = [7,8,9,10,11];
                    include_days.G19early = [1,2];include_days.G19late = [6,8,9];
                    include_days.G21early = [1,2];include_days.G21late = [8,10];
                    include_days.G22early = [1];include_days.G22late = [5,6,7,8];
                    include_days.G23early = [1,2,3,4];include_days.G23late = [9,10];
                    include_days.G24early = [1];include_days.G24late = [8,9,10];
                    
                    include_days.G25early = [1:4];include_days.G25late = [7:10];
                    include_days.G26early = [1:7];include_days.G26late = [15];
                    include_days.G27early = [1:6];include_days.G27late = [8];
                case 1 % omit early, put all pavday into late
                    include_days.G12late = [1:12];
                    include_days.G13late = [1:12];
                    include_days.G15late = [1:7];
                    include_days.G17late = [1:8];
                    include_days.G19late = [1:9];
                    include_days.G21late = [1:10];
                    include_days.G22late = [1:8];
                    include_days.G23late = [1:10];
                    include_days.G24late = [1:10];
                    include_days.G25late = [1:10];
                    include_days.G26late = [1:15];
                    include_days.G27late = [1:9];
            end
            
            include_days.G12LEDomi = [2:5]; include_days.G12LEDomi = include_days.G12LEDomi + tmp_get_day(training_info,"G12",2);
            include_days.G15LEDomi = [3:6]; include_days.G15LEDomi = include_days.G15LEDomi + tmp_get_day(training_info,"G15",2);
            
            include_days.G17LEDomi = [2:6]; include_days.G17LEDomi = include_days.G17LEDomi + tmp_get_day(training_info,"G17",2);
            include_days.G19LEDomi = [2:4]; include_days.G19LEDomi = include_days.G19LEDomi + tmp_get_day(training_info,"G19",1);
            include_days.G21LEDomi = [2:5]; include_days.G21LEDomi = include_days.G21LEDomi + tmp_get_day(training_info,"G21",2);
            include_days.G22LEDomi = [2:5]; include_days.G22LEDomi = include_days.G22LEDomi + tmp_get_day(training_info,"G22",1);
            include_days.G23LEDomi = [2:5]; include_days.G23LEDomi = include_days.G23LEDomi + tmp_get_day(training_info,"G23",1);
            include_days.G24LEDomi = [2:5]; include_days.G24LEDomi = include_days.G24LEDomi + tmp_get_day(training_info,"G24",2);
            
            include_days.G25LEDomi = [2:4]; include_days.G25LEDomi = include_days.G25LEDomi + tmp_get_day(training_info,"G25",1);
            include_days.G26LEDomi = [5:6]; include_days.G26LEDomi = include_days.G26LEDomi + tmp_get_day(training_info,"G26",1);
            include_days.G27LEDomi = [4]; include_days.G27LEDomi = include_days.G27LEDomi + tmp_get_day(training_info,"G27",1);
            
            include_days.G12Toneomi = [2:4]; include_days.G12Toneomi = include_days.G12Toneomi + tmp_get_day(training_info,"G12",1);
            include_days.G15Toneomi = [2:5]; include_days.G15Toneomi = include_days.G15Toneomi + tmp_get_day(training_info,"G15",1);
            
            include_days.G17Toneomi = [3:7]; include_days.G17Toneomi = include_days.G17Toneomi + tmp_get_day(training_info,"G17",1);
            include_days.G19Toneomi = [2:5]; include_days.G19Toneomi = include_days.G19Toneomi + tmp_get_day(training_info,"G19",2);
            include_days.G21Toneomi = [2:5]; include_days.G21Toneomi = include_days.G21Toneomi + tmp_get_day(training_info,"G21",1);
            include_days.G22Toneomi = [2:5]; include_days.G22Toneomi = include_days.G22Toneomi + tmp_get_day(training_info,"G22",2);
            include_days.G23Toneomi = [3:5]; include_days.G23Toneomi = include_days.G23Toneomi + tmp_get_day(training_info,"G23",2);
            include_days.G24Toneomi = [3:5]; include_days.G24Toneomi = include_days.G24Toneomi + tmp_get_day(training_info,"G24",1);

            function out=tmp_get_day(info,mname,i)
                out=info{[info{:,1}]==mname,i+1};
            end
        end
        
        function output = get_include_days_DA(m_name)
            % divide learning phases into pre, post and omission learnings
            output = dictionary;
            switch m_name
                case "DL18"
                    output("pre") = {2:3};
                    output("post") = {[13:17,19]};
                    output("LED_omi") = {29:32};
                    output("Tone_omi") = {21:24};
                    output("LED_omi_complete") = {25:32};
                    output("Tone_omi_complete") = {20:24};
                case "DL20"
                    output("pre") = {1:2};
                    output("post") = {6:14};
                    output("LED_omi") = {21:24};
                    output("Tone_omi") = {16:19};
                    output("LED_omi_complete") = {20:24};
                    output("Tone_omi_complete") = {15:19};
                case "DL21"
                    output("pre") = {1:4};
                    output("post") = {7:14};
                    output("LED_omi") = {16:19};
                    output("Tone_omi") = {22:26};
                    output("LED_omi_complete") = {15:19};
                    output("Tone_omi_complete") = {20:26};
                case "DL23"
                    output("pre") = {1:11};
                    output("post") = {12:14};
                    output("LED_omi") = {21:24};
                    output("Tone_omi") = {16:19};
                    output("LED_omi_complete") = {20:24};
                    output("Tone_omi_complete") = {15:19};
            end
        end

        function main_sr = get_main_samplerate()
            main_sr = dictionary(...
                ["G12","G15","G17","G19","G22","G21","G23","G24","G25","G26","G27","glu924","glu926"],...
                [30,30,repmat(18,[1,11])]);
        end
        
        function field_windows = get_comp_timewindow()
            % this component window is determined based on triggered
            % avearage and timing distribution histograms
            field_windows = struct;
            field_windows.unpred = [0.8,1.3;1.2,1.9;1.4,2.4];
            
            field_windows.rew1early = [0.8,1.4;1.2,1.8;1.5,2.2];
            field_windows.rew1late = [0.8,1.4;1.2,1.8;1.5,2.2];
            field_windows.rew1LEDomi = [0.8,1.4;1.2,1.8;1.5,2.2];
            field_windows.rew1Toneomi = [0.8,1.4;1.2,1.8;1.5,2.2];
            
            field_windows.rew2early = [0.8,1.4;1.2,1.8;1.5,2.2];
            field_windows.rew2late = [0.8,1.4;1.2,1.8;1.5,2.2];
            field_windows.rew2LEDomi = [0.8,1.5;1.2,1.8;1.5,2.2];
            field_windows.rew2Toneomi = [0.8,1.4;1.1,1.8;1.5,2.2];
            
            field_windows.cue2Toneomi = [19,26;25,33;30,38]/18;
            field_windows.cue2LEDomi = [19,26;25,33;30,38]/18;
            field_windows.cue2late = [19,26;25,33;30,38]/18;
            field_windows.cue2early = [19,26;25,33;30,38]/18;
            
            field_windows.cue1Toneomi = [19,25;24,31;27,40]/18;
            field_windows.cue1LEDomi = [19,25;24,31;27,40]/18;
            field_windows.cue1late = [19,25;24,31;27,40]/18;
            field_windows.cue1early = [19,25;24,31;27,40]/18;
        end
        
        function comp_time_window = get_ITI_licking_comp_timewindow()
            comp_time_window = [-0.27,0.23;0,0.6;0.2,0.8];
        end
        
        function out = DA_get_comp_timewindow()
            % get components window of DA
            out = [0,1;0.2,1.5;]; % this dip window up to 1.5 is very loose already
        end

        function out = DA_get_comp_timewindow_ITI()
            % get components window of DA for ITI licking
            out = [-0.2,0.5;0,0.5;];
        end

    end % Static methods end

end