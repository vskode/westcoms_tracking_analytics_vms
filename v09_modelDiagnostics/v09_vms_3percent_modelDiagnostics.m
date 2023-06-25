% v09_02_beachedPartsOnly_subarea_3percent

% clear 
% close all

%% ** set up data ** 
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_data');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms/analytic_data');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms/v09_modelDiagnostics');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_trackingOutput/VMS/v09_3percent_vms');

% load mesh 
load westcoms2_mesh
folder_dir = ('E:/2020_particle_track_models/2021_WestComs_na/westcoms_trackingOutput/VMS/v09_3percent_vms/'); 
input_folder_dir = ('E:/2020_particle_track_models/2021_WestComs_na/westcoms_trackingOutput/VMS/v09_3percent_vms');

% do for jan only
% folder_dir = ('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/01_ModelRuns_ClydePaper/v09_02_3percent_jan/'); 
% input_folder_dir = ('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/01_ModelRuns_ClydePaper/v09_02_3percent_jan');

% first create the polygon for the sub area 
sub_area = polyshape([-7.8 -4.23 -4.23 -7.8], [55.06 55.06 58.6 58.6]); % north and extends to north channel 
north_sub_area = polyshape([-7.9 -4.23 -4.23 -7.9], [58.64 58.64 59.3 59.3]);
open_bound = polyshape([-9 -7.9 -7.9 -4.23 -4.23 -9], [55.06 55.06 58.64 58.64 59.3 59.3]);

%% ******* concatenate the files together by x and y positions *******

% set the directory to the input folder 
mat_files = dir(input_folder_dir);

% create the file list in a cell array 
matfile_list = {mat_files.name};
matfile_list(1:2) = []; 

% loop through each file in the directory from 1 to N
    for i = 1:(length(matfile_list))
    
    % get the file name from the struct (need curly brackets) 
    file_1 = matfile_list{i}; % read first file in the list 
    
    % load data file 
    data_1 = load(file_1);

    % count how many ashore at the end of the model run 

    onbeach = numel(data_1.ashore(data_1.ashore));

    % concatenate the positions and logical variable of ashore 
    fate = [data_1.ashore, data_1.x, data_1.y]; 
    
    % concatenate the positions and logical variable of ashore 
    total = [data_1.x, data_1.y];
    
    % select only points that equal 1 (beached)
    onlybeached = fate(fate == 1,:); 

    % use only x and y 
    onlybeached = onlybeached(:,[2 3]); 

    %%%% find total particles (active and ashore) 
    % save x and y individually for the total data set (here i need to extract particles for the 1st year only) 

    totalx_active = data_1.x(data_1.active);
    totaly_active = data_1.y(data_1.active);
    totalx_ashore = data_1.x(data_1.ashore);
    totaly_ashore = data_1.y(data_1.ashore);
    
    totalx = vertcat(totalx_active, totalx_ashore);
    totaly = vertcat(totaly_active, totaly_ashore);

    totalACTIVE = numel(totalx);
    
    %%%%%%%%% how many active particles are in the new sub-area  
    sub_areax = sub_area.Vertices(:,1); 
    sub_areay = sub_area.Vertices(:,2); 

    sub_areaxqTOT = totalx;
    sub_areayqTOT = totaly; 
    
    [insub_areaTOT] = inpolygon(sub_areaxqTOT,sub_areayqTOT,sub_areax,sub_areay);

    sub_xTOT = sub_areaxqTOT(insub_areaTOT);
    sub_yTOT = sub_areayqTOT(insub_areaTOT);
    
    NTOTsub_x3(:,i) = numel(sub_xTOT);
    NTOTsub_y3(:,i) = numel(sub_yTOT);
    
    totalsub_area = numel(sub_xTOT);
    
    %%%%%%%% total out of northern sub-area boundary 
    north_sub_areax = north_sub_area.Vertices(:,1); 
    north_sub_areay = north_sub_area.Vertices(:,2); 

    north_sub_areaxqTOT = totalx;
    north_sub_areayqTOT = totaly; 
    
    [innorth_sub_areaTOT] = inpolygon(north_sub_areaxqTOT,north_sub_areayqTOT,north_sub_areax,north_sub_areay);
    
    north_sub_xTOT = north_sub_areaxqTOT(innorth_sub_areaTOT);
    north_sub_yTOT = north_sub_areayqTOT(innorth_sub_areaTOT);
    
    north_NTOTsub_x3(:,i) = numel(north_sub_xTOT);
    north_NTOTsub_y3(:,i) = numel(north_sub_yTOT);
    
    totalNorth = numel(north_sub_xTOT);
    
    %%%%%%%% total out of all sub-area boundary (north and west) 
    open_bound_sub_areax = open_bound.Vertices(:,1); 
    open_bound_sub_areay = open_bound.Vertices(:,2); 

    open_bound_sub_areaxqTOT = totalx;
    open_bound_sub_areayqTOT = totaly; 
    
    [inopen_bound_sub_areaTOT] = inpolygon(open_bound_sub_areaxqTOT,open_bound_sub_areayqTOT,open_bound_sub_areax,open_bound_sub_areay);
    
    open_bound_sub_xTOT = open_bound_sub_areaxqTOT(inopen_bound_sub_areaTOT);
    open_bound_sub_yTOT = open_bound_sub_areayqTOT(inopen_bound_sub_areaTOT);
    
    open_bound_NTOTsub_x3(:,i) = numel(open_bound_sub_xTOT);
    open_bound_NTOTsub_y3(:,i) = numel(open_bound_sub_yTOT);
    
    totalopen_bound = numel(open_bound_sub_xTOT);
    
    %%%%%% now look at beached in sub-area
    beachedx = onlybeached(:,1);
    beachedy = onlybeached(:,2); 
    
    Nbeachedx(:,i) = numel(beachedx); 
    Nbeachedy(:,i) = numel(beachedy); 
    
    %%%%%%%% now how many are beached with-in the sub-area 
    sub_areaxq = beachedx;
    sub_areayq = beachedy; 

    [insub_area] = inpolygon(sub_areaxq,sub_areayq,sub_areax,sub_areay);
    
    sub_x = sub_areaxq(insub_area);
    sub_y = sub_areayq(insub_area);
    
    Nsubx_1min_vms(:,i) = numel(sub_x);
    Nsuby_1min_vms(:,i) = numel(sub_y);
    
    totalBEACHED_tot = numel(beachedx); 
    totalBEACHED_sub = numel(sub_x);
    
    totalSUB = totalsub_area; 
    end 
    
    for i = 2:length(open_bound_NTOTsub_x3)
    M_daily_open_bound_1percent(i,1) = open_bound_NTOTsub_x3(1,i) - open_bound_NTOTsub_x3(1,i-1);
    end 
    
    M_daily_open_bound_1percent = abs(M_daily_open_bound_1percent);
%     totalopen_bound_ACTIVE = sum(M_daily_open_bound_1percent);
    
    %%%% calculate percentages 
    percentNorth_bound = totalNorth/totalACTIVE*100;
    percentOpen_bound = totalopen_bound/totalACTIVE*100; 
    percent_beached_sub = totalBEACHED_sub/totalACTIVE*100;
    percent_beached_tot = totalBEACHED_tot/totalACTIVE*100; 
    
    % double check method by calculating beached per day for full region 
    for i = 2:length(Nbeachedx)
    M_dailybeached_3percent_vms(i,1) = Nbeachedx(1,i) - Nbeachedx(1,i-1);
    end 
    
    % now use the newly created variable to calculate how many beach per day for the new sub-area 
    for i = 2:length(Nsubx_1min_vms)
    M_dailybeached_3percent_vms(i,1) = Nsubx_1min_vms(1,i) - Nsubx_1min_vms(1,i-1);
    end 

    figure(1) 
    plot(Nsubx_1min_vms, 'b-')
    hold on 
    xlabel('Time (days)');
    ylabel('Total Beached Particles Over Time');
%     title('Resuspending - Classified 2020');
    lgd = legend;
    legend({'Resuspending - Classified 2020'},'Location','northwest','Orientation','vertical')
    lgd.FontSize = 14;
    fig_name = ('totalBeached_3percent');
    print('-f1', '-dpng', '-loose', '-r500',...
        ['v09_vms_plot_' fig_name '.png']);
    
        figure(2) 
    plot(M_dailybeached_3percent_vms, 'b-')
    hold on 
    xlabel('Time (days)');
    ylabel('Daily Beached Particles Over Time');
%     title('Resuspending - Classified 2020');
    lgd = legend;
    legend({'Resuspending - Classified 2020'},'Location','northwest','Orientation','vertical')
    lgd.FontSize = 14;
    fig_name = ('dailyBeached_3percent');
    print('-f2', '-dpng', '-loose', '-r500',...
        ['v09_vms_plot_' fig_name '.png']);
    
% now save as a new .mat file containing new data 
save v09_vms_3percent_modelDiagnostics Nsubx_1min_vms M_dailybeached_3percent_vms