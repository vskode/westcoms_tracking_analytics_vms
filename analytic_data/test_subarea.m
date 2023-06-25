% v09_02_beachedPartsOnly_subarea_3percent

clear 
close all

% load mesh 
load ../data/westcoms2_mesh

%% ** set up data ** 
addpath('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/01_ModelRuns_ClydePaper/v09_02_3percent');

folder_dir = ('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/01_ModelRuns_ClydePaper/v09_02_3percent'); 
input_folder_dir = ('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/01_ModelRuns_ClydePaper/v09_02_3percent');

% first create the polygon for the sub area 
sub_area1 = polyshape([-7.9 -4.366 -4.366 -7.9], [54.1 54.1 58.64 58.64]); % full coast no boundary 
sub_area2 = polyshape([-7.9 -4.3665 -4.3665 -7.9], [54.625 54.625 58.64 58.64]); % north and extends to solway coast 
sub_area3 = polyshape([-7.9 -4.23 -4.23 -7.9], [55.06 55.06 58.64 58.64]); % north and extends to north channel 
sub_area4 = polyshape([-7.1 -4.3665 -4.3665 -7.1], [55.06 55.06 57.77 57.77]); % skye to north channel 

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

    % select only points that equal 1 (beached)
    onlybeached = fate(fate == 1,:); 

    % use only x and y 
    onlybeached = onlybeached(:,[2 3]); 

    % save x and y individually 

    beachedx = onlybeached(:,1);
    beachedy = onlybeached(:,2); 
    
    Nbeachedx(:,i) = numel(beachedx); 
    Nbeachedy(:,i) = numel(beachedy); 
    
    % Clyde Sea 
    sub_areax1 = sub_area1.Vertices(:,1); 
    sub_areay1 = sub_area1.Vertices(:,2); 

    sub_areaxq1 = beachedx;
    sub_areayq1 = beachedy; 

    [insub_area1] = inpolygon(sub_areaxq1,sub_areayq1,sub_areax1,sub_areay1);
    
    sub_x1 = sub_areaxq1(insub_area1);
    sub_y1 = sub_areayq1(insub_area1);
    
    Nsub_x1(:,i) = numel(sub_x1);
    Nsub_y1(:,i) = numel(sub_y1);
   
    sub_areax2 = sub_area2.Vertices(:,1); 
    sub_areay2 = sub_area2.Vertices(:,2); 
    
    sub_areaxq2 = beachedx;
    sub_areayq2 = beachedy; 

    [insub_area2] = inpolygon(sub_areaxq2,sub_areayq2,sub_areax2,sub_areay2);
    
    sub_x2 = sub_areaxq2(insub_area2);
    sub_y2 = sub_areayq2(insub_area2);
    
    Nsub_x2(:,i) = numel(sub_x2);
    Nsub_y2(:,i) = numel(sub_y2);
    
    sub_areax3 = sub_area3.Vertices(:,1); 
    sub_areay3 = sub_area3.Vertices(:,2); 
    
    sub_areaxq3 = beachedx;
    sub_areayq3 = beachedy; 

    [insub_area3] = inpolygon(sub_areaxq3,sub_areayq3,sub_areax3,sub_areay3);
    
    sub_x3 = sub_areaxq3(insub_area3);
    sub_y3 = sub_areayq3(insub_area3);
    
    Nsub_x3(:,i) = numel(sub_x3);
    Nsub_y3(:,i) = numel(sub_y3);
    
    sub_areax4 = sub_area4.Vertices(:,1); 
    sub_areay4 = sub_area4.Vertices(:,2); 
    
    sub_areaxq4 = beachedx;
    sub_areayq4 = beachedy; 

    [insub_area4] = inpolygon(sub_areaxq4,sub_areayq4,sub_areax4,sub_areay4);
    
    sub_x4 = sub_areaxq4(insub_area4);
    sub_y4 = sub_areayq4(insub_area4);
    
    Nsub_x4(:,i) = numel(sub_x4);
    Nsub_y4(:,i) = numel(sub_y4);
    end 

    % double check method by calculating beached per day for full region 
    for i = 2:length(Nbeachedx)
    M_dailybeached_3percent(i,1) = Nbeachedx(1,i) - Nbeachedx(1,i-1);
    end 
    
    % now use the newly created variable to calculate how many beach per day for the new sub-area 
    for i = 2:length(Nsub_x1)
    M_dailybeached_3percent_SUB(i,1) = Nsub_x1(1,i) - Nsub_x1(1,i-1);
    end 

    figure(1) 
    plot(Nsub_x1, 'g-')
    hold on 
    plot(Nsub_x2, 'm-')
    hold on 
    plot(Nsub_x3, 'b-')
    hold on 
    plot(Nsub_x4, 'c-')
    xlabel('Time (days)');
    ylabel('No. of Beached Particles per Day');
    title('Resuspending - Classified 2020');
    lgd = legend;
    legend({'Full Coast - without boundary','North - Solway Coast', 'North to North Channel, incl outer hebs', 'Skye to North Channel, no outer hebs'},'Location','northwest','Orientation','vertical')
    lgd.FontSize = 14;
    fig_name = ('subarea_TOTbeached_3percent');
    print('-f1', '-dpng', '-loose', '-r500',...
        ['v09_02_plot_' fig_name '.png']);
    
% now save as a new .mat file containing new data 
save v09_02_3percent_xyBeachedYear Nbeachedx Nbeachedy Nsub_x Nsub_y