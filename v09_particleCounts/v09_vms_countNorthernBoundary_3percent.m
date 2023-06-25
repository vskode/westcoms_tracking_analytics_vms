%% - ** count how many particles go north within a model run - **

% - in this scenario the boundary is the Clyde Sea 
% - so how many make it out and go north 
% set boundaries of ocean regions (same as MCS) 
% then count how many particles inside each 

close all
clear 

addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_data');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms/analytic_data');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms/v09_modelDiagnostics');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_trackingOutput/VMS/v09_3percent_vms');

% load beached particle positions 
v09_vms_3percent_modelDiagnostics
 
%load polygons 
polygons_northernBoundary

westcoast = load('E:\2020_particle_track_models\2021_WestComs_na\westcoms_data\ukireland_coastline.NaN.dat');
l = find(isnan(westcoast(:,1)));

lonlim = [-10 -4];
latlim = [54 59.2]; 

% iactive = find(active);   % a list of the active particles
% ibeached = find(ashore); % a list of beached particles 
% nactive = numel(iactive); % the number of active particles
% nashore = numel(iashore); % the number of beached particles

%% find total particles in each polygon

% Clyde Sea 
clydeTOTx = Clyde.Vertices(:,1); 
clydeTOTy = Clyde.Vertices(:,2); 

clydeTOTxq = totalx; % total_x is floating and beached particles 
clydeTOTyq = totaly; 

[inClyde_total] = inpolygon(clydeTOTxq,clydeTOTyq,clydeTOTx,clydeTOTy);

% North of Clyde Sea 
northTOTx = GoingNorth.Vertices(:,1); 
northTOTy = GoingNorth.Vertices(:,2); 

northTOTxq = totalx;
northTOTyq = totaly; 

[inNorth_total] = inpolygon(northTOTxq,northTOTyq,northTOTx,northTOTy);

%% find beached particles in each polygon
% Clyde Sea 
clydex = Clyde.Vertices(:,1); 
clydey = Clyde.Vertices(:,2); 

clydexq = sub_x; % sub_x is beached particles only 
clydeyq = sub_y; 

[inClyde_beached] = inpolygon(clydexq,clydeyq,clydex,clydey);

% North of Clyde Sea 
northx = GoingNorth.Vertices(:,1); 
northy = GoingNorth.Vertices(:,2); 

northxq = sub_x;
northyq = sub_y; 

[inNorth_beached] = inpolygon(northxq,northyq,northx,northy);

%% count how many particles in each polygon 
inclyde_total = numel(clydeTOTxq(inClyde_total)); 
innorth_total = numel(northTOTxq(inNorth_total)); 

inclyde_beached = numel(clydexq(inClyde_beached)); 
innorth_beached = numel(northxq(inNorth_beached)); 

outclyde = numel(clydeTOTxq(~inClyde_total));
totalSouth = outclyde - innorth_total;

% clyde = totalACTIVE-outclyde;

%% PLOT FIGURE 
figure(1)

clf;
hold on

% first draw the coast 
for i = 1:size(l)-1
    plot(westcoast(l(i) + 1:l(i+1)-1, 1), westcoast(l(i) + 1:l(i+1)-1, 2), 'k')
end 

%Clyde
plot(clydex,clydey)
axis equal 

hold on 
plot(clydexq(inClyde_beached),clydeyq(inClyde_beached), 'g.')
plot(clydexq(~inClyde_beached),clydeyq(~inClyde_beached), 'b.')

set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

hold on

%North of Clyde Sea 
plot(northx,northy)
axis equal 

hold on 
plot(northxq(inNorth_beached),northyq(inNorth_beached), 'r.')
%plot(northxq(~inNorth),northyq(~inNorth), 'g.')
set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

%% calculate % in each polygon 
% percentleftClyde = outclyde/totalACTIVE*100;
percentinClydeTOT = inclyde_total/totalACTIVE*100;
percentNorthTOT = innorth_total/totalACTIVE*100;
percentinClyde = inclyde_beached/totalACTIVE*100;
percentNorth = innorth_beached/totalACTIVE*100;
percentbeached = totalBEACHED_sub/totalACTIVE*100;
percentfloating = (totalSUB-totalBEACHED_sub)/totalACTIVE*100;
percentClydefloating = (inclyde_total-inclyde_beached)/totalACTIVE*100;
percentNorthfloating = (innorth_total-innorth_beached)/totalACTIVE*100;
percent_inSUB = totalSUB/totalACTIVE*100;

% percent beached in sub area of total particles in sub area 
percent_beachedSUB = totalBEACHED_sub/totalSUB*100;
percentinClydeSUB = inclyde_beached/totalSUB*100;
percentNorthSUB = innorth_beached/totalSUB*100;
percent_floatingSUB = (totalSUB-totalBEACHED_sub)/totalSUB*100;