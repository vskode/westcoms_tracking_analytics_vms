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

%load model output 
% load ../cluster_output/01_ModelRuns_ClydePaper/V09_01_3percent/v09_01_3percent_20200616.mat  
% load ../cluster_output/01_ModelRuns_ClydePaper/V09_02_3percent/v09_02_3percent_20201230.mat  

%load polygons 
polygons_seaRegion

westcoast = load('E:\2020_particle_track_models\2021_WestComs_na\westcoms_data\ukireland_coastline.NaN.dat');
l = find(isnan(westcoast(:,1)));

lonlim = [-10 -4];
latlim = [54 59.2]; 

% iactive = find(active);   % a list of the active particles
% ibeached = find(ashore); % a list of beached particles 
% nactive = numel(iactive); % the number of active particles
% nashore = numel(iashore); % the number of beached particles

%% find particles in each polygon
% Clyde Sea 
clydex = Clyde.Vertices(:,1); 
clydey = Clyde.Vertices(:,2); 

clydexq = sub_x;
clydeyq = sub_y; 

[inClyde] = inpolygon(clydexq,clydeyq,clydex,clydey);

outclyde = numel(clydexq(~inClyde));
clyde = totalACTIVE-outclyde;

% Firth of Lorn 
lornx = Lorn.Vertices(:,1); 
lorny = Lorn.Vertices(:,2); 

lornxq = sub_x;
lornyq = sub_y; 

[inLorn] = inpolygon(lornxq,lornyq,lornx,lorny);

%Little Minch 
lminchx = LMinch.Vertices(:,1); 
lminchy = LMinch.Vertices(:,2); 

lminchxq = sub_x;
lminchyq = sub_y; 

[inLMinch] = inpolygon(lminchxq,lminchyq,lminchx,lminchy);

%Upper Minch 
uminchx = UMinch.Vertices(:,1); 
uminchy = UMinch.Vertices(:,2); 

uminchxq = sub_x;
uminchyq = sub_y; 

[inUMinch] = inpolygon(uminchxq,uminchyq,uminchx,uminchy);

% Outer Hebs 
hebsx = Hebs.Vertices(:,1); 
hebsy = Hebs.Vertices(:,2); 

hebsxq = sub_x;
hebsyq = sub_y; 

[inHebs] = inpolygon(hebsxq,hebsyq,hebsx,hebsy);

% North Channel 
channelx = Channel.Vertices(:,1); 
channely = Channel.Vertices(:,2); 

channelxq = sub_x;
channelyq = sub_y; 

[inChannel] = inpolygon(channelxq,channelyq,channelx,channely);

% % Solway Coast
% solcolx = SolCol.Vertices(:,1); 
% solcoly = SolCol.Vertices(:,2); 
% 
% solcolxq = beachedx;
% solcolyq = beachedy; 
% 
% [inSolCol] = inpolygon(solcolxq,solcolyq,solcolx,solcoly);

%% count how many particles in each polygon 
inclyde = numel(clydexq(inClyde)); 
inlorn = numel(lornxq(inLorn)); 
inlminch = numel(lminchxq(inLMinch)); 
inuminch = numel(uminchxq(inUMinch)); 
inhebs = numel(hebsxq(inHebs)); 
inchannel = numel(channelxq(inChannel));
% insolcol = numel(solcolxq(inSolCol));
inireland = totalBEACHED_sub - inclyde - inlorn - inlminch - inhebs - inchannel; 

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
plot(clydexq(inClyde),clydeyq(inClyde), 'g.')
plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')

set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

hold on

%Lorn
plot(lornx,lorny)
axis equal 

hold on 
plot(lornxq(inLorn),lornyq(inLorn), 'r.')
% plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')

set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

hold on

%LittleMinch
plot(lminchx,lminchy)
axis equal 

hold on 
plot(lminchxq(inLMinch),lminchyq(inLMinch), 'm.')
% plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')

set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

hold on

%UpperMinch
plot(uminchx,uminchy)
axis equal 

hold on 
plot(uminchxq(inUMinch),uminchyq(inUMinch), 'c.')
% plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')

set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

hold on

%OuterHebs
plot(hebsx,hebsy)
axis equal 

hold on 
plot(hebsxq(inHebs),hebsyq(inHebs), 'g.')
% plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')

set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

hold on

%Channel
plot(channelx,channely)
axis equal 

hold on 
plot(channelxq(inChannel),channelyq(inChannel), 'y.')
% plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')

set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)

hold on

%SolwayCoast
% plot(solcolx,solcoly)
% axis equal 
% 
% hold on 
% plot(solcolxq(inSolCol),solcolyq(inSolCol), 'b.')
% plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')

set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
xlabel('Longitude', 'FontSize', 20); 
ylabel('Latitude', 'FontSize', 20); 

%% calculate % beached in each polygon of the total released particles 
percentClyde = inclyde/totalACTIVE*100;
percentUMinch = inuminch/totalACTIVE*100;
percentHebs = inhebs/totalACTIVE*100;
percentLorn = inlorn/totalACTIVE*100; 
percentLMinch = inlminch/totalACTIVE*100;
percentChannel = inchannel/totalACTIVE*100;
percentireland = inireland/totalACTIVE*100;

%% calculate % beached in each polygon of the total beached in the sub-area 
percentClyde_sub = inclyde/totalBEACHED_sub*100;
percentUMinch_sub = inuminch/totalBEACHED_sub*100;
percentHebs_sub = inhebs/totalBEACHED_sub*100;
percentLorn_sub = inlorn/totalBEACHED_sub*100; 
percentLMinch_sub = inlminch/totalBEACHED_sub*100;
percentChannel_sub = inchannel/totalBEACHED_sub*100;
percentireland_sub = inireland/totalBEACHED_sub*100;

%% calculate % beached in each polygon of the total beached particles  
percentClyde_tot = inclyde/totalBEACHED_tot*100;
percentUMinch_tot = inuminch/totalBEACHED_tot*100;
percentHebs_tot = inhebs/totalBEACHED_tot*100;
percentLorn_tot = inlorn/totalBEACHED_tot*100; 
percentLMinch_tot = inlminch/totalBEACHED_tot*100;
percentChannel_tot = inchannel/totalBEACHED_tot*100;
percentireland_tot = inireland/totalBEACHED_tot*100;

% %load polygons 
% polygons_seaRegion
% 
% westcoast = load('../data/ukireland_coastline.NaN.dat');
% l = find(isnan(westcoast(:,1)));
% 
% lonlim = [-10 -4];
% latlim = [54 59.2]; 
% 
% totalACTIVE = nactive + nbeached; 
% 
% %% find particles in each polygon
% % Clyde Sea 
% clydex = Clyde.Vertices(:,1); 
% clydey = Clyde.Vertices(:,2); 
% 
% clydexq = x;
% clydeyq = y; 
% 
% [inClyde] = inpolygon(clydexq,clydeyq,clydex,clydey);
% 
% outclyde = numel(clydexq(~inClyde));
% 
% % Firth of Lorn 
% lornx = Lorn.Vertices(:,1); 
% lorny = Lorn.Vertices(:,2); 
% 
% lornxq = x;
% lornyq = y; 
% 
% [inLorn] = inpolygon(lornxq,lornyq,lornx,lorny);
% 
% %Little Minch 
% lminchx = LMinch.Vertices(:,1); 
% lminchy = LMinch.Vertices(:,2); 
% 
% lminchxq = x;
% lminchyq = y; 
% 
% [inLMinch] = inpolygon(lminchxq,lminchyq,lminchx,lminchy);
% 
% %Upper Minch 
% uminchx = UMinch.Vertices(:,1); 
% uminchy = UMinch.Vertices(:,2); 
% 
% uminchxq = x;
% uminchyq = y; 
% 
% [inUMinch] = inpolygon(uminchxq,uminchyq,uminchx,uminchy);
% 
% % Outer Hebs 
% hebsx = Hebs.Vertices(:,1); 
% hebsy = Hebs.Vertices(:,2); 
% 
% hebsxq = x;
% hebsyq = y; 
% 
% [inHebs] = inpolygon(hebsxq,hebsyq,hebsx,hebsy);
% 
% % North Channel 
% channelx = Channel.Vertices(:,1); 
% channely = Channel.Vertices(:,2); 
% 
% channelxq = x;
% channelyq = y; 
% 
% [inChannel] = inpolygon(channelxq,channelyq,channelx,channely);
% 
% % Solway Coast
% solcolx = SolCol.Vertices(:,1); 
% solcoly = SolCol.Vertices(:,2); 
% 
% solcolxq = x;
% solcolyq = y; 
% 
% [inSolCol] = inpolygon(solcolxq,solcolyq,solcolx,solcoly);
% 
% %% count how many particles in each polygon 
% inclyde = numel(clydexq(inClyde)); 
% inlorn = numel(lornxq(inLorn)); 
% inlminch = numel(lminchxq(inLMinch)); 
% inuminch = numel(uminchxq(inUMinch)); 
% inhebs = numel(hebsxq(inHebs)); 
% inchannel = numel(channelxq(inChannel));
% insolcol = numel(solcolxq(inSolCol));
% 
% %% PLOT FIGURE 
% figure(1)
% 
% clf;
% hold on
% 
% % first draw the coast 
% for i = 1:size(l)-1
%     plot(westcoast(l(i) + 1:l(i+1)-1, 1), westcoast(l(i) + 1:l(i+1)-1, 2), 'k')
% end 
% 
% %Clyde
% plot(clydex,clydey)
% axis equal 
% 
% hold on 
% plot(clydexq(inClyde),clydeyq(inClyde), 'g.')
% plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')
% 
% set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
% set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
% 
% hold on
% 
% %Lorn
% plot(lornx,lorny)
% axis equal 
% 
% hold on 
% plot(lornxq(inLorn),lornyq(inLorn), 'r.')
% % plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')
% 
% set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
% set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
% 
% hold on
% 
% %LittleMinch
% plot(lminchx,lminchy)
% axis equal 
% 
% hold on 
% plot(lminchxq(inLMinch),lminchyq(inLMinch), 'm.')
% % plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')
% 
% set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
% set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
% 
% hold on
% 
% %UpperMinch
% plot(uminchx,uminchy)
% axis equal 
% 
% hold on 
% plot(uminchxq(inUMinch),uminchyq(inUMinch), 'c.')
% % plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')
% 
% set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
% set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
% 
% hold on
% 
% %OuterHebs
% plot(hebsx,hebsy)
% axis equal 
% 
% hold on 
% plot(hebsxq(inHebs),hebsyq(inHebs), 'c.')
% % plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')
% 
% set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
% set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
% 
% hold on
% 
% %Channel
% plot(channelx,channely)
% axis equal 
% 
% hold on 
% plot(channelxq(inChannel),channelyq(inChannel), 'y.')
% % plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')
% 
% set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
% set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
% 
% hold on
% 
% %SolwayCoast
% plot(solcolx,solcoly)
% axis equal 
% 
% hold on 
% plot(solcolxq(inSolCol),solcolyq(inSolCol), 'b.')
% % plot(clydexq(~inClyde),clydeyq(~inClyde), 'b.')
% 
% set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
% set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
% 
% hold on
% 
% %% calculate % in each polygon 
% percentClyde = inclyde/totalACTIVE*100;
% percentUMinch = inuminch/totalACTIVE*100;
% percentHebs = inhebs/totalACTIVE*100;
% percentLorn = inlorn/totalACTIVE*100; 
% percentLMinch = inlminch/totalACTIVE*100;
% percentChannel = inchannel/totalACTIVE*100;
% percentSolCol = insolcol/totalACTIVE*100;