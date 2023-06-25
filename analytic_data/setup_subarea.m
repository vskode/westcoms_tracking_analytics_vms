%% create a script to select a sub-area from westcoms2 mesh 
%% - include Clyde Sea and Islands - Mull, Tiree, Islay 
%% - avoid the open boundaries 

clear 
close all 

% first load in full model output 
load ../'cluster_output/01_ModelRuns_ClydePaper/v09_03_3percent/v09_03_3percent_20210913'.mat  
addpath('E:/2020_particle_track_models/2021_WestComs_na/cluster_output/01_ModelRuns_ClydePaper/v09_03_3percent');

v09_03_beachedPartsOnly_subarea_3percent

westcoast = load('../data/ukireland_coastline.NaN.dat');
l = find(isnan(westcoast(:,1)));

lonlim = [-10 -4];
latlim = [54 59.2]; 

sub_area = polyshape([-7.1 -4.3665 -4.3665 -7.1], [55.06 55.06 57.77 57.77]);
% sub_area = polyshape([-7.9 -4.36 -4.36 -7.9], [54.1 54.1 58.64 58.64]);

% Clyde Sea 
sub_areax = sub_area.Vertices(:,1); 
sub_areay = sub_area.Vertices(:,2); 

sub_areaxq = beachedx;
sub_areayq = beachedy; 

[insub_area] = inpolygon(sub_areaxq,sub_areayq,sub_areax,sub_areay);

%% PLOT FIGURE 
figure(1)

clf;
hold on

% first draw the coast 
for i = 1:size(l)-1
    plot(westcoast(l(i) + 1:l(i+1)-1, 1), westcoast(l(i) + 1:l(i+1)-1, 2), 'k')
end 

%Clyde
plot(sub_area)
hold on 
% plot(sub_areaxq(insub_area),sub_areayq(insub_area), 'r.')
set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
