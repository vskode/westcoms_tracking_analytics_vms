%% now try by plotting a circle with a radius round the main point 
% load ../cluster_output/01_ModelRuns_ClydePaper/V08_01_3percent20days/v08_01_output_edited_3percent.mat 
% load ../cluster_output/01_ModelRuns_ClydePaper/V08_01_3percent20days/v08_01_3percent_20days20201230.mat 
% load ../cluster_output/01_ModelRuns_ClydePaper/V07_01_3percent/v07_01_output_edited_3percent.mat 

close all 
clear 

addpath('E:/2020_particle_track_models/2021_WestComs_na');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms/analytic_data');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms/v09_modelDiagnostics');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_analytics/westcoms_velocity');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_data');

% load beached particle positions 
v09_vms_3percent_modelDiagnostics

% load required data
load('E:/2020_particle_track_models/2021_WestComs_na/westcoms_analytics/westcoms_velocity/westcoms_uv_AVERAGE2020.mat');
load westcoms2_mesh.mat
westcoast = load('E:\2020_particle_track_models\2021_WestComs_na\westcoms_data\ukireland_coastline.NaN.dat');
l = find(isnan(westcoast(:,1)));

x = westcoast(:,1);
y = westcoast(:,2);
lonlim = [-7.9 -4.23];
latlim = [55.06 58.64]; 

% lonlim = [-8.15 -4];
% latlim = [54.02 59.18]; 

theta = linspace(0,2*pi,100);
for i = 1:length(sub_x)
cx = sub_x(i,:); % centre points 
cy = sub_y(i,:); 
R = 500; % in metres 
R_km = R/1000; % convert to km 
r = R_km/111.111; % in degrees
X = cx + r*cos(theta);
Y = cy + r*sin(theta); 

% check the measurement of the radius in meters 
[dist, angle] = sw_dist([X(1,1) X(1,50)], [Y(1,1) Y(1,50)], 'km'); % calculates the diameter of the circle 
in_meters = dist*1000; % convert to metres 

%% find particles that lie inside the circle (only beached ones)

% now count how many are inside the circle 
part_count(i,:) = (sum((sub_x-cx).^2+(sub_y-cy).^2<=r^2)); % need to convert to no. of plastic items 
part_idx = (sub_x-cx).^2+(sub_y-cy).^2<=r^2;

end

% standardise to 100m 
part_count = part_count/10; 
h = 0.09; % grid size

%%% create positions for currents 
mesh_x = Mesh.uvnode(:,1);
mesh_y = Mesh.uvnode(:,2);  

u = double(u(:,1)); 
v = double(v(:,1)); 

x1 = -3;
y1 = 60;
u1=1;
v1=0; 

figure(1)
clf;
hold on

% first draw the coast 
for i = 1:size(l)-1
    plot(westcoast(l(i) + 1:l(i+1)-1, 1), westcoast(l(i) + 1:l(i+1)-1, 2), 'k')
end 

%% use scatplot

dens = scatter(sub_x, sub_y, [], part_count, 'fill');
hold on 
% cquiver([x1;mesh_x],[y1;mesh_y],[u1;u],[v1;v], 'sampling', h, 'scale', 0.55, 'HeadScale', 0.55);
ax = gca;
ax.FontSize = 20;
dens.SizeData = 16;
colormap jet;
a = colorbar;
set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
set(gcf, 'position', [53 53 650 800]);
view(0,90); 
set(gcf, 'PaperPositionMode', 'auto')
% caxis([0 10]);
caxis([0.05 0.5]);
% caxis([0 5*10^2]);
title('Fishing Release', 'FontSize', 20); 
% ylabel('Latitude', 'FontSize', 20); 
ylabel(a, 'Beached Model Particles / 100m' ,'FontSize', 20, 'Rotation',90);
fig_name = ('beachDens_3percent_2020');
print('-f1', '-dpng', '-loose', '-r500',...
    ['v09_plot_vms_' fig_name '.png']);

% figure(2)
% clf;
% hold on
% % first draw the coast 
% for i = 1:size(l)-1
%     plot(westcoast(l(i) + 1:l(i+1)-1, 1), westcoast(l(i) + 1:l(i+1)-1, 2), 'k')
% end 
% map = colormap(jet);
% colorbar;
% ind = fix((part_count-min(part_count))/(max(part_count)-min(part_count))*(size(map,1)-1))+1;
% h = [];
% %much more efficient than matlab's scatter plot
% for k=1:size(map,1) 
%     if any(ind==k)
%         h(end+1) = line('Xdata',beachedx(ind==k),'Ydata',beachedy(ind==k), ...
%             'LineStyle','none','Color',map(k,:), ...
%             'Marker','.','MarkerSize',7);
%     end
% end
% set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
% set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
% set(gcf, 'position', [53 53 480 850]);
% view(0,90); 
% set(gcf, 'PaperPositionMode', 'auto')
% caxis([0 1*10^3]);
% title('1% wind drift', 'FontSize', 20); 
% % ylabel('Latitude', 'FontSize', 20); 
% fig_name = ('ModelDens_1percent20day_2020');
% print('-f2', '-dpng', '-loose', '-r500',...
%     ['v08_plot_' fig_name '.png']);

%%%%%%% FLOATING PARTICLES 

for j = 1:length(totalx_active)
tot_cx = totalx_active(j,:); % centre points 
tot_cy = totaly_active(j,:); 
tot_R = 5000; % in metres 
tot_R_km = tot_R/1000; % convert to km 
tot_r = tot_R_km/111.111; % in degrees
tot_X = tot_cx + tot_r*cos(theta);
tot_Y = tot_cy + tot_r*sin(theta); 

% check the measurement of the radius in meters 
[tot_dist, tot_angle] = sw_dist([tot_X(1,1) tot_X(1,50)], [tot_Y(1,1) tot_Y(1,50)], 'km'); % calculates the diameter of the circle 
tot_in_meters = tot_dist*1000; % convert to metres 

%% find particles that lie inside the circle (only beached ones)

% now count how many are inside the circle 
tot_part_count(j,:) = (sum((totalx_active-tot_cx).^2+(totaly_active-tot_cy).^2<=tot_r^2)); % need to convert to no. of plastic items 
tot_part_idx = (totalx_active-tot_cx).^2+(totaly_active-tot_cy).^2<=tot_r^2;

end

% standardise to 100m 
tot_part_count = tot_part_count/100; 

figure(2)
clf;
hold on

% first draw the coast 
for i = 1:size(l)-1
    plot(westcoast(l(i) + 1:l(i+1)-1, 1), westcoast(l(i) + 1:l(i+1)-1, 2), 'k')
end 

%% use scatplot

float_dens = scatter(totalx_active, totaly_active, [], tot_part_count, 'fill');
hold on 
ax = gca;
ax.FontSize = 20;
dens.SizeData = 16;
colormap jet;
b = colorbar;
set(gca,'ydir','normal','dataaspectratio',[1 cosd(mean(ylim)) 1])
set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1],'xlim',lonlim,'ylim',latlim)
set(gcf, 'position', [53 53 650 800]);
view(0,90); 
set(gcf, 'PaperPositionMode', 'auto')
caxis([0.01 0.05]);
title('3% wind drift', 'FontSize', 20);
ylabel(b, 'Floating Model Particles / 100m' ,'FontSize', 20, 'Rotation',90);
% ylabel('Latitude', 'FontSize', 20);
fig_name = ('floatDens_3percent_2020');
print('-f2', '-dpng', '-loose', '-r500',...
    ['v09_plot_vms_' fig_name '.png']);