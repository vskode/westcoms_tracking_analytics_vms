clear
close all 

addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms/v09_modelDiagnostics');
addpath('E:/2020_particle_track_models/2021_WestComs_na/westcoms_tracking_analytics_vms/v09_windDirection');
% v08_plastic_conversion3percent 
load beachLocation_windAV.mat
% load ../wind_speed/av_wind_speed_daily_year.mat
load v09_vms_3percent_modelDiagnostics.mat

av_wind_speed = sqrt(clydeU.^2 + clydeV.^2); % sqrt(u.^2 + v.^2) 
wind_dir_to = 180/pi * atan2(clydeU, clydeV);
wind_dir_to = mod(wind_dir_to,360);

wind_dir_from = wind_dir_to + 180;
wind_dir_from = mod(wind_dir_from,360);

% convert wind_dir_from from degrees to radians 
wind_dir_from_R = deg2rad(wind_dir_from); 
% wind_dir_from_R(366) = [];
% av_wind_speed(366) = [];
wind_dir_from_R = wind_dir_from_R(301:365);
av_wind_speed = av_wind_speed(301:365);

% M_dailybeached_3percent(366) = [];

% beach_percent = v08_L_daily_beachPercent_3percent;
% 
% beach_percent(366) = []; 

% restore daily beached variable 
dailybeached_3percent_vms = M_dailybeached_3percent(301:365); 

% dailybeached_3percent(367) = [];

% convert to plastic items 
% dailybeached_3percent_vms = dailybeached_3percent_vms*1104;

figure(1)
plot_one = polarscatter(wind_dir_from_R, av_wind_speed, 200, dailybeached_3percent_vms, 'filled', 'o'); 
hold on
datarange_one = [-2*10^6 2*10^6];
breakpoint = 0;
ctlen = 3;
colors = [0 0 1; 0 0.0625 1; 0 0.156 1; 0 0.2500 1; 0 0.3440 1.0000; 0 0.4380 1; 0 0.5310 1; 0 0.6250 1; 0 0.7190 1; 0 0.8130 1; 0 0.9060 1; 0 1 1; 1 1 0; 1 0.9380 0; 1 0.8440 0; 1 0.7500 0; 1 0.6560 0; 1 0.5630 0; 1 0.4690 0; 1 0.3750 0; 1 0.2810 0; 1 0.1880 0; 1 0.0938 0; 1 0 0];
pos = (max(datarange_one)-breakpoint)/diff(datarange_one);
n = [floor(ctlen*(1-pos)) ceil(ctlen*pos)];
cmap = [repmat(colors(1:12,:),[n(1),1]); repmat(colors(13:24,:),[n(1),1])];
colormap(cmap);
a = colorbar;
pax = gca;
pax.FontSize = 14;
pax.LineWidth = 1.5;
caxis([-1*10^3,1*10^3]);
ylabel(a, 'Number of Daily Beaching Plastic Items' ,'FontSize', 14, 'Rotation',90);
hold on
% title('3% wind drift');
fig_name = ('windscatter_3percent_2020_vms');
print('-f1', '-dpng', '-loose', '-r500',...
    ['v09_plot_' fig_name '.png']);

% figure(2)
% plot_two = polarscatter(wind_dir_from_R, av_wind_speed, 75, beach_percent, 'filled', 'o'); 
% hold on
% datarange_two = [-100 100];
% breakpoint = 0;
% ctlen = 3;
% colors = [0 0 1; 0 0.0625 1; 0 0.156 1; 0 0.2500 1; 0 0.3440 1.0000; 0 0.4380 1; 0 0.5310 1; 0 0.6250 1; 0 0.7190 1; 0 0.8130 1; 0 0.9060 1; 0 1 1; 1 1 0; 1 0.9380 0; 1 0.8440 0; 1 0.7500 0; 1 0.6560 0; 1 0.5630 0; 1 0.4690 0; 1 0.3750 0; 1 0.2810 0; 1 0.1880 0; 1 0.0938 0; 1 0 0];
% pos = (max(datarange_two)-breakpoint)/diff(datarange_two);
% n = [floor(ctlen*(1-pos)) ceil(ctlen*pos)];
% cmap = [repmat(colors(1:12,:),[n(1),1]); repmat(colors(13:24,:),[n(1),1])];
% colormap(cmap);
% colorbar;
% pax = gca;
% pax.FontSize = 14;
% caxis([-100 100]);
% hold on
% title('3% wind drift');
% fig_name = ('windscatterPercent_3percent_2020');
% print('-f2', '-dpng', '-loose', '-r500',...
%     ['v08_plot_' fig_name '.png']);