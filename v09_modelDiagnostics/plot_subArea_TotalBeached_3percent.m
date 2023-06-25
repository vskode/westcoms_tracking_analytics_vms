clear
close all

load v07_01_3percent_sub.mat
load v08_01_3percent_sub.mat
load v09_02_3percent_sub.mat
load ../wind_speed/beachLocation_windAV.mat
load ../wind_speed/Av_windspeed_year.mat
load ../wind_speed/westcoms_wind_year2020.mat
load ../wind_speed/av_wind_speed_daily_year.mat
load ../wind_speed/av_wind_speed.mat

% Nsub_x1 = Nsub_x1(:,[300:366]);
% Nsub_x2 = Nsub_x2(:,[300:366]);
% Nsub_x3 = Nsub_x3(:,[300:366]);

% find the direction 
wind_dir_to = 180/pi * atan2(clydeU, clydeV);
wind_dir_to = mod(wind_dir_to,360);

wind_dir_from = wind_dir_to + 180;
wind_dir_from = mod(wind_dir_from,360);
clydeU = clydeU';
clydeV = clydeV';

% convert time 
wdatenum = datenum(timevec_day); 
wdatenum(:,2) = wdatenum(:,1);

% change degree to radian 
wradian = deg2rad(wind_dir_from);
wradian = wradian';

[u,v] = pol2cart(wradian, [clydeU,clydeV]); 

q = quiver(wdatenum, wdatenum.*0,u,v,0,'b');
q.ShowArrowHead='off';   
xlabel('Time (days)','fontsize',15)  
ylabel('Wind Speed (ms^-1)','fontsize',15);

M_dailybeached_3percent_SUB1(365) = [];
M_dailybeached_3percent_SUB2(365) = [];
av_wind_speed(366) = [];
av_wind_speed(365) = []; 

Nsub_x1 = Nsub_x1*1104;
Nsub_x2 = Nsub_x2*1104;
Nsub_x3 = Nsub_x3*1104;
% Nsub_x4 = Nsub_x4*1104;

M_dailybeached_3percent_SUB1 = M_dailybeached_3percent_SUB1*1104;
M_dailybeached_3percent_SUB2 = M_dailybeached_3percent_SUB2*1104;
M_dailybeached_3percent_SUB3 = M_dailybeached_3percent_SUB3*1104;

M_dailybeached_3percent_SUB4 = M_dailybeached_3percent_SUB1(300:364,:);
M_dailybeached_3percent_SUB5 = M_dailybeached_3percent_SUB2(300:364,:);
M_dailybeached_3percent_SUB6 = M_dailybeached_3percent_SUB3(300:364,:);
av_wind_speed_65 = av_wind_speed(:,300:364);

u1 = 1; 
v1 = 0;

figure(1) 
subplot(3,1,1)
q = quiver(wdatenum, wdatenum.*0,u,v,0,'k');
q.ShowArrowHead='off';   
xlim([min(wdatenum(:,1)), max(wdatenum(:,1))]);
set(gca,'xtick', []);
ax = gca;
ax.FontSize = 14;
% xlabel('Time (days)','fontsize',14)  
ylabel('Wind Speed (ms^-^1)','fontsize',14);

subplot(3,1,2)
plot(Nsub_x1, 'g-', 'LineWidth', 1.5)
hold on 
plot(Nsub_x2, 'm-', 'LineWidth', 1.5)
hold on 
plot(Nsub_x3, 'b-', 'LineWidth', 1.5)
% hold on 
% plot(Nsub_x4, 'r-', 'LineWidth', 1.5)
% xlabel('Time (days)');
ylabel('Total Beached Litter / Day', 'fontsize', 14);
xlim([0 365]);
ax = gca;
ax.FontSize = 14;
set(gca,'xtick', []);
% title('Resuspending - Classified 2020');
lgd = legend;
legend({'Sticky Coast','Resuspending Unclassified', 'Resuspending Classified'},'Location','northwest','Orientation','vertical')
lgd.FontSize = 14;

subplot(3,1,3)
plot(M_dailybeached_3percent_SUB1, 'g-', 'LineWidth', 1.5)
hold on 
plot(M_dailybeached_3percent_SUB2, 'm-', 'LineWidth', 1.5)
hold on 
plot(M_dailybeached_3percent_SUB3, 'b-', 'LineWidth', 1.5)
xlabel('Time (days)');
ylabel('Daily Beaching Litter Items', 'FontSize', 14);
xlim([0 365]);
% xlim([300 365]);
ylim([-2*10^6 3*10^6]);
hold on
yyaxis 'right'
plot(av_wind_speed, 'k--', 'LineWidth', 0.5)
ylabel('Vector Averaged Daily Wind Speed (ms^-^1)', 'fontsize', 14, 'color', 'k');
ax = gca;
ax.YColor = 'k';
ax.FontSize = 14;

% subplot(3,1,3)
% plot(M_dailybeached_3percent_SUB6, 'b-')
% xlim([0 65])
% xlabel('Time (days)');
% ylabel('Number of Daily Beaching Litter Items / Day)');
% yyaxis 'right'
% plot(av_wind_speed_65, 'k--', 'LineWidth', 0.5)
% ylabel('Vector Averaged Daily Wind Speed (ms^-^1)', 'fontsize', 12, 'color', 'k');
% ax=gca;
% ax.YColor = 'k';
% xlim([300 365]);
% title('Resuspending - Classified 2020');
% lgd = legend;
% legend({'Sticky Coast','Resuspending Unclassified', 'Resuspending Classified'},'Location','northwest','Orientation','vertical')
% lgd.FontSize = 14;
set(gcf,'units','points','position',[10,10,1100,800])
fig_name = ('subarea_TOTbeached_3percent');
print('-f1', '-dpng', '-loose', '-r500',...
    ['plot_' fig_name '.png']); 

figure(2) 
plot(M_dailybeached_3percent_SUB1, 'g-', 'LineWidth', 1.5)
hold on 
plot(M_dailybeached_3percent_SUB2, 'm-', 'LineWidth', 1.5)
hold on 
plot(M_dailybeached_3percent_SUB3, 'b-', 'LineWidth', 1.5)
xlabel('Time (days)');
ylabel('Number of Daily Beaching Litter Items / Day');
xlim([300 365]);
% title('Resuspending - Classified 2020');
hold on
yyaxis 'right'
plot(av_wind_speed, 'k--', 'LineWidth', 0.5)
ylabel('Vector Averaged Daily Wind Speed (ms^-^1)', 'fontsize', 14, 'color', 'k');
ax = gca;
ax.YColor = 'k';
ax.FontSize = 14;
lgd = legend;
legend({'Sticky Coast','Resuspending Unclassified', 'Resuspending Classified'},'Location','northwest','Orientation','vertical')
lgd.FontSize = 14;
fig_name = ('subarea_Dailybeached_3percent');
print('-f2', '-dpng', '-loose', '-r500',...
    ['plot_' fig_name '.png']); 

figure(3)
plot(Nsub_x1, 'g-', 'LineWidth', 1.5)
hold on 
plot(Nsub_x2, 'm-', 'LineWidth', 1.5)
hold on 
plot(Nsub_x3, 'b-', 'LineWidth', 1.5)
% hold on 
% plot(Nsub_x4, 'r-', 'LineWidth', 1.5)
% xlabel('Time (days)');
ylabel('Total Beached Litter / Day', 'fontsize', 14);
xlim([0 365]);
ax = gca;
ax.FontSize = 14;
set(gca,'xtick', []);
% title('Resuspending - Classified 2020');
lgd = legend;
legend({'Sticky Coast','Resuspending Unclassified', 'Resuspending Classified', 'Resuspending Classified - 5 min'},'Location','northwest','Orientation','vertical')
lgd.FontSize = 14;
fig_name = ('subarea_TOTbeached_3percent_5min');
print('-f3', '-dpng', '-loose', '-r500',...
    ['plot_' fig_name '.png']); 

figure(4)
plot(M_dailybeached_3percent_SUB1, 'g-', 'LineWidth', 1.5)
hold on 
plot(M_dailybeached_3percent_SUB2, 'm-', 'LineWidth', 1.5)
hold on 
plot(M_dailybeached_3percent_SUB3, 'b-', 'LineWidth', 1.5)
xlabel('Time (days)');
ylabel('Daily Beaching Litter Items', 'FontSize', 14);
xlim([300 360]);
hold on
yyaxis 'right'
plot(av_wind_speed, 'k--', 'LineWidth', 0.5)
ylabel('Vector Averaged Daily Wind Speed (ms^-^1)', 'fontsize', 14, 'color', 'k');
ax = gca;
ax.YColor = 'k';
ax.FontSize = 14;

figure(5) 
subplot(3,1,1)
plot(Nsub_x1, 'g-', 'LineWidth', 1.5)
hold on 
plot(Nsub_x2, 'm-', 'LineWidth', 1.5)
hold on 
plot(Nsub_x3, 'b-', 'LineWidth', 1.5)
% hold on 
% plot(Nsub_x4, 'r-', 'LineWidth', 1.5)
% xlabel('Time (days)');
ylabel('Total Beached Litter / Day', 'fontsize', 14);
xlim([0 365]);
ax = gca;
ax.FontSize = 14;
set(gca,'xtick', []);
% title('Resuspending - Classified 2020');
lgd = legend;
legend({'Sticky Coast','Resuspending Unclassified', 'Resuspending Classified'},'Location','northwest','Orientation','vertical')
lgd.FontSize = 14;

subplot(3,1,2)
q = quiver(wdatenum, wdatenum.*0,u,v,0,'k');
q.ShowArrowHead='off';   
xlim([min(wdatenum(:,1)), max(wdatenum(:,1))]);
% xlim([300 365]);
set(gca,'xtick', []);
ax = gca;
ax.FontSize = 14;
% xlabel('Time (days)','fontsize',14)  
ylabel('Wind Speed (ms^-^1)','fontsize',14);

subplot(3,1,3)
plot(M_dailybeached_3percent_SUB1, 'g-', 'LineWidth', 1.5)
hold on 
plot(M_dailybeached_3percent_SUB2, 'm-', 'LineWidth', 1.5)
hold on 
plot(M_dailybeached_3percent_SUB3, 'b-', 'LineWidth', 1.5)
xlabel('Time (days)');
ylabel('Daily Beaching Litter Items', 'FontSize', 14);
% xlim([0 365]);
xlim([300 365]);
ylim([-2*10^6 3*10^6]);
hold on
yyaxis 'right'
plot(av_wind_speed, 'k--', 'LineWidth', 0.5)
ylabel('Vector Averaged Daily Wind Speed (ms^-^1)', 'fontsize', 14, 'color', 'k');
ax = gca;
ax.YColor = 'k';
ax.FontSize = 14;

% subplot(3,1,3)
% plot(M_dailybeached_3percent_SUB6, 'b-')
% xlim([0 65])
% xlabel('Time (days)');
% ylabel('Number of Daily Beaching Litter Items / Day)');
% yyaxis 'right'
% plot(av_wind_speed_65, 'k--', 'LineWidth', 0.5)
% ylabel('Vector Averaged Daily Wind Speed (ms^-^1)', 'fontsize', 12, 'color', 'k');
% ax=gca;
% ax.YColor = 'k';
% xlim([300 365]);
% title('Resuspending - Classified 2020');
% lgd = legend;
% legend({'Sticky Coast','Resuspending Unclassified', 'Resuspending Classified'},'Location','northwest','Orientation','vertical')
% lgd.FontSize = 14;
set(gcf,'units','points','position',[10,10,1100,800])
fig_name = ('subarea_TOTbeached_3percent_60days');
print('-f5', '-dpng', '-loose', '-r500',...
    ['plot_' fig_name '.png']); 
